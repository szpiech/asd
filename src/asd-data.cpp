#include "asd-data.h"

const char DEL = ' ';
const string EMPTY_STRING = " ";

double **DIST_MAT;
int **NUM_LOCI;
int **IBS_0_MAT;
int **IBS_1_MAT;
int **IBS_2_MAT;
pthread_mutex_t mutex_dist_mat = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_loci_mat = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ibs_0 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ibs_1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_ibs_2 = PTHREAD_MUTEX_INITIALIZER;


bool readData_Check(igzstream &fin, structure_data &data,
                    int sort, int ndcols, int ndrows,
                    int nrows, int ncols)
{
    bool FILE_STATUS = true;

    string line;
    int nind = nrows / 2;
    data.nind = nind;
    data.data = new short*[nind];
    data.ind_names = new string[nind];

    for (int i = 0; i < nind; i++)
    {
        data.ind_names[i] = EMPTY_STRING;
    }

    getline(fin, line);
    int size;
    data.locus_names = split_str_str(size, line.c_str(), DEL);

    size = ncols;
    data.nloci = size;

    for (int i = 1; i < ndrows; i++)
    {
        getline(fin, line);
    }

    string key;
    string field;
    short int *tmp;
    //int **block;
    //double *block;
    short tmp_dbl;
    int index;


    short *ind_count = new short[nind];

    for (int i = 0; i < nind; i++) ind_count[i] = 0;

    for (int row = 0; row < nrows; row++)
    {
        tmp = NULL;
        //block = NULL;
        for (int i = 1; i <= ndcols; i++)
        {
            fin >> field;
            if (i == sort) key = field;
        }

        index = search(data.ind_names, nind, key);
        tmp = split_int(fin, ncols);

        if (index >= 0)
        {
            ind_count[index]++;

            if (ind_count[index] > 2)
            {
                FILE_STATUS &= false;
                cerr << "Individual " << key << " seen "
                     << ind_count[index] << "times.\n";
            }

            for (int i = 0; i < size; i++)
            {

                if (tmp[i] != 1 && tmp[i] != 0 && tmp[i] != -9)
                {
                    FILE_STATUS &= false;
                    cerr << "LINE " << row + ndrows << " COL "
                         << i + ndcols << ": Allele '" << tmp[i]
                         << "' is not 0/1/-9.\n";
                }
            }
        }
        else
        {
            index = put(data.ind_names, nind, key);
            if (index > nind - 1)
            {
                cerr << "Found more than " << nind << " individuals.";
                exit(-1);
            }
            ind_count[index]++;

            for (int i = 0; i < size; i++)
            {
                if (tmp[i] != 1 && tmp[i] != 0 && tmp[i] != -9)
                {
                    FILE_STATUS &= false;
                    cerr << "LINE " << row + ndrows << " COL "
                         << i << ": Allele " << tmp[i] << " is not 0/1/-9.\n";
                }
            }
        }
        delete [] tmp;
    }

    return FILE_STATUS;
}

bool checkFile(param_t *params)
{

    string outfile = params->getStringFlag(ARG_OUTFILE);
    string filename = params->getStringFlag(ARG_FILENAME);
    int nrows = params->getIntFlag(ARG_NROWS);
    int ndrows = params->getIntFlag(ARG_NDROWS);
    int ncols = params->getIntFlag(ARG_NCOLS);
    int ndcols = params->getIntFlag(ARG_NDCOLS);
    int sort = params->getIntFlag(ARG_SORT);
    int num_threads = params->getIntFlag(ARG_THREAD);
    int nind = nrows / 2;
    bool PRINT_FULL = params->getBoolFlag(ARG_FULL);
    bool PRINT_FULL_LOG = params->getBoolFlag(ARG_FULL_LOG);
    bool quit = false;
    bool CALC_ALL_IBS = params->getBoolFlag(ARG_CALC_IBS);
    bool CHECK_FILE = params->getBoolFlag(ARG_CHECK_FILE);
    bool CHECK_FILE_DEEP = params->getBoolFlag(ARG_CHECK_FILE_DEEP);
    bool FILE_STATUS = true;

    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "Could not open " << filename << " for reading.'n";
        return 0;
    }

    string junk;
    int fields;
    int counter = 1;
    int ndr_counter = ndrows;
    int obs_cols;
    int obs_rows = 0;
    do
    {
        getline(fin, junk);
        fields = countFields(junk);
        if (fields == 0 && counter < nrows + ndrows)
        {
            FILE_STATUS &= false;
            cerr << "LINE " << counter << ": "
                 << "Blank line found.\n";
        }
        else if (fields == 0)
        {
        }
        else if (ndr_counter > 0)
        {
            if (counter == 1)
            {
                obs_cols = fields;
                if (obs_cols != ncols)
                {
                    FILE_STATUS &= false;
                    cerr << "LINE " << counter << ": ";
                    cerr << "Found " << obs_cols << " loci. Expected "
                         << ncols << ".\n";;
                }
            }
            else
            {
                if (fields != ncols)
                {
                    FILE_STATUS &= false;
                    cerr << "LINE " << counter << ": ";
                    cerr << "Found " << fields << " fields in an assumed "
                         << "header row. Expected " << ncols << ".\n";;
                }
            }
            ndr_counter--;
        }
        else
        {
            if (fields - ncols != ndcols)
            {
                FILE_STATUS &= false;
                cerr << "LINE " << counter << ": ";
                cerr << "Found " << fields << " fields. Expected "
                     << ndcols << " headers + "
                     << ncols << " data = "
                     << ndcols + ncols << " fields.\n";
            }
            obs_rows++;
        }

        counter++;
    }
    while (fin.good());

    if (obs_rows != nrows)
    {
        FILE_STATUS &= false;
        cerr << "Found " << obs_rows << " lines of data. "
             << "Expected " << nrows << ".\n";
    }

    fin.close();

    short **data;

    if (FILE_STATUS && CHECK_FILE_DEEP)
    {
        fin.open(filename.c_str());
        structure_data data;
        FILE_STATUS &= readData_Check(fin, data, sort, ndcols, ndrows, nrows, ncols);
        fin.close();
    }


    return FILE_STATUS;
}


int countFields(string str)
{
    string::iterator it;
    bool inword = false;
    bool space = false;
    int counter = 0;
    for ( it = str.begin() ; it < str.end(); it++ )
    {
        space = isspace(*it);
        if (!inword && !space)
        {
            inword = true;
            counter++;
        }
        else if (inword && space)
        {
            inword = false;
        }
    }

    return counter;
}


//Yes, this is really dumb for a lot of reasons.
void output(void *order)
{
    output_order_t *p = (output_order_t *)order;
    ostream *out = p->out;
    bool PRINT_FULL = p->PRINT_FULL;
    bool PRINT_FULL_LOG = p->PRINT_FULL_LOG;
    string *ind_names = p->ind_names;
    int ncols = p->ncols;
    string type = p->type;
    int nind = p->nind;

    if (type.compare("dist") == 0)
    {
        if (!PRINT_FULL && !PRINT_FULL_LOG)
        {
            for (int i = 0; i < nind; i++)
            {
                for (int j = 0; j < nind; j++)
                {
                    *out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
                }
                *out << endl;
            }
            *out << endl;
        }

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
        }
        *out << endl;

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
            for (int j = 0; j < nind ; j++)
            {
                if (PRINT_FULL_LOG)
                {
                    *out << 0 - log(double(DIST_MAT[i][j]) /
                                    (double(ncols) + double(NUM_LOCI[i][j])))
                         << " ";
                }
                else if (PRINT_FULL)
                {
                    *out << 1 - (double(DIST_MAT[i][j]) /
                                 (double(ncols) + double(NUM_LOCI[i][j])))
                         << " ";
                }
                else *out << DIST_MAT[i][j] << " ";
            }
            *out << endl;
        }
    }
    else if (type.compare("ibs0") == 0)
    {
        if (!PRINT_FULL && !PRINT_FULL_LOG)
        {
            for (int i = 0; i < nind; i++)
            {
                for (int j = 0; j < nind; j++)
                {
                    *out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
                }
                *out << endl;
            }
            *out << endl;
        }

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
        }
        *out << endl;

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
            for (int j = 0; j < nind ; j++)
            {
                if (PRINT_FULL_LOG || PRINT_FULL)
                {
                    *out << double(IBS_0_MAT[i][j]) /
                         (double(ncols) + double(NUM_LOCI[i][j])) << " ";
                }
                else *out << IBS_0_MAT[i][j] << " ";

            }
            *out << endl;
        }
    }
    else if (type.compare("ibs1") == 0)
    {
        if (!PRINT_FULL && !PRINT_FULL_LOG)
        {
            for (int i = 0; i < nind; i++)
            {
                for (int j = 0; j < nind; j++)
                {
                    *out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
                }
                *out << endl;
            }
            *out << endl;
        }

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
        }
        *out << endl;

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
            for (int j = 0; j < nind ; j++)
            {
                if (PRINT_FULL_LOG || PRINT_FULL)
                {
                    *out << double(IBS_1_MAT[i][j]) /
                         (double(ncols) + double(NUM_LOCI[i][j])) << " ";
                }
                else *out << IBS_1_MAT[i][j] << " ";

            }
            *out << endl;
        }
    }
    else if (type.compare("ibs2") == 0)
    {
        if (!PRINT_FULL && !PRINT_FULL_LOG)
        {
            for (int i = 0; i < nind; i++)
            {
                for (int j = 0; j < nind; j++)
                {
                    *out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
                }
                *out << endl;
            }
            *out << endl;
        }

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
        }
        *out << endl;

        for (int i = 0; i < nind; i++)
        {
            *out << ind_names[i] << " ";
            for (int j = 0; j < nind ; j++)
            {
                if (PRINT_FULL_LOG || PRINT_FULL)
                {
                    *out << double(IBS_2_MAT[i][j]) /
                         (double(ncols) + double(NUM_LOCI[i][j])) << " ";
                }
                else *out << IBS_2_MAT[i][j] << " ";

            }
            *out << endl;
        }
    }
    else
    {
        cerr << "UNKNOWN TYPE: " << type << endl;
        exit(-1);
    }
    return;
}

void readData_ind_asd(igzstream &fin, structure_data &data,
                      int sort, int ndcols, int ndrows,
                      int nrows, int ncols, int STRU_MISSING)
{
    string line;
    int nind = nrows / 2;
    data.nind = nind;
    data.data = new short*[nind];
    data.ind_names = new string[nind];

    for (int i = 0; i < nind; i++)
    {
        data.ind_names[i] = EMPTY_STRING;
    }
    int size;

    for (int i = 0; i < ndrows; i++)
    {
        getline(fin, line);
        //this would load loci names, but we don't actualy use them, so why bother.
        //if(i==0) data.locus_names = split_str_str(size,line.c_str(),DEL);
    }

    size = ncols;
    data.nloci = size;

    string key;
    string field;
    short int *tmp;
    //int **block;
    //double *block;
    short tmp_dbl;
    int index;
    for (int row = 0; row < nrows; row++)
    {
        tmp = NULL;
        //block = NULL;
        for (int i = 1; i <= ndcols; i++)
        {
            fin >> field;
            if (i == sort) key = field;
        }

        index = search(data.ind_names, nind, key);

        //getline(fin,line);
        //cerr << line << endl;
        tmp = split_int(fin, ncols);

        if (index >= 0)
        {
            for (int i = 0; i < size; i++)
            {
                tmp_dbl = data.data[index][i];
                //(*data.data)[key][1][i] = tmp[i];
                if (tmp_dbl != STRU_MISSING && tmp[i] != STRU_MISSING)
                {
                    data.data[index][i] += tmp[i];
                }
                else
                {
                    data.data[index][i] = -9;
                }
            }
        }
        else
        {
            index = put(data.ind_names, nind, key);
            data.data[index] = new short[size];
            //block = new double[size];
            /*
            block = new int[2];
            block[0] = new int[size];
            block[1] = new int[size];
            */
            for (int i = 0; i < size; i++)
            {
                if (tmp[i] != STRU_MISSING) data.data[index][i] = tmp[i];
                else data.data[index][i] = -9;
            }
            //data.data[index] = block;
        }
        delete [] tmp;
    }

    return;
}

void readData_ind_asd_tped_tfam(string tped_filename, string tfam_filename, structure_data &data,
                                int &nrow, int &nloci, string TPED_MISSING)
{
    string junk;

    igzstream famin, pedin;

    pedin.open(tped_filename.c_str());
    if (pedin.fail())
    {
        cerr << "Could not open " << tped_filename << " for reading1.\n";
        throw - 1;
    }

    famin.open(tfam_filename.c_str());
    if (famin.fail())
    {
        cerr << "Could not open " << tfam_filename << " for reading1.\n";
        throw - 1;
    }

    //int start = famin.tellg();
    while (getline(famin, junk)) nrow += 2;
    famin.close();
    famin.clear();
    //famin.seekg(start);

    //start = pedin.tellg();
    while (getline(pedin, junk)) nloci++;
    pedin.close();
    pedin.clear();
    //pedin.seekg(start);

    pedin.open(tped_filename.c_str());
    if (pedin.fail())
    {
        cerr << "Could not open " << tped_filename << " for reading2.\n";
        throw - 1;
    }

    famin.open(tfam_filename.c_str());
    if (famin.fail())
    {
        cerr << "Could not open " << tfam_filename << " for reading2.\n";
        throw - 1;
    }

    int nind = nrow / 2;

    cerr << "Reading " << nind << " diploid individuals at " << nloci << " loci.\n";

    data.nind = nind;
    data.data = new short*[nind];
    for (int i = 0; i < nind; i++) data.data[i] = new short[nloci];
    data.ind_names = new string[nind];

    for (int i = 0; i < nind; i++)
    {
        famin >> junk;
        famin >> data.ind_names[i];
        getline(famin, junk);
    }

    data.nloci = nloci;
    //data.locus_names = new string[nloci];
    /*
    string *zeroAllele = new string[nloci];
    for (int i = 0; i < nloci; i++) zeroAllele[i] = TPED_MISSING;
    */
    string zeroAllele;
    string allele1, allele2;
    short alleleCount = 0;
    for (int locus = 0; locus < nloci; locus++)
    {
        zeroAllele = TPED_MISSING;
        pedin >> junk;
        //pedin >> data.locus_names[locus];
        pedin >> junk;
        pedin >> junk;
        pedin >> junk;

        for (int ind = 0; ind < nind; ind++)
        {
            alleleCount = 0;
            pedin >> allele1;
            pedin >> allele2;
            if (allele1.compare(TPED_MISSING) == 0 || allele2.compare(TPED_MISSING) == 0)
            {
                data.data[ind][locus] = -9;
            }
            else if (zeroAllele.compare(TPED_MISSING) == 0)
            {
                zeroAllele = allele1;
                if (allele2.compare(zeroAllele) != 0) alleleCount++;
                data.data[ind][locus] = alleleCount;
            }
            else
            {
                if (allele1.compare(zeroAllele) != 0) alleleCount++;
                if (allele2.compare(zeroAllele) != 0) alleleCount++;
                data.data[ind][locus] = alleleCount;
            }
        }
    }

    //delete [] zeroAllele;

    return;
}


void readData_pop_freq(igzstream &fin, structure_data &data,
                       int sort, int ndcols, int ndrows,
                       int nrows, int ncols)
{
    string line;
    int nind = nrows / 2;
    data.nind = nind;
    data.data = new short*[nind];
    data.ind_names = new string[nind];

    for (int i = 0; i < nind; i++)
    {
        data.ind_names[i] = EMPTY_STRING;
    }

    getline(fin, line);
    int size;
    data.locus_names = split_str_str(size, line.c_str(), DEL);

    size = ncols;
    data.nloci = size;

    for (int i = 1; i < ndrows; i++)
    {
        getline(fin, line);
    }

    string key;
    string field;
    short int *tmp;
    //int **block;
    //double *block;
    short tmp_dbl;
    int index;
    for (int row = 0; row < nrows; row++)
    {
        tmp = NULL;
        //block = NULL;
        for (int i = 1; i <= ndcols; i++)
        {
            fin >> field;
            if (i == sort) key = field;
        }

        index = search(data.ind_names, nind, key);

        //getline(fin,line);
        //cerr << line << endl;
        tmp = split_int(fin, ncols);

        if (index >= 0)
        {
            for (int i = 0; i < size; i++)
            {
                tmp_dbl = data.data[index][i];
                //(*data.data)[key][1][i] = tmp[i];
                if (tmp_dbl >= 0 && tmp[i] >= 0)
                {
                    data.data[index][i] += tmp[i];
                }
                else
                {
                    data.data[index][i] = -9;
                }
            }
        }
        else
        {
            index = put(data.ind_names, nind, key);
            data.data[index] = new short[size];
            //block = new double[size];
            /*
            block = new int[2];
            block[0] = new int[size];
            block[1] = new int[size];
            */
            for (int i = 0; i < size; i++)
            {
                data.data[index][i] = tmp[i];
            }
            //data.data[index] = block;
        }
        delete [] tmp;
    }

    return;
}


int search(string *s, int size, string key)
{
    for (int i = 0; i < size; i++)
    {
        if (s[i].compare(key) == 0) return i;
    }
    return -9;
}

int put(string *s, int size, string key)
{
    for (int i = 0; i < size; i ++)
    {
        if (s[i].compare(EMPTY_STRING) == 0)
        {
            s[i] = key;
            return i;
        }
    }
    cerr << "ERROR: This shouldn't be possible.  Try --check-file or --check-deep.\n";
    exit(0);
}

// Split a NULL-terminated character array on a given character into
// a vector of strings
// The vector is passed by reference and cleared each time
// The number of strings split out is returned
short int *split_int(igzstream &fin, int fields)
{
    /*
    vector<string> v;
    while (true)
      {
        const char* begin = s;

        while (*s != c && *s) { ++s; }

        v.push_back(string(begin, s));

        if (!*s)
    {
      break;
    }

        if (!*++s)
    {
      //v.push_back("");
      break;
    }
      }
    */
    short int *a = new short int[fields];
    for (int i = 0; i < fields; i++)
    {
        fin >> a[i];
        //cout << a[i] << endl;
    }
    /*
    for(int i = 0; i < v.size();i++)
      {
        a[i]=atoi(v[i].c_str());
      }
    */
    return a;
}

string *split_str_str(int &size, const char *s, char c)
{
    vector<string> v;
    while (true)
    {
        const char *begin = s;

        while (*s != c && *s)
        {
            ++s;
        }

        v.push_back(string(begin, s));

        if (!*s)
        {
            break;
        }

        if (!*++s)
        {
            v.push_back("");
            break;
        }
    }
    string *x = new string[v.size()];

    for (int i = 0; i < v.size(); i++)
    {
        x[i] = v[i];
    }

    size = v.size();
    return x;
}

