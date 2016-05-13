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

bool init_storage(int nind, bool CALC_ALL_IBS) {
	DIST_MAT = new double*[nind];
	for (int i = 0; i < nind; i++)
	{
		DIST_MAT[i] = new double[nind];
		for (int j = 0; j < nind; j++) DIST_MAT[i][j] = 0;
	}

	NUM_LOCI = new int *[nind];
	for (int i = 0; i < nind; i++)
	{
		NUM_LOCI[i] = new int[nind];
		for (int j = 0; j < nind; j++) NUM_LOCI[i][j] = 0;
	}

	if (CALC_ALL_IBS)
	{
		IBS_0_MAT = new int *[nind];
		for (int i = 0; i < nind; i++)
		{
			IBS_0_MAT[i] = new int[nind];
			for (int j = 0; j < nind; j++) IBS_0_MAT[i][j] = 0;
		}
		IBS_1_MAT = new int *[nind];
		for (int i = 0; i < nind; i++)
		{
			IBS_1_MAT[i] = new int[nind];
			for (int j = 0; j < nind; j++) IBS_1_MAT[i][j] = 0;
		}
		IBS_2_MAT = new int *[nind];
		for (int i = 0; i < nind; i++)
		{
			IBS_2_MAT[i] = new int[nind];
			for (int j = 0; j < nind; j++) IBS_2_MAT[i][j] = 0;
		}
	}

	return true;
}

bool finalize_calculations(int nind, int ncols, bool CALC_ALL_IBS) {
	for (int i = 0; i < nind; i++)
	{
		for (int j = i; j < nind; j++)
		{
			if (i == j)
			{
				DIST_MAT[i][j] = ncols + NUM_LOCI[i][j];
				//NUM_LOCI[i][j] = ncols+NUM_LOCI[i][j];
				if (CALC_ALL_IBS)
				{
					IBS_0_MAT[i][j] = 0;
					IBS_1_MAT[i][j] = 0;
					IBS_2_MAT[i][j] = ncols + NUM_LOCI[i][j];
				}
			}
			else
			{
				DIST_MAT[j][i] = DIST_MAT[i][j];
				NUM_LOCI[j][i] = NUM_LOCI[i][j];
				if (CALC_ALL_IBS)
				{
					IBS_0_MAT[j][i] = IBS_0_MAT[i][j];
					IBS_1_MAT[j][i] = IBS_1_MAT[i][j];
					IBS_2_MAT[j][i] = IBS_2_MAT[i][j];
				}
			}
		}
	}

	return true;
}

void write_ibs_matrices(string outfile, int nind, int ncols, string *ind_names, bool PRINT_PARTIAL) {
	string ibs_fname[3];
	ibs_fname[0] = outfile + ".ibs0";
	ibs_fname[1] = outfile + ".ibs1";
	ibs_fname[2] = outfile + ".ibs2";

	ofstream out;
	for (int ibs = 0; ibs <= 2; ibs++) {

		int **mat;
		if (ibs == 0) {
			mat = IBS_0_MAT;
		}
		else if (ibs == 1) {
			mat = IBS_1_MAT;
		}
		else if (ibs == 2) {
			mat = IBS_2_MAT;
		}

		out.open(ibs_fname[ibs].c_str());
		if (out.fail()) {
			LOG.err("ERROR: Could not open", ibs_fname[ibs]);
			throw 0;
		}
		if (PRINT_PARTIAL) {
			out << nind << endl;
			for (int i = 0; i < nind; i++) {
				for (int j = 0; j < nind; j++) {
					out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
				}
				out << endl;
			}
			out << endl;
		}

		for (int i = 0; i < nind; i++)
		{
			out << ind_names[i] << " ";
		}
		out << endl;

		for (int i = 0; i < nind; i++)
		{
			out << ind_names[i] << " ";
			for (int j = 0; j < nind ; j++)
			{
				if (!PRINT_PARTIAL)
				{
					out << double(mat[i][j]) /
					    (double(ncols) + double(NUM_LOCI[i][j])) << " ";
				}
				else out << mat[i][j] << " ";

			}
			out << endl;
		}
		out.close();
		out.clear();
	}

	return;
}

void write_dist_matrix(string outfile, int nind, int ncols, string *ind_names, bool PRINT_PARTIAL, bool PRINT_LOG) {
	ofstream out;
	if (!PRINT_PARTIAL) {
		outfile += ".asd.dist";
	}
	else {
		outfile += ".asd.partial";
	}

	out.open(outfile.c_str());
	if (out.fail()) {
		LOG.err("ERROR: Could not open", outfile);
		throw 0;
	}
	if (PRINT_PARTIAL)
	{
		out << nind << endl;
		for (int i = 0; i < nind; i++)
		{
			for (int j = 0; j < nind; j++)
			{
				out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
			}
			out << endl;
		}
		out << endl;
	}

	for (int i = 0; i < nind; i++)
	{
		out << ind_names[i] << " ";
	}
	out << endl;

	for (int i = 0; i < nind; i++)
	{
		out << ind_names[i] << " ";
		for (int j = 0; j < nind ; j++)
		{
			if (PRINT_LOG)
			{
				out << 0 - log(double(DIST_MAT[i][j]) /
				               (double(ncols) + double(NUM_LOCI[i][j])))
				    << " ";
			}
			else if (!PRINT_PARTIAL)
			{
				out << 1 - (double(DIST_MAT[i][j]) /
				            (double(ncols) + double(NUM_LOCI[i][j])))
				    << " ";
			}
			else out << DIST_MAT[i][j] << " ";
		}
		out << endl;
	}
	out.close();
	return;
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

structure_data *readData_stru(string infile, int sort, int &nrows, int &ncols, string STRU_MISSING) {
	int ndcols = -1;
	int ndrows = -1;
	nrows = -1;
	ncols = -1;

	string line;
	igzstream fin;

	fin.open(infile.c_str());

	if (fin.fail()) {
		LOG.err("ERROR: Could not open", infile);
		throw - 1;
	}

	int totalRows = 0;
	//int nloci = 0;
	int currentCols = 0;
	while (getline(fin, line)) {
		totalRows++;
		currentCols = countFields(line);
		if (totalRows == 1) ncols = currentCols;

		if (ndrows < 0 && currentCols != ncols) {
			ndcols = currentCols - ncols;
			ndrows = totalRows - 1;
		}
	}
	nrows = totalRows - ndrows;

	bool err = false;
	if (!check_int_gt_0(nrows)) {
		err = true;
		LOG.err("ERROR: Number of chr must be > 0. Found", nrows);
	}

	if (!check_int_gt_0(ndrows)) {
		err = true;
		LOG.err("ERROR: Non-data rows must be > 0. Found", ndrows);
	}
	LOG.log("Non-data header rows:", ndrows);

	if (!check_int_gt_0(ncols)) {
		err = true;
		LOG.err("ERROR: Number of loci must be > 0. Found", ncols);
	}

	if (!check_int_gt_0(ndcols)) {
		err = true;
		LOG.err("ERROR: Non-data columns must be > 0. Found", ndcols);
	}
	LOG.log("Non-data header columns:", ndcols);

	if (!check_sort_ge_ndcols(sort, ndcols)) {
		err = true;
		LOG.err("ERROR: Individual ID column must be in [ 1, ", ndcols, false);
		LOG.err("].");
	}

	if (err) {
		throw 0;
	}

	fin.close();
	fin.clear();

	fin.open(infile.c_str());


	structure_data *data = new structure_data;
	int nind = nrows / 2;
	data->nind = nind;
	data->data = new short*[nind];
	for (int i = 0; i < nind; i++) {
		data->data[i] = new short[ncols];
	}
	data->ind_names = new string[nind];

	//toss out non-data rows
	for (int i = 0; i < ndrows; i++) getline(fin, line);

	data->nloci = ncols;

	map<string, short> *allele2code = new map<string, short>[ncols];
	short *lastAlleleCode = new short[ncols];
	for (int i = 0; i < ncols; i++) {
		allele2code[i][STRU_MISSING] = -9;
		lastAlleleCode[i] = -1;
	}

	int index;
	int indCount = -1;
	string field, key, allele1, allele2, line1, line2;
	stringstream sin1, sin2;
	short genotypeCode;
	for (int ind = 0; ind < nind; ind++) {
		getline(fin, line1);
		sin1 << line1;
		getline(fin, line2);
		sin2 << line2;

		for (int i = 1; i <= ndcols; i++) {
			sin1 >> field;
			sin2 >> field;
			if (i == sort) {
				key = field;
				data->ind_names[ind] = key;
			}
		}

		for (int locus = 0; locus < ncols; locus++) {
			sin1 >> allele1;
			sin2 >> allele2;

			if (allele2code[locus].count(allele1) == 0) {
				lastAlleleCode[locus]++;
				if (lastAlleleCode[locus] > 1) {
					LOG.err("ERORR: --biallelic flag set, but found more than 2 alleles at locus", locus + 1);
					throw - 1;
				}
				allele2code[locus][allele1] = lastAlleleCode[locus];
			}

			if (allele2code[locus].count(allele2) == 0) {
				lastAlleleCode[locus]++;
				if (lastAlleleCode[locus] > 1) {
					LOG.err("ERORR: --biallelic flag set, but found more than 2 alleles at locus", locus + 1);
					throw - 1;
				}
				allele2code[locus][allele2] = lastAlleleCode[locus];
			}

			genotypeCode = allele2code[locus][allele1] + allele2code[locus][allele2];
			if (genotypeCode < 0) {
				data->data[ind][locus] = allele2code[locus][STRU_MISSING];
			}
			else {
				data->data[ind][locus] = genotypeCode;
			}
		}
	}
	fin.close();
	return data;
}

structure_data *readData_stru2(string infile, int sort, int &nrows, int &ncols, string STRU_MISSING) {
	int ndcols = -1;
	int ndrows = -1;
	nrows = -1;
	ncols = -1;

	string line;
	igzstream fin;

	fin.open(infile.c_str());

	if (fin.fail()) {
		LOG.err("ERROR: Could not open", infile);
		throw - 1;
	}

	int totalRows = 0;
	//int nloci = 0;
	int currentCols = 0;
	while (getline(fin, line)) {
		totalRows++;
		currentCols = countFields(line);
		if (totalRows == 1) ncols = currentCols;

		if (ndrows < 0 && currentCols != ncols) {
			ndcols = currentCols - ncols;
			ndrows = totalRows - 1;
		}
	}
	nrows = totalRows - ndrows;

	bool err = false;
	if (!check_int_gt_0(nrows)) {
		err = true;
		LOG.err("ERROR: Number of chr must be > 0. Found", nrows);
	}

	if (!check_int_gt_0(ndrows)) {
		err = true;
		LOG.err("ERROR: Non-data rows must be > 0. Found", ndrows);
	}
	LOG.log("Non-data header rows:", ndrows);

	if (!check_int_gt_0(ncols)) {
		err = true;
		LOG.err("ERROR: Number of loci must be > 0. Found", ncols);
	}

	if (!check_int_gt_0(ndcols)) {
		err = true;
		LOG.err("ERROR: Non-data columns must be > 0. Found", ndcols);
	}
	LOG.log("Non-data header columns:", ndcols);

	if (!check_sort_ge_ndcols(sort, ndcols)) {
		err = true;
		LOG.err("ERROR: Individual ID column must be in [ 1, ", ndcols, false);
		LOG.err("].");
	}

	if (err) {
		throw 0;
	}

	fin.close();
	fin.clear();

	fin.open(infile.c_str());


	structure_data *data = new structure_data;
	int nind = nrows / 2;
	data->nind = nind;
	data->data = new short*[nrows];
	for (int i = 0; i < nrows; i++) {
		data->data[i] = new short[ncols];
	}
	data->ind_names = new string[nind];

	//toss out non-data rows
	for (int i = 0; i < ndrows; i++) getline(fin, line);

	data->nloci = ncols;
	string field, key, allele;
	map<string, short> *allele2code = new map<string, short>[ncols];
	short *lastAlleleCode = new short[ncols];
	for (int i = 0; i < ncols; i++) {
		allele2code[i][STRU_MISSING] = -9;
		lastAlleleCode[i] = -1;
	}
	int index;
	int indCount = -1;
	for (int row = 0; row < nrows; row++) {
		for (int i = 1; i <= ndcols; i++) {
			fin >> field;
			if (i == sort) key = field;
			if (row % 2 == 0) {
				data->ind_names[row / 2] = key;
			}
		}

		for (int locus = 0; locus < ncols; locus++) {
			fin >> allele;
			if (allele2code[locus].count(allele) == 0) {
				lastAlleleCode[locus]++;
				allele2code[locus][allele] = lastAlleleCode[locus];
				data->data[row][locus] = lastAlleleCode[locus];
			}
			else {
				data->data[row][locus] = allele2code[locus][allele];
			}
		}
	}
	fin.close();
	return data;
}

/*
structure_data *readData_stru(string infile, int sort, int ndcols, int ndrows, int nrows, int ncols, string STRU_MISSING) {
	igzstream fin;

	fin.open(infile.c_str());

	if (fin.fail())
	{
		cerr << "Could not open " << infile << " for reading.\n";
		throw - 1;
	}

	structure_data *data = new structure_data;
	string line;
	int nind = nrows / 2;
	data->nind = nind;
	data->data = new short*[nrows];
	for (int i = 0; i < nrows; i++) {
		data->data[i] = new short[ncols];
	}
	data->ind_names = new string[nind];

	//toss out non-data rows
	for (int i = 0; i < ndrows; i++) getline(fin, line);

	data->nloci = ncols;
	string field, key, allele;
	//map<string, int> ind2index;
	map<string, short> *allele2code = new map<string, short>[ncols];
	short *lastAlleleCode = new short[ncols];
	for (int i = 0; i < ncols; i++) {
		allele2code[i][STRU_MISSING] = -9;
		lastAlleleCode[i] = -1;
	}
	int index;
	int indCount = -1;
	for (int row = 0; row < nrows; row++) {
		for (int i = 1; i <= ndcols; i++) {
			fin >> field;
			if (i == sort) key = field;
			if (row % 2 == 0) {
				data->ind_names[row / 2] = key;
			}
		}

		for (int locus = 0; locus < ncols; locus++) {
			fin >> allele;
			if (allele2code[locus].count(allele) == 0) {
				lastAlleleCode[locus]++;
				allele2code[locus][allele] = lastAlleleCode[locus];
				data->data[row][locus] = lastAlleleCode[locus];
			}
			else {
				data->data[row][locus] = allele2code[locus][allele];
			}
		}
	}
	fin.close();
	return data;
}
*/
/*
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
*/

structure_data *readData_tped_tfam2(string tped_filename, string tfam_filename, int &nrow, int &nloci, string TPED_MISSING) {
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;

	pedin.open(tped_filename.c_str());
	if (pedin.fail())
	{
		LOG.err("ERROR: Coult not open", tped_filename);
		throw - 1;
	}

	famin.open(tfam_filename.c_str());
	if (famin.fail())
	{
		LOG.err("ERROR: Coult not open", tfam_filename);
		throw - 1;
	}

	while (getline(famin, junk)) nrow += 2;
	famin.close();
	famin.clear();

	while (getline(pedin, junk)) nloci++;
	pedin.close();
	pedin.clear();

	pedin.open(tped_filename.c_str());
	famin.open(tfam_filename.c_str());

	structure_data *data = new structure_data;
	int nind = nrow / 2;
	data->nind = nind;
	data->data = new short*[2 * nind];
	for (int i = 0; i < nrow; i++) {
		data->data[i] = new short[nloci];
	}
	data->ind_names = new string[nind];

	for (int i = 0; i < nind; i++)
	{
		famin >> junk;
		famin >> data->ind_names[i];
		getline(famin, junk);
	}
	famin.close();

	data->nloci = nloci;

	map<string, short> allele2code;
	short lastAlleleCode = -1;

	string allele;
	for (int locus = 0; locus < nloci; locus++)
	{
		allele2code.clear();
		allele2code[TPED_MISSING] = -9;
		lastAlleleCode = -1;

		pedin >> junk;
		pedin >> junk;
		pedin >> junk;
		pedin >> junk;

		for (int i = 0; i < 2 * nind; i++)
		{
			pedin >> allele;
			if (allele2code.count(allele) == 0) {
				lastAlleleCode++;
				allele2code[allele] = lastAlleleCode;
				data->data[i][locus] = lastAlleleCode;
			}
			else {
				data->data[i][locus] = allele2code[allele];
			}
		}
	}

	pedin.close();
	return data;
}


structure_data *readData_tped_tfam(string tped_filename, string tfam_filename, int &nrow, int &nloci, string TPED_MISSING)
{
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;

	pedin.open(tped_filename.c_str());
	if (pedin.fail())
	{
		LOG.err("ERROR: Coult not open", tped_filename);
		throw - 1;
	}

	famin.open(tfam_filename.c_str());
	if (famin.fail())
	{
		LOG.err("ERROR: Coult not open", tfam_filename);
		throw - 1;
	}

	while (getline(famin, junk)) nrow += 2;
	famin.close();
	famin.clear();

	while (getline(pedin, junk)) nloci++;
	pedin.close();
	pedin.clear();

	pedin.open(tped_filename.c_str());
	famin.open(tfam_filename.c_str());

	structure_data *data = new structure_data;

	int nind = nrow / 2;

	//cerr << "Reading " << nind << " diploid individuals at " << nloci << " loci.\n";

	data->nind = nind;
	data->data = new short*[nind];
	for (int i = 0; i < nind; i++) data->data[i] = new short[nloci];
	data->ind_names = new string[nind];

	for (int i = 0; i < nind; i++)
	{
		famin >> junk;
		famin >> data->ind_names[i];
		getline(famin, junk);
	}

	data->nloci = nloci;

	map<string, short> allele2code;
	short lastAlleleCode = -1;

	string allele1, allele2;
	short genotypeCode;
	for (int locus = 0; locus < nloci; locus++)
	{
		allele2code.clear();
		allele2code[TPED_MISSING] = -9;
		lastAlleleCode = -1;

		pedin >> junk;
		pedin >> junk;
		pedin >> junk;
		pedin >> junk;

		for (int ind = 0; ind < nind; ind++)
		{
			//alleleCount = 0;
			pedin >> allele1;
			pedin >> allele2;

			if (allele2code.count(allele1) == 0) {
				lastAlleleCode++;
				if (lastAlleleCode > 1) {
					LOG.err("ERORR: --biallelic flag set, but found more than 2 alleles at locus", locus + 1);
					throw - 1;
				}
				allele2code[allele1] = lastAlleleCode;
			}

			if (allele2code.count(allele2) == 0) {
				lastAlleleCode++;
				if (lastAlleleCode > 1) {
					LOG.err("ERORR: --biallelic flag set, but found more than 2 alleles at locus", locus + 1);
					throw - 1;
				}
				allele2code[allele2] = lastAlleleCode;
			}

			genotypeCode = allele2code[allele1] + allele2code[allele2];
			if (genotypeCode < 0) {
				data->data[ind][locus] = allele2code[TPED_MISSING];
			}
			else {
				data->data[ind][locus] = genotypeCode;
			}
		}
	}

	return data;
}


void readData_pop_freq(igzstream & fin, structure_data & data,
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

int countFields(const string & str)
{
	string::const_iterator it;
	int result;
	int numFields = 0;
	int seenChar = 0;
	for (it = str.begin() ; it < str.end(); it++)
	{
		result = isspace(*it);
		if (result == 0 && seenChar == 0)
		{
			numFields++;
			seenChar = 1;
		}
		else if (result != 0)
		{
			seenChar = 0;
		}
	}
	return numFields;
}

int search(string * s, int size, string key)
{
	for (int i = 0; i < size; i++)
	{
		if (s[i].compare(key) == 0) return i;
	}
	return -9;
}

int put(string * s, int size, string key)
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
short int *split_int(igzstream & fin, int fields)
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

