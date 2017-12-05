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

void combine_partial_files(param_t *params) {

	string outfile = params->getStringFlag(ARG_OUTFILE);
	bool PRINT_LOG = params->getBoolFlag(ARG_LOG);
	bool PRINT_LONG = params->getBoolFlag(ARG_LONG_FORMAT);
	bool PRINT_LONG_IBS = params->getBoolFlag(ARG_IBS_LONG);
	vector<string> comboFiles = params->getStringListFlag(ARG_COMBINE);

	int nind;
	string type;
	string *ind_names;

	LOG.log("Combining", int(comboFiles.size()), false);
	LOG.log(" partial files.");
	for (int i = 0; i < comboFiles.size(); i++) {
		igzstream fin;
		istringstream ssin;

		int fileNind;
		string fileType;
		string junk, line;
		double dnum;
		double inum;

		fin.open(comboFiles[i].c_str());
		if (fin.fail()) {
			LOG.err("ERROR: Can not open", comboFiles[i]);
			throw 0;
		}

		LOG.log("Reading", comboFiles[i], false);

		ssin.clear();
		getline(fin, line);
		ssin.str(line);
		ssin >> fileType;
		ssin >> fileNind;
		if (fileNind <= 0) {
			LOG.err("ERROR: Bad number of individuals:", fileNind);
			LOG.err("ERROR: Check", comboFiles[i]);
			throw 0;
		}

		LOG.log(" of type", fileType, false);
		LOG.log(" with", fileNind, false);
		LOG.log(" individuals.");

		if (i == 0) {
			type = fileType;
			nind = fileNind;
			init_storage(nind, false);
			ind_names = new string[nind];
		}

		if (nind != fileNind) {
			LOG.err("ERROR: Number of individuals did not match across files. Check", comboFiles[i]);
			throw 0;
		}

		if (fileType.compare(type) != 0) {
			LOG.err("ERROR: Can not combine files of different types:", type, false);
			LOG.err(" ", fileType);
			throw 0;
		}

		for (int row = 0; row < nind; row++) {
			ssin.clear();
			getline(fin, line);
			ssin.str(line);
			for (int col = 0; col < nind; col++) {
				ssin >> inum;
				NUM_LOCI[row][col] += inum;
			}
		}

		ssin.clear();
		getline(fin, line);
		ssin.str(line);
		if ( i == 0 ) {
			for (int ind = 0; ind < nind; ind++) {
				ssin >> ind_names[ind];
			}
		}

		for (int row = 0; row < nind; row++) {
			ssin.clear();
			getline(fin, line);
			ssin.str(line);
			for (int col = 0; col < nind + 1; col++) {
				if (col == 0) ssin >> junk;
				else {
					ssin >> dnum;
					DIST_MAT[row][col - 1] += dnum;
				}
			}
		}

		fin.close();
		fin.clear();
	}

	ofstream out;
	if (type.compare("dist") == 0) outfile += ".asd.dist";
	else outfile += "." + type;

	out.open(outfile.c_str());
	if (out.fail()) {
		LOG.err("ERROR: Could not open", outfile);
		throw 0;
	}

	if (!PRINT_LONG && !PRINT_LONG_IBS) {
		for (int i = 0; i < nind; i++) out << ind_names[i] << " ";
		out << endl;

		for (int i = 0; i < nind; i++) {
			out << ind_names[i] << " ";
			for (int j = 0; j < nind; j++) {
				if (type.compare("dist") == 0) {
					if (PRINT_LOG) out << 0 - log(double(DIST_MAT[i][j]) / double(NUM_LOCI[i][j])) << " ";
					else out << 1 - (double(DIST_MAT[i][j]) / double(NUM_LOCI[i][j])) << " ";
				}
				else {
					out << double(DIST_MAT[i][j]) / double(NUM_LOCI[i][j]) << " ";
				}
			}
			out << endl;
		}
	}
	else {
		for (int i = 0; i < nind; i++) {
			for (int j = i; j < nind; j++) {
				out << ind_names[i] << " " << ind_names[j] << " ";
				if (type.compare("dist") == 0) {
					if (PRINT_LOG) out << 0 - log(double(DIST_MAT[i][j]) / double(NUM_LOCI[i][j])) << endl;
					else out << 1 - (double(DIST_MAT[i][j]) / double(NUM_LOCI[i][j])) << endl;
				}
				else {
					out << double(DIST_MAT[i][j]) / double(NUM_LOCI[i][j]) << endl;
				}
			}
		}
	}
	out.close();
	LOG.log("Wrote results to", outfile);
	return;
}

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

bool finalize_calculations(int nind, int ncols, bool CALC_ALL_IBS, bool GRM) {
	for (int i = 0; i < nind; i++)
	{
		for (int j = i; j < nind; j++)
		{
			if (i == j && !GRM)
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

void write_ibs_matrices(string outfile, int nind, int ncols, string *ind_names, bool PRINT_PARTIAL, bool PRINT_LONG_IBS) {
	string ibs_fname[3];
	ibs_fname[0] = outfile + ".ibs0";
	ibs_fname[1] = outfile + ".ibs1";
	ibs_fname[2] = outfile + ".ibs2";

	if (PRINT_PARTIAL) {
		for (int i = 0; i < 3; i++) ibs_fname[i] += ".partial";
	}

	ofstream out;
	string type;
	for (int ibs = 0; ibs <= 2; ibs++) {

		int **mat;
		if (ibs == 0) {
			mat = IBS_0_MAT;
			type = "ibs0";
		}
		else if (ibs == 1) {
			mat = IBS_1_MAT;
			type = "ibs1";
		}
		else if (ibs == 2) {
			mat = IBS_2_MAT;
			type = "ibs2";
		}

		out.open(ibs_fname[ibs].c_str());
		if (out.fail()) {
			LOG.err("ERROR: Could not open", ibs_fname[ibs]);
			throw 0;
		}
		out << std::fixed;
		if (PRINT_PARTIAL) {
			out << std::setprecision(0);
			out << type << " " << nind << endl;
			for (int i = 0; i < nind; i++) {
				for (int j = 0; j < nind; j++) {
					out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
				}
				out << endl;
			}
			//out << endl;
		}

		if (!PRINT_LONG_IBS) {
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
		}
		else {
			for (int i = 0; i < nind; i++) {
				for (int j = i; j < nind ; j++) {
					out << ind_names[i] << " " << ind_names[j] << " ";
					out << double(mat[i][j]) / (double(ncols) + double(NUM_LOCI[i][j])) << endl;
				}
			}
		}
		out.close();
		out.clear();
		LOG.log("Wrote results to", ibs_fname[ibs]);
	}

	return;
}

void write_dist_matrix(string outfile, int nind, int ncols, string *ind_names, bool PRINT_PARTIAL, bool PRINT_LOG, bool PRINT_LONG, bool GRM) {
	ofstream out;
	string type = (GRM ? "grm" : "dist");
	if (!PRINT_PARTIAL) {
		outfile += (GRM ? ".grm" : ".asd.dist");
	}
	else {
		outfile += (GRM ? ".grm.partial" : ".asd.partial");
	}

	out.open(outfile.c_str());
	if (out.fail()) {
		LOG.err("ERROR: Could not open", outfile);
		throw 0;
	}
	out << std::fixed;
	if (PRINT_PARTIAL)
	{
		out << std::setprecision(0);
		out << type << " " << nind << endl;
		for (int i = 0; i < nind; i++)
		{
			for (int j = 0; j < nind; j++)
			{
				out << double(ncols) + double(NUM_LOCI[i][j]) << " ";
			}
			out << endl;
		}
		out << std::setprecision(6);
		//out << endl;
	}

	if (!PRINT_LONG) {
		for (int i = 0; i < nind; i++) out << ind_names[i] << " ";
		out << endl;

		for (int i = 0; i < nind; i++){
			out << ind_names[i] << " ";
			for (int j = 0; j < nind ; j++){
				if (PRINT_LOG){
					out << 0 - log(double(DIST_MAT[i][j]) /
					               (double(ncols) + double(NUM_LOCI[i][j])))
					    << " ";
				}
				else if (!PRINT_PARTIAL && !GRM){
					out << 1 - (double(DIST_MAT[i][j]) /
					            (double(ncols) + double(NUM_LOCI[i][j])))
					    << " ";
				}
				else if (!PRINT_PARTIAL && GRM){
					out << (double(DIST_MAT[i][j]) /
					            (double(ncols) + double(NUM_LOCI[i][j])))
					    << " ";
				}
				else out << DIST_MAT[i][j] << " ";
			}
			out << endl;
		}
	}
	else {
		for (int i = 0; i < nind; i++) {
			for (int j = i; j < nind ; j++) {
				out << ind_names[i] << " " << ind_names[j] << " ";
				if (PRINT_LOG) {
					out << 0 - log(double(DIST_MAT[i][j]) /
					               (double(ncols) + double(NUM_LOCI[i][j])))
					    << endl;
				}
				else if(GRM){
					out << (double(DIST_MAT[i][j]) /
					            (double(ncols) + double(NUM_LOCI[i][j])))
					    << endl;
				}
				else {
					out << 1 - (double(DIST_MAT[i][j]) /
					            (double(ncols) + double(NUM_LOCI[i][j])))
					    << endl;
				}
			}
		}
	}
	out.close();
	LOG.log("Wrote results to", outfile);
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

vector<double>* biallelic_maf_stru(string infile, int sort, int &nrows, int &ncols, int &ndrows, int &ndcols, int &nkeep, double MAF, string STRU_MISSING) {
	string line;
	igzstream fin;
	ndrows = -1;
	ndcols = -1;
	nrows = 0;
	ncols = 0;
	fin.open(infile.c_str());

	if (fin.fail()) {
		LOG.err("ERROR: Could not open", infile);
		throw - 1;
	}

	int totalRows = 0;
	int currentCols = 0;
	stringstream ss;
	map<string, int> *allele2count;
	string allele;
	while (getline(fin, line)) {
		currentCols = countFields(line);
		if (totalRows == 0) {
			ncols = currentCols;

			if (!check_int_gt_0(ncols)) {
				LOG.err("ERROR: Number of loci must be > 0. Found", ncols);
				throw 0;
			}
			LOG.log("Number of loci:", ncols);

			allele2count = new map<string, int>[ncols];
			for (int i = 0; i < ncols; i++) allele2count[i][STRU_MISSING] = 0;
		}
		if (ndrows < 0 && currentCols != ncols) {
			ndcols = currentCols - ncols;
			ndrows = totalRows;

			if (!check_int_gt_0(ndrows)) {
				LOG.err("ERROR: Non-data rows must be > 0. Found", ndrows);
				throw 0;
			}
			LOG.log("Non-data header rows:", ndrows);

			if (!check_int_gt_0(ndcols)) {
				LOG.err("ERROR: Non-data columns must be > 0. Found", ndcols);
				throw 0;
			}
			LOG.log("Non-data header columns:", ndcols);

			if (!check_sort(sort, ndcols)) {
				LOG.err("ERROR: Individual ID column must be in [ 1, ", ndcols, false);
				LOG.err("].");
				throw 0;
			}
		}
		if (ndrows > 0 || (ndrows < 0 && currentCols != ncols) ) {
			ss.str(line);
			for (int i = 0; i < ndcols; i++) ss >> line;
			for (int i = 0; i < ncols; i++) {
				ss >> allele;
				if (allele2count[i].count(allele) == 0) {
					allele2count[i][allele] = 1;
				}
				else {
					allele2count[i][allele]++;
				}
				if (allele2count[i].size() > 3) {
					LOG.err("ERROR: --biallelic flag set but found more than two alleles at locus ", i);
					throw 0;
				}
			}
			ss.clear();
		}
		totalRows++;
	}
	nrows = totalRows - ndrows;

	if (!check_int_gt_0(nrows)) {
		LOG.err("ERROR: Number of samples must be > 0. Found", nrows);
		throw 0;
	}

	vector<double> *AF0 = new vector<double>;
	//double MAF_HOM = biallelic_ehom(MAF);
	for (int i = 0; i < ncols; i++) {
		int n = 0;
		double maf = 0;
		for (std::map<string, int>::iterator it = allele2count[i].begin(); it != allele2count[i].end(); ++it) {
			if ((it->first).compare(STRU_MISSING) != 0) {
				n += it->second;
				maf = it->second;
			}
		}
		maf /= n;
		if (maf <= 1 - MAF && maf >= MAF) nkeep++;
		AF0->push_back(maf);
	}

	fin.close();

	LOG.log("Loci filtered:", ncols - nkeep);
	delete [] allele2count;
	return AF0;
}

structure_data *readData_stru(string infile, int sort, int &nrows, int &ncols, double MAF, string STRU_MISSING) {
	int ndcols = -1;
	int ndrows = -1;
	nrows = -1;
	ncols = -1;
	string line;
	igzstream fin;
	int nkeep = 0;

	vector<double> *AF0 = biallelic_maf_stru(infile, sort, nrows, ncols, ndrows, ndcols, nkeep, MAF, STRU_MISSING);

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

	data->nloci = nkeep;

	map<string, short> *allele2code = new map<string, short>[nkeep];
	short *lastAlleleCode = new short[nkeep];
	for (int i = 0; i < nkeep; i++) {
		allele2code[i][STRU_MISSING] = -9;
		lastAlleleCode[i] = -1;
	}

	int index;
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

		int i = 0;
		for (int locus = 0; locus < ncols; locus++) {
			sin1 >> allele1;
			sin2 >> allele2;
			if (AF0->at(locus) >= MAF && AF0->at(locus) <= 1 - MAF) {
				if (allele2code[i].count(allele1) == 0) {
					lastAlleleCode[i]++;
					allele2code[i][allele1] = lastAlleleCode[i];
				}

				if (allele2code[i].count(allele2) == 0) {
					lastAlleleCode[i]++;
					allele2code[i][allele2] = lastAlleleCode[i];
				}

				genotypeCode = allele2code[i][allele1] + allele2code[i][allele2];
				if (genotypeCode < 0) {
					data->data[ind][i] = allele2code[i][STRU_MISSING];
				}
				else {
					data->data[ind][i] = genotypeCode;
				}
				i++;
			}
		}
	}
	fin.close();
	ncols = nkeep;
	delete AF0;
	return data;
}

structure_data *readData_stru2(string infile, int sort, int &nrows, int &ncols, double MAF, string STRU_MISSING) {
	string line;
	int ndcols = -1;
	int ndrows = -1;
	nrows = -1;
	ncols = -1;
	int nkeep = 0;
	double MAF_HOM = biallelic_ehom(MAF);
	vector<double>* EHOM = multiallelic_maf_stru(infile, sort, nrows, ncols, ndrows, ndcols, nkeep, MAF, STRU_MISSING);

	igzstream fin;
	fin.open(infile.c_str());

	structure_data *data = new structure_data;
	int nind = nrows / 2;
	data->nind = nind;
	data->data = new short*[nrows];
	for (int i = 0; i < nrows; i++) {
		data->data[i] = new short[nkeep];
	}
	data->ind_names = new string[nind];

	//toss out non-data rows
	for (int i = 0; i < ndrows; i++) getline(fin, line);

	data->nloci = nkeep;
	string field, key, allele;
	map<string, short> *allele2code = new map<string, short>[nkeep];
	short *lastAlleleCode = new short[nkeep];
	for (int i = 0; i < nkeep; i++) {
		allele2code[i][STRU_MISSING] = -9;
		lastAlleleCode[i] = -1;
	}
	int index;
	for (int row = 0; row < nrows; row++) {

		for (int i = 1; i <= ndcols; i++) {
			fin >> field;
			if (i == sort) key = field;
			if (row % 2 == 0) {
				data->ind_names[row / 2] = key;
			}
		}

		index = 0;
		for (int locus = 0; locus < ncols; locus++) {
			fin >> allele;
			if (EHOM->at(locus) <= MAF_HOM) {
				if (allele2code[index].count(allele) == 0) {
					lastAlleleCode[index]++;
					allele2code[index][allele] = lastAlleleCode[index];
					data->data[row][index] = lastAlleleCode[index];
				}
				else {
					data->data[row][index] = allele2code[index][allele];
				}
				index++;
			}
		}
	}
	fin.close();
	ncols = nkeep;
	delete EHOM;
	return data;
}

vector<double>* multiallelic_maf_stru(string infile, int sort, int &nrows, int &ncols, int &ndrows, int &ndcols, int &nkeep, double MAF, string STRU_MISSING) {
	string line;
	igzstream fin;
	ndrows = -1;
	ndcols = -1;
	nrows = 0;
	ncols = 0;
	fin.open(infile.c_str());

	if (fin.fail()) {
		LOG.err("ERROR: Could not open", infile);
		throw - 1;
	}

	int totalRows = 0;
	int currentCols = 0;
	stringstream ss;
	map<string, int> *allele2count;
	string allele;
	while (getline(fin, line)) {
		currentCols = countFields(line);
		if (totalRows == 0) {
			ncols = currentCols;

			if (!check_int_gt_0(ncols)) {
				LOG.err("ERROR: Number of loci must be > 0. Found", ncols);
				throw 0;
			}
			LOG.log("Number of loci:", ncols);

			allele2count = new map<string, int>[ncols];
		}
		if (ndrows < 0 && currentCols != ncols) {
			ndcols = currentCols - ncols;
			ndrows = totalRows;

			if (!check_int_gt_0(ndrows)) {
				LOG.err("ERROR: Non-data rows must be > 0. Found", ndrows);
				throw 0;
			}
			LOG.log("Non-data header rows:", ndrows);

			if (!check_int_gt_0(ndcols)) {
				LOG.err("ERROR: Non-data columns must be > 0. Found", ndcols);
				throw 0;
			}
			LOG.log("Non-data header columns:", ndcols);

			if (!check_sort(sort, ndcols)) {
				LOG.err("ERROR: Individual ID column must be in [ 1, ", ndcols, false);
				LOG.err("].");
				throw 0;
			}
		}
		if (ndrows > 0 || (ndrows < 0 && currentCols != ncols) ) {
			ss.str(line);
			for (int i = 0; i < ndcols; i++) ss >> line;
			for (int i = 0; i < ncols; i++) {
				ss >> allele;
				if (allele2count[i].count(allele) == 0) {
					allele2count[i][allele] = 1;
				}
				else {
					allele2count[i][allele]++;
				}
			}
			ss.clear();
		}
		totalRows++;
	}
	nrows = totalRows - ndrows;

	if (!check_int_gt_0(nrows)) {
		LOG.err("ERROR: Number of samples must be > 0. Found", nrows);
		throw 0;
	}

	vector<double> *EHOM = new vector<double>;
	double MAF_HOM = biallelic_ehom(MAF);
	for (int i = 0; i < ncols; i++) {
		int n = 0;
		double ehom = 0;
		for (std::map<string, int>::iterator it = allele2count[i].begin(); it != allele2count[i].end(); ++it) {
			if ((it->first).compare(STRU_MISSING) != 0) {
				n += it->second;
			}
		}
		for (std::map<string, int>::iterator it = allele2count[i].begin(); it != allele2count[i].end(); ++it) {
			if ((it->first).compare(STRU_MISSING) != 0) {
				ehom += (double(it->second) / double(n)) * (double(it->second) / double(n));
			}
		}
		if (ehom <= MAF_HOM) nkeep++;
		EHOM->push_back(ehom);
	}

	fin.close();

	LOG.log("Loci filtered:", ncols - nkeep);

	return EHOM;
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

vector<double>* multiallelic_maf_tped(string tped_filename, string tfam_filename, int &nrow, int &nloci, int &nkeep, double MAF, string TPED_MISSING) {
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;
	nkeep = 0;
	double MAF_HOM = biallelic_ehom(MAF);

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

	vector<double> *EHOM = new vector<double>;
	map<string, int> allele2count;
	stringstream ss;
	string allele;
	int n = 0;
	double ehom = 0;

	while (getline(pedin, junk)) {
		ss.str(junk);
		ss >> junk;
		ss >> junk;
		ss >> junk;
		ss >> junk;
		n = 0;
		allele2count.clear();
		for (int i = 0; i < nrow; i++) {
			ss >> allele;
			if (allele.compare(TPED_MISSING) == 0) {
				continue;
			}

			if (allele2count.count(allele) == 0) {
				allele2count[allele] = 1;
			}
			else {
				allele2count[allele]++;
			}

			n++;
		}

		ss.clear();
		ehom = 0;
		for (std::map<string, int>::iterator it = allele2count.begin(); it != allele2count.end(); ++it) {
			ehom += (double(it->second) / double(n)) * (double(it->second) / double(n));
		}

		EHOM->push_back(ehom);
		if (EHOM->at(nloci) <= MAF_HOM) nkeep++;

		nloci++;
	}
	pedin.close();

	LOG.log("Loci filtered:", nloci - nkeep);

	return EHOM;
}

structure_data *readData_tped_tfam2(string tped_filename, string tfam_filename, int &nrow, int &nloci, double MAF, string TPED_MISSING) {
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;
	int nkeep = 0;
	double MAF_HOM = biallelic_ehom(MAF);

	vector<double>* EHOM = multiallelic_maf_tped(tped_filename, tfam_filename, nrow, nloci, nkeep, MAF, TPED_MISSING);

	pedin.open(tped_filename.c_str());
	famin.open(tfam_filename.c_str());

	structure_data *data = new structure_data;
	int nind = nrow / 2;
	data->nind = nind;
	data->data = new short*[2 * nind];
	for (int i = 0; i < nrow; i++) {
		data->data[i] = new short[nkeep];
	}
	data->ind_names = new string[nind];

	for (int i = 0; i < nind; i++) {
		famin >> junk;
		famin >> data->ind_names[i];
		getline(famin, junk);
	}
	famin.close();

	data->nloci = nkeep;

	map<string, short> allele2code;
	short lastAlleleCode = -1;

	string allele;
	int index = 0;
	for (int locus = 0; locus < nloci; locus++) {
		getline(pedin, junk);
		if (EHOM->at(locus) <= MAF_HOM) {
			stringstream ss;
			ss.str(junk);
			allele2code.clear();
			allele2code[TPED_MISSING] = -9;
			lastAlleleCode = -1;

			ss >> allele;
			ss >> allele;
			ss >> allele;
			ss >> allele;

			for (int i = 0; i < 2 * nind; i++) {
				ss >> allele;
				if (allele2code.count(allele) == 0) {
					lastAlleleCode++;
					allele2code[allele] = lastAlleleCode;
					data->data[i][index] = lastAlleleCode;
				}
				else {
					data->data[i][index] = allele2code[allele];
				}
			}
			index++;
		}
	}

	pedin.close();
	nloci = nkeep;
	delete EHOM;
	return data;
}

vector<double>* biallelic_maf_tped(string tped_filename, string tfam_filename, int &nrow, int &nloci, int &nkeep, double MAF, string TPED_MISSING) {
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;
	nkeep = 0;

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

	vector<double> *AF0 = new vector<double>;
	map<string, int> allele2count;
	stringstream ss;
	string allele;
	int n = 0;

	while (getline(pedin, junk)) {
		ss.str(junk);
		ss >> junk;
		ss >> junk;
		ss >> junk;
		ss >> junk;
		n = 0;
		allele2count.clear();
		for (int i = 0; i < nrow; i++) {
			ss >> allele;
			if (allele.compare(TPED_MISSING) == 0) {
				continue;
			}

			if (allele2count.count(allele) == 0) {
				allele2count[allele] = 1;
			}
			else {
				allele2count[allele]++;
			}

			n++;
		}
		ss.clear();
		map<string, int>::iterator it = allele2count.begin();

		AF0->push_back(double(it->second) / double(n));

		if (AF0->at(nloci) >= MAF && AF0->at(nloci) <= 1 - MAF) nkeep++;

		nloci++;
	}
	pedin.close();

	LOG.log("Loci filtered:", nloci - nkeep);

	return AF0;
}

structure_data *readData_tped_tfam(string tped_filename, string tfam_filename, int &nrow, int &nloci, double MAF, string TPED_MISSING) {
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;
	int nkeep = 0;
	vector<double> *AF0 = biallelic_maf_tped(tped_filename, tfam_filename, nrow, nloci, nkeep, MAF, TPED_MISSING);

	pedin.open(tped_filename.c_str());
	famin.open(tfam_filename.c_str());

	structure_data *data = new structure_data;

	int nind = nrow / 2;

	data->nind = nind;
	data->data = new short*[nind];
	for (int i = 0; i < nind; i++) data->data[i] = new short[nkeep];
	data->ind_names = new string[nind];

	for (int i = 0; i < nind; i++)
	{
		famin >> junk;
		famin >> data->ind_names[i];
		getline(famin, junk);
	}

	data->nloci = nkeep;

	map<string, short> allele2code;
	short lastAlleleCode = -1;

	string allele1, allele2;
	short genotypeCode;
	int index = 0;
	for (int locus = 0; locus < nloci; locus++)
	{
		getline(pedin, junk);
		if (AF0->at(locus) >= MAF && AF0->at(locus) <= 1 - MAF) {
			stringstream ss;
			ss.str(junk);
			allele2code.clear();
			allele2code[TPED_MISSING] = -9;
			lastAlleleCode = -1;

			ss >> allele1;
			ss >> allele1;
			ss >> allele1;
			ss >> allele1;

			for (int ind = 0; ind < nind; ind++)
			{
				//alleleCount = 0;
				ss >> allele1;
				ss >> allele2;

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
					data->data[ind][index] = allele2code[TPED_MISSING];
				}
				else {
					data->data[ind][index] = genotypeCode;
				}
			}
			index++;
		}
	}

	nloci = nkeep;
	pedin.close();
	delete AF0;
	return data;
}

vector<double>* biallelic_maf_vcf(string vcf_filename, int &nrow, int &nloci, int &nkeep, int &numComments, double MAF) {
	string line;
	nrow = 0;
	nloci = 0;
	numComments = 0;
	igzstream fin;

	nkeep = 0;

	fin.open(vcf_filename.c_str());
	if (fin.fail())
	{
		LOG.err("ERROR: Coult not open", vcf_filename);
		throw - 1;
	}

	vector<double> *AF0 = new vector<double>;

	int previous_nind = -1;
	int current_nind = 0;
	stringstream ss;
	string gt;
	int n = 0;
	int n0 = 0;
	string a1, a2;
	while (getline(fin, line)) {
		if (line[0] == '#') {
			numComments++;
			continue;
		}
		current_nind = countFields(line);
		if (previous_nind >= 0 && previous_nind != current_nind)
		{
			cerr << "ERROR: line " << nloci + 1 << " of " << vcf_filename << " has " << current_nind
			     << " fields, but the previous line has " << previous_nind << " fields.\n";
			throw 0;
		}
		previous_nind = current_nind;
		nrow = 2 * (current_nind - 9);
		ss.str(line);
		n = 0;
		n0 = 0;
		for (int i = 0; i < current_nind; i++) {
			ss >> gt;
			if (i < 9) continue;
			if (gt.compare("./.") == 0 || gt.compare(".|.") == 0 || gt.compare(".") == 0) {
				continue;
			}

			a1 = gt[0];
			a2 = gt[2];

			if (atoi(a1.c_str()) != 0 && atoi(a1.c_str()) != 1 && atoi(a2.c_str()) != 1 && atoi(a2.c_str()) != 0) {
				LOG.err("ERORR: --biallelic flag set, but found more than 2 alleles at locus", nloci + 1);
				throw - 1;
			}

			if (atoi(a1.c_str()) == 0) n0++;
			if (atoi(a2.c_str()) == 0) n0++;
			n += 2;
		}
		ss.clear();
		AF0->push_back(double(n0) / double(n));
		if (AF0->at(nloci) >= MAF && AF0->at(nloci) <= 1 - MAF) nkeep++;
		nloci++;
	}

	fin.close();
	fin.clear();

	LOG.log("Loci filtered:", nloci - nkeep);

	return AF0;
}

vector<double> *multiallelic_maf_vcf(string vcf_filename, int &nrow, int &nloci, int &nkeep, int &numComments, double MAF) {
	string line;
	nrow = 0;
	nloci = 0;
	numComments = 0;
	igzstream fin;

	nkeep = 0;
	double MAF_HOM = biallelic_ehom(MAF);

	fin.open(vcf_filename.c_str());
	if (fin.fail())
	{
		LOG.err("ERROR: Coult not open", vcf_filename);
		throw - 1;
	}

	vector<double> *EHOM = new vector<double>;
	map<string, int> allele2count;
	int previous_nind = -1;
	int current_nind = 0;
	stringstream ss;
	string gt;
	int n = 0;
	double ehom = 0;
	string a1, a2;
	while (getline(fin, line)) {
		if (line[0] == '#') {
			numComments++;
			continue;
		}
		current_nind = countFields(line);
		if (previous_nind >= 0 && previous_nind != current_nind)
		{
			cerr << "ERROR: line " << nloci + 1 << " of " << vcf_filename << " has " << current_nind
			     << " fields, but the previous line has " << previous_nind << " fields.\n";
			throw 0;
		}
		previous_nind = current_nind;
		nrow = 2 * (current_nind - 9);
		ss.str(line);
		n = 0;
		allele2count.clear();
		for (int i = 0; i < current_nind; i++) {
			ss >> gt;
			if (i < 9) continue;
			if (gt.compare("./.") == 0 || gt.compare(".|.") == 0 || gt.compare(".") == 0) {
				continue;
			}

			a1 = gt[0];
			a2 = gt[2];

			if (allele2count.count(a1) == 0) {
				allele2count[a1] = 1;
			}
			else {
				allele2count[a1]++;
			}

			if (allele2count.count(a2) == 0) {
				allele2count[a2] = 1;
			}
			else {
				allele2count[a2]++;
			}

			n += 2;
		}
		ss.clear();
		ehom = 0;
		for (std::map<string, int>::iterator it = allele2count.begin(); it != allele2count.end(); ++it) {
			ehom += (double(it->second) / double(n)) * (double(it->second) / double(n));
		}

		EHOM->push_back(ehom);

		if (EHOM->at(nloci) <= MAF_HOM) nkeep++;

		nloci++;
	}

	fin.close();
	fin.clear();

	LOG.log("Loci filtered:", nloci - nkeep);

	return EHOM;
}

double biallelic_ehom(double MAF) {
	return MAF * MAF + (1 - MAF) * (1 - MAF);
}

structure_data *readData_vcf(string vcf_filename, int &nrow, int &nloci, double MAF) {
	int numComments = 0;
	string line;
	int nkeep = 0;
	vector<double> *AF0 = biallelic_maf_vcf(vcf_filename, nrow, nloci, nkeep, numComments, MAF);
	int nind = nrow / 2;

	igzstream fin;
	fin.open(vcf_filename.c_str());

	structure_data *data = new structure_data;
	data->nloci = nkeep;
	data->nind = nind;
	data->data = new short*[nind];
	for (int i = 0; i < nind; i++) data->data[i] = new short[nkeep];
	data->ind_names = new string[nind];

	for (int i = 0; i < numComments; i++) getline(fin, line);
	stringstream ss;
	ss.str(line);
	for (int i = 0; i < 9; i++) ss >> line;
	for (int i = 0; i < nind; i++) ss >> data->ind_names[i];

	string gt, a1, a2;
	int i = 0;
	for (int locus = 0; locus < nloci; locus++)
	{
		getline(fin, line);
		if (AF0->at(locus) >= MAF && AF0->at(locus) <= 1 - MAF) {
			ss.clear();
			ss.str(line);
			for (int field = 0; field < 9; field++) ss >> gt;
			for (int ind = 0; ind < nind; ind++) {
				ss >> gt;
				if (gt.compare("./.") == 0 || gt.compare(".|.") == 0 || gt.compare(".") == 0) {
					data->data[ind][i] = -9;
					i++;
					continue;
				}
				else {
					a1 = gt[0];
					a2 = gt[2];
				}
				data->data[ind][i] = atoi(a1.c_str()) + atoi(a2.c_str());
			}
			i++;
		}
	}

	nloci = nkeep;
	fin.close();
	delete AF0;
	return data;
}

structure_data *readData_vcf2(string vcf_filename, int &nrow, int &nloci, double MAF) {
	int numComments = 0;
	string line;
	double MAF_HOM = biallelic_ehom(MAF);
	int nkeep = 0;
	vector<double> *EHOM = multiallelic_maf_vcf(vcf_filename, nrow, nloci, nkeep, numComments, MAF);
	int nind = nrow / 2;

	igzstream fin;
	fin.open(vcf_filename.c_str());

	structure_data *data = new structure_data;
	data->nloci = nkeep;
	data->nind = nind;
	data->data = new short*[2 * nind];
	for (int i = 0; i < 2 * nind; i++) data->data[i] = new short[nkeep];
	data->ind_names = new string[nind];

	for (int i = 0; i < numComments; i++) getline(fin, line);
	stringstream ss;
	ss.str(line);
	for (int i = 0; i < 9; i++) ss >> line;
	for (int i = 0; i < nind; i++) ss >> data->ind_names[i];

	string gt, a1, a2;
	int i = 0;
	for (int locus = 0; locus < nloci; locus++)
	{
		getline(fin, line);
		if (EHOM->at(locus) <= MAF_HOM) {
			ss.clear();
			ss.str(line);
			for (int field = 0; field < 9; field++) ss >> gt;
			int index = 0;
			for (int ind = 0; ind < nind; ind++) {
				ss >> gt;
				if (gt.compare("./.") == 0 || gt.compare(".|.") == 0 || gt.compare(".") == 0) {
					data->data[index][i] = -9;
					index++;
					data->data[index][i] = -9;
					index++;
				}
				else {
					a1 = gt[0];
					a2 = gt[2];
					data->data[index][i] = atoi(a1.c_str());
					index++;
					data->data[index][i] = atoi(a2.c_str());
					index++;
				}
			}
			i++;
		}
	}
	fin.close();
	nloci = nkeep;
	delete EHOM;
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

