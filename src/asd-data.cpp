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

map<string,bool>* readSubsetFile(string infile){
	igzstream fin;
	fin.open(infile.c_str());

	if(fin.fail()){
		LOG.err("ERROR: Failed to open",infile);
		throw 0;
	}

	string id;
	map<string,bool>* subset = new map<string,bool>;
	
	while(fin >> id) {
		subset->operator[](id) = true;
	}

	fin.close();
	return subset;
}

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

bool finalize_calculations(int nind, int ncols, bool CALC_ALL_IBS, bool GRM, bool ASD, bool WEIGHTED_ASD) {
	for (int i = 0; i < nind; i++)
	{
		for (int j = i; j < nind; j++)
		{
			if (i == j && ASD)
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

void write_dist_matrix(string outfile, int nind, int ncols, string *ind_names, bool PRINT_PARTIAL, bool PRINT_LOG, bool PRINT_LONG, bool GRM, bool ASD, bool WEIGHTED_ASD) {
	ofstream out;
	string type;

	if(GRM){
		type = "grm";
	}
	else if (ASD){
		type = "dist";
	}
	else if (WEIGHTED_ASD){
		type = "wsim";
	}

	if (!PRINT_PARTIAL) {
		if(GRM){
			outfile += ".grm";
		}
		else if (ASD){
			outfile += ".asd.dist";
		}
		else if (WEIGHTED_ASD){
			outfile += ".wsim";
		}
	}
	else {
		if(GRM){
			outfile += ".grm.partial";
		}
		else if (ASD){
			outfile += ".asd.dist.partial";
		}
		else if (WEIGHTED_ASD){
			outfile += ".wsim.partial";
		}
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
				if (!PRINT_PARTIAL && PRINT_LOG && ASD){
					out << 0 - log(double(DIST_MAT[i][j]) /
					               (double(ncols) + double(NUM_LOCI[i][j])))
					    << " ";
				}
				else if (!PRINT_PARTIAL && ASD){
					out << 1 - (double(DIST_MAT[i][j]) /
					            (double(ncols) + double(NUM_LOCI[i][j])))
					    << " ";
				}
				else if (!PRINT_PARTIAL && (GRM || WEIGHTED_ASD)){
					if(i == j && WEIGHTED_ASD){
						out << 1 << " ";
					}
					else{
						out << (double(DIST_MAT[i][j]) /
						            (double(ncols) + double(NUM_LOCI[i][j])))
						    << " ";
					}
				}
				else out << DIST_MAT[i][j] << " ";
			}
			out << endl;
		}
	}
	else if (PRINT_LONG) {
		for (int i = 0; i < nind; i++) {
			for (int j = i+1; j < nind ; j++) {
				out << ind_names[i] << " " << ind_names[j] << " ";
				if (PRINT_LOG && ASD) {
					out << 0 - log(double(DIST_MAT[i][j]) /
					               (double(ncols) + double(NUM_LOCI[i][j])))
					    << endl;
				}
				else if(GRM || WEIGHTED_ASD){
					out << (double(DIST_MAT[i][j]) /
					            (double(ncols) + double(NUM_LOCI[i][j])))
					    << endl;
				}
				else if (ASD){
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

			if (ndrows < 2) {
				LOG.err("ERROR: Non-data rows must be >= 2. Found", ndrows);
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

	fin.close();

	nkeep = ncols;
	delete [] allele2count;
	return NULL;
}

/*
structure_data *readData_stru(string infile, int sort, int &nrows, int &ncols, double MAF, string STRU_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS) {
	int ndcols = -1;
	int ndrows = -1;
	nrows = -1;
	ncols = -1;
	string line;
	igzstream fin;
	int nkeep = 0;

	vector<double> *AF0 = biallelic_maf_stru(infile, sort, nrows, ncols, ndrows, ndcols, nkeep, MAF, STRU_MISSING);

	int expectednindFilter;

	if(KEEP_IND) expectednindFilter = keepIND->size();
	else expectednindFilter = nrows/2;	

	fin.open(infile.c_str());

	bool *keep_locus;
	if(KEEP_SNP || KEEP_POS) keep_locus = new bool[ncols];

	string locus_name;
	int totalLoci = 0;
	//Get locus names
	if(KEEP_SNP){
		for (int locus = 0; locus < ncols; locus++){
			fin >> locus_name;
			if(keepSNP->count(locus_name) > 0){
				keep_locus[locus] = true;
				totalLoci++;
			}
			else{
				keep_locus[locus] = false;
			}
		}
	}
	else{
		getline(fin, line);
	}
	//Get locus positions
	if(KEEP_POS){
		for (int locus = 0; locus < ncols; locus++){
			fin >> locus_name;
			if(keepPOS->count(locus_name) > 0){
				keep_locus[locus] = true;
				totalLoci++;
			}
			else{
				keep_locus[locus] = false;
			}
		}
	}
	else{
		getline(fin, line);
	}
	//toss out other non-data rows
	for (int i = 2; i < ndrows; i++) getline(fin, line);

	if(totalLoci == 0) totalLoci = ncols;

	structure_data *data = new structure_data;
	int nind = nrows / 2;
	data->nindAllocated = expectednindFilter;
	data->data = new short*[expectednindFilter];
	for (int i = 0; i < expectednindFilter; i++) {
		data->data[i] = new short[totalLoci];
	}
	data->ind_names = new string[expectednindFilter];

	map<string, short> *allele2code = new map<string, short>[totalLoci];
	short *lastAlleleCode = new short[totalLoci];
	for (int i = 0; i < totalLoci; i++) {
		allele2code[i][STRU_MISSING] = -9;
		lastAlleleCode[i] = -1;
	}

	int index, actualInd = 0;
	string field, indID, allele1, allele2, line1, line2;
	short genotypeCode;
	for (int ind = 0; ind < nind; ind++) {
		stringstream sin1, sin2;
		getline(fin, line1);
		sin1 << line1;
		getline(fin, line2);
		sin2 << line2;

		for (int i = 1; i <= ndcols; i++) {
			sin1 >> field;
			sin2 >> field;
			if (i == sort) indID = field;
		}

		if(KEEP_IND){
			if(keepIND->count(indID) == 0) continue;
		}

		data->ind_names[actualInd] = indID;

		int i = 0;
		for (int locus = 0; locus < ncols; locus++) {

			sin1 >> allele1;
			sin2 >> allele2;

			if(KEEP_POS || KEEP_SNP){
				if(!keep_locus[locus]) continue;
			}

			//if (AF0->at(locus) >= MAF && AF0->at(locus) <= 1 - MAF) {
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
					data->data[actualInd][i] = allele2code[i][STRU_MISSING];
				}
				else {
					data->data[actualInd][i] = genotypeCode;
				}
				i++;
			//}
		}
		actualInd++;
	}
	fin.close();
	ncols = totalLoci;
	nrows = 2*actualInd;

	data->nloci = totalLoci;
	data->nind = actualInd;
	if(KEEP_SNP || KEEP_POS) delete [] keep_locus;
	return data;
}

structure_data *readData_stru2(string infile, int sort, int &nrows, int &ncols, double MAF, string STRU_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS) {
	string line;
	int ndcols = -1;
	int ndrows = -1;
	nrows = -1;
	ncols = -1;
	int nkeep = 0;
	//double MAF_HOM = biallelic_ehom(MAF);
	vector<double>* EHOM = multiallelic_maf_stru(infile, sort, nrows, ncols, ndrows, ndcols, nkeep, MAF, STRU_MISSING);

	int expectednindFilter;

	if(KEEP_IND) expectednindFilter = keepIND->size();
	else expectednindFilter = nrows/2;	

	igzstream fin;
	fin.open(infile.c_str());

	bool *keep_locus;
	if(KEEP_SNP || KEEP_POS) keep_locus = new bool[ncols];

	string locus_name;
	int totalLoci = 0;
	//Get locus names
	if(KEEP_SNP){
		for (int locus = 0; locus < ncols; locus++){
			fin >> locus_name;
			if(keepSNP->count(locus_name) > 0){
				keep_locus[locus] = true;
				totalLoci++;
			}
			else{
				keep_locus[locus] = false;
			}
		}
	}
	else{
		getline(fin, line);
	}
	//Get locus positions
	if(KEEP_POS){
		for (int locus = 0; locus < ncols; locus++){
			fin >> locus_name;
			if(keepPOS->count(locus_name) > 0){
				keep_locus[locus] = true;
				totalLoci++;
			}
			else{
				keep_locus[locus] = false;
			}
		}
	}
	else{
		getline(fin, line);
	}
	//toss out other non-data rows
	for (int i = 2; i < ndrows; i++) getline(fin, line);

	if(totalLoci == 0) totalLoci = ncols;


	structure_data *data = new structure_data;
	int nind = nrows / 2;
	data->nindAllocated = expectednindFilter;
	data->data = new short*[2*expectednindFilter];
	for (int i = 0; i < 2*expectednindFilter; i++) {
		data->data[i] = new short[totalLoci];
	}
	data->ind_names = new string[expectednindFilter];


	string field, indID, allele;
	map<string, short> *allele2code = new map<string, short>[totalLoci];
	short *lastAlleleCode = new short[totalLoci];
	for (int i = 0; i < totalLoci; i++) {
		allele2code[i][STRU_MISSING] = -9;
		lastAlleleCode[i] = -1;
	}
	int index, actualRow = 0;
	for (int row = 0; row < nrows; row++) {

		for (int i = 1; i <= ndcols; i++) {
			fin >> field;
			if (i == sort) indID = field;
		}

		if(KEEP_IND){
			if(keepIND->count(indID) == 0){
				getline(fin,field);
				continue;
			}
		}

		if (actualRow % 2 == 0) data->ind_names[actualRow / 2] = indID;

		index = 0;
		for (int locus = 0; locus < ncols; locus++) {
			fin >> allele;

			if(KEEP_POS || KEEP_SNP){
				if(!keep_locus[locus]) continue;
			}

			if (allele2code[index].count(allele) == 0) {
				lastAlleleCode[index]++;
				allele2code[index][allele] = lastAlleleCode[index];
				data->data[actualRow][index] = lastAlleleCode[index];
			}
			else {
				data->data[actualRow][index] = allele2code[index][allele];
			}
			index++;
		}
		actualRow++;
	}
	fin.close();
	ncols = totalLoci;
	nrows = actualRow;

	data->nloci = totalLoci;
	data->nind = actualRow/2;
	if(KEEP_SNP || KEEP_POS) delete [] keep_locus;

	return data;
}
*/

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

	if (!check_sort(sort, ndcols)) {
		err = true;
		LOG.err("ERROR: Individual ID column must be in [ 1,", ndcols, false);
		LOG.err(" ].");
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

	if (!check_sort(sort, ndcols)) {
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
	//stringstream ss;
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

			if (ndrows < 2) {
				LOG.err("ERROR: Non-data rows must be >= 2. Found", ndrows);
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
		totalRows++;
	}
	nrows = totalRows - ndrows;

	if (!check_int_gt_0(nrows)) {
		LOG.err("ERROR: Number of samples must be > 0. Found", nrows);
		throw 0;
	}

	fin.close();

	nkeep = nrows;
	return NULL;
}


/*
structure_data *readData_tped_tfam2(string tped_filename, string tfam_filename, int &nrow, int &nloci, double MAF, string TPED_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS) {
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;
	
	famin.open(tfam_filename.c_str());
	if (famin.fail())
	{
		LOG.err("ERROR: Coult not open", tfam_filename);
		throw - 1;
	}

	vector<string> names;
	vector<bool> keep_ind;
	nrow = 0;
	int actualInd = 0;
	while(getline(famin, junk)){
		stringstream ss;
		ss.str(junk);
		ss >> junk;
		ss >> junk;
		names.push_back(junk);
		nrow+=2;
		if(KEEP_IND){
			if(keepIND->count(junk) > 0){
				keep_ind.push_back(true);
				actualInd++;
			}
			else{
				keep_ind.push_back(false);
			}
		}
	}

	if(actualInd == 0) actualInd = nrow/2;

	structure_data *data = new structure_data;
	data->nind = actualInd;
	data->data = new short*[2*actualInd];

	pedin.open(tped_filename.c_str());
	if (pedin.fail())
	{
		LOG.err("ERROR: Coult not open", tped_filename);
		throw - 1;
	}

	int actualLoci = 0;
	string snpid, pos;
	vector<bool> keep_locus;
	while(getline(pedin,junk)){
		if(KEEP_SNP || KEEP_POS){
			stringstream ss;
			ss.str(junk);
			ss >> junk >> snpid >> junk >> pos;
			if(KEEP_SNP){
				if(keepSNP->count(snpid) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
			else if (KEEP_POS){
				if(keepPOS->count(pos) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
		}
		nloci++;
	}
	pedin.close();
	pedin.clear();

	if(actualLoci == 0) actualLoci = nloci;

	data->nloci = actualLoci;
	for (int i = 0; i < 2*actualInd; i++) data->data[i] = new short[actualLoci];
	data->ind_names = new string[actualInd];

	pedin.open(tped_filename.c_str());

	map<string, short> allele2code;
	short lastAlleleCode = -1;

	string allele;
	int index = 0;
	for (int locus = 0; locus < nloci; locus++) {
		getline(pedin, junk);

		if(KEEP_POS || KEEP_SNP){
			if(!keep_locus[locus]) continue;
		}

		stringstream ss;
		ss.str(junk);
		allele2code.clear();
		allele2code[TPED_MISSING] = -9;
		lastAlleleCode = -1;

		ss >> allele;
		ss >> allele;
		ss >> allele;
		ss >> allele;

		int currI = 0;
		for (int i = 0; i < nrow; i++) {
			ss >> allele;
			if(KEEP_IND){
				if(!keep_ind[i/2]) continue;
			}

			data->ind_names[currI/2] = names[i/2];

			if (allele2code.count(allele) == 0) {
				lastAlleleCode++;
				allele2code[allele] = lastAlleleCode;
				data->data[currI][index] = lastAlleleCode;
			}
			else {
				data->data[currI][index] = allele2code[allele];
			}
			currI++;
		}
		index++;
	}

	nloci = actualLoci;
	nrow = 2*actualInd;
	pedin.close();
	return data;
}

structure_data *readData_tped_tfam(string tped_filename, string tfam_filename, int &nrow, int &nloci, double MAF, string TPED_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS) {
	string junk;
	nrow = 0;
	nloci = 0;
	igzstream famin, pedin;
	//int nkeep = 0;
	//vector<double> *AF0 = biallelic_maf_tped(tped_filename, tfam_filename, nrow, nloci, nkeep, MAF, TPED_MISSING);
	
	famin.open(tfam_filename.c_str());
	if (famin.fail())
	{
		LOG.err("ERROR: Coult not open", tfam_filename);
		throw - 1;
	}

	vector<string> names;
	vector<bool> keep_ind;
	nrow = 0;
	int actualInd = 0;
	while(getline(famin, junk)){
		stringstream ss;
		ss.str(junk);
		ss >> junk;
		ss >> junk;
		names.push_back(junk);
		nrow+=2;
		if(KEEP_IND){
			if(keepIND->count(junk) > 0){
				keep_ind.push_back(true);
				actualInd++;
			}
			else{
				keep_ind.push_back(false);
			}
		}
	}

	if(actualInd == 0) actualInd = nrow/2;

	structure_data *data = new structure_data;
	data->nind = actualInd;
	data->data = new short*[actualInd];

	pedin.open(tped_filename.c_str());
	if (pedin.fail()){
		LOG.err("ERROR: Coult not open", tped_filename);
		throw - 1;
	}

	int actualLoci = 0;
	string snpid, pos;
	vector<bool> keep_locus;
	while(getline(pedin,junk)){
		if(KEEP_SNP || KEEP_POS){
			stringstream ss;
			ss.str(junk);
			ss >> junk >> snpid >> junk >> pos;
			if(KEEP_SNP){
				if(keepSNP->count(snpid) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
			else if (KEEP_POS){
				if(keepPOS->count(pos) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
		}
		nloci++;
	}
	pedin.close();
	pedin.clear();

	if(actualLoci == 0) actualLoci = nloci;

	data->nloci = actualLoci;
	for (int i = 0; i < actualInd; i++) data->data[i] = new short[actualLoci];
	data->ind_names = new string[actualInd];

	pedin.open(tped_filename.c_str());

	map<string, short> allele2code;
	short lastAlleleCode = -1;

	string allele1, allele2;
	short genotypeCode;
	int index = 0;
	for (int locus = 0; locus < nloci; locus++){
		getline(pedin, junk);
			
		if(KEEP_POS || KEEP_SNP){
			if(!keep_locus[locus]) continue;
		}

		stringstream ss;
		ss.str(junk);
		allele2code.clear();
		allele2code[TPED_MISSING] = -9;
		lastAlleleCode = -1;

		ss >> allele1;
		ss >> allele1;
		ss >> allele1;
		ss >> allele1;

		int currInd = 0;
		for (int ind = 0; ind < nrow/2; ind++){
			ss >> allele1;
			ss >> allele2;

			if(KEEP_IND){
				if(!keep_ind[ind]) continue;
			}

			data->ind_names[currInd] = names[ind];

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
				data->data[currInd][index] = allele2code[TPED_MISSING];
			}
			else {
				data->data[currInd][index] = genotypeCode;
			}

			currInd++;
		}
		index++;
	}

	nloci = actualLoci;
	nrow = 2*actualInd;
	pedin.close();
	//delete AF0;
	return data;
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

/*
structure_data *readData_vcf(string vcf_filename, int &nrow, int &nloci, double MAF,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS) {
	int numComments = 0;
	string line, junk;

	igzstream fin;

	fin.open(vcf_filename.c_str());
	if (fin.fail()){
		LOG.err("ERROR: Coult not open", vcf_filename);
		throw - 1;
	}

	while(getline(fin, line)){
		numComments++;
		if(line[0] == '#' && line[1] == 'C') break;
	}

	stringstream ss;
	ss.str(line);
	for (int i = 0; i < 9; i++) ss >> line;

	int actualInd = 0;
	vector<string> names;
	vector<bool> keep_ind;
	while(ss >> line){
		names.push_back(line);
		if(KEEP_IND){
			if(keepIND->count(line)>0){
				keep_ind.push_back(true);
				actualInd++;
			}
			else keep_ind.push_back(false);
		}
		nrow++;
	}

	nloci = 0;
	int actualLoci = 0;
	vector<bool> keep_locus;
	string pos, id;
	while(getline(fin, line)){
		if(KEEP_SNP || KEEP_POS){
			//stringstream ss;
			ss.clear();
			ss.str(line);
			ss >> junk >> pos >> id;
			if(KEEP_SNP){
				if(keepSNP->count(id) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
			else if (KEEP_POS){
				if(keepPOS->count(pos) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
		}
		nloci++;
	}
	fin.close();
	fin.clear();

	fin.open(vcf_filename.c_str());

	structure_data *data = new structure_data;
	data->nloci = actualLoci;
	data->nind = actualInd;
	data->data = new short*[actualInd];
	for (int i = 0; i < actualInd; i++) data->data[i] = new short[actualLoci];
	data->ind_names = new string[actualInd];

	for (int i = 0; i < numComments; i++) getline(fin, line);
	
	string gt, a1, a2;
	int i = 0;
	for (int locus = 0; locus < nloci; locus++){
		getline(fin, line);

		if(KEEP_POS || KEEP_SNP){
			if(!keep_locus[locus]) continue;
		}

		ss.clear();
		ss.str(line);
		for (int field = 0; field < 9; field++) ss >> gt;
		int currInd = 0;
		for (int ind = 0; ind < nrow; ind++) {
			ss >> gt;

			if(KEEP_IND){
				if(!keep_ind[ind]) continue;
			}

			data->ind_names[currInd] = names[ind];

			if (gt.compare("./.") == 0 || gt.compare(".|.") == 0 || gt.compare(".") == 0) {
				data->data[currInd][i] = -9;
				currInd++;
				continue;
			}
			else {
				a1 = gt[0];
				a2 = gt[2];

				if( (a1.compare("0") != 0 && a1.compare("1") != 0) || (a2.compare("0") != 0 && a2.compare("1") != 0) ){
					LOG.err("ERORR: --biallelic flag set, but found more than 2 alleles at locus", locus + 1);
					throw - 1;
				}
			}
			data->data[currInd][i] = atoi(a1.c_str()) + atoi(a2.c_str());
			currInd++;
		}
		i++;
	}

	nloci = actualLoci;
	nrow = 2*actualInd;
	fin.close();
	return data;
}

structure_data *readData_vcf2(string vcf_filename, int &nrow, int &nloci, double MAF,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS) {
	int numComments = 0;
	string line, junk;
	double MAF_HOM = biallelic_ehom(MAF);

	igzstream fin;

	fin.open(vcf_filename.c_str());
	if (fin.fail()){
		LOG.err("ERROR: Coult not open", vcf_filename);
		throw - 1;
	}

	while(getline(fin, line)){
		numComments++;
		if(line[0] == '#' && line[1] == 'C') break;
	}

	stringstream ss;
	ss.str(line);
	for (int i = 0; i < 9; i++) ss >> line;

	int actualInd = 0;
	vector<string> names;
	vector<bool> keep_ind;
	while(ss >> line){
		names.push_back(line);
		if(KEEP_IND){
			if(keepIND->count(line)>0){
				keep_ind.push_back(true);
				actualInd++;
			}
			else keep_ind.push_back(false);
		}
		nrow++;
	}

	nloci = 0;
	int actualLoci = 0;
	vector<bool> keep_locus;
	string pos, id;
	while(getline(fin, line)){
		if(KEEP_SNP || KEEP_POS){
			//stringstream ss;
			ss.clear();
			ss.str(line);
			ss >> junk >> pos >> id;
			if(KEEP_SNP){
				if(keepSNP->count(id) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
			else if (KEEP_POS){
				if(keepPOS->count(pos) > 0){
					actualLoci++;
					keep_locus.push_back(true);
				}
				else keep_locus.push_back(false);
			}
		}
		nloci++;
	}
	fin.close();
	fin.clear();

	fin.open(vcf_filename.c_str());
	structure_data *data = new structure_data;
	data->nloci = actualLoci;
	data->nind = actualInd;
	data->data = new short*[2*actualInd];
	for (int i = 0; i < 2*actualInd; i++) data->data[i] = new short[actualLoci];
	data->ind_names = new string[actualInd];

	for (int i = 0; i < numComments; i++) getline(fin, line);

	string gt, a1, a2;
	int i = 0;
	for (int locus = 0; locus < nloci; locus++)
	{
		getline(fin, line);

		if(KEEP_POS || KEEP_SNP){
			if(!keep_locus[locus]) continue;
		}

		ss.clear();
		ss.str(line);
		for (int field = 0; field < 9; field++) ss >> gt;
		int index = 0;
		int currInd = 0;
		for (int ind = 0; ind < nrow; ind++) {
			ss >> gt;

			if(KEEP_IND){
				if(!keep_ind[ind]) continue;
			}

			data->ind_names[currInd] = names[ind];

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
			currInd++;
		}
		i++;
	}

	nloci = actualLoci;
	nrow = 2*actualInd;

	fin.close();
	return data;
}
*/

structure_data *readData_vcf(string vcf_filename, int &nrow, int &nloci) {
	int numComments = 0;
	string line, junk;

	igzstream fin;

	fin.open(vcf_filename.c_str());
	if (fin.fail()){
		LOG.err("ERROR: Coult not open", vcf_filename);
		throw - 1;
	}

	while(getline(fin, line)){
		numComments++;
		if(line[0] == '#' && line[1] == 'C') break;
	}

	stringstream ss;
	ss.str(line);
	for (int i = 0; i < 9; i++) ss >> line;

	//int actualInd = 0;
	vector<string> names;
	//vector<bool> keep_ind;
	while(ss >> line){
		names.push_back(line);
		nrow++;
	}

	nloci = 0;
	//int actualLoci = 0;
	vector<bool> keep_locus;
	string pos, id;
	while(getline(fin, line)){
		nloci++;
	}
	fin.close();
	fin.clear();

	fin.open(vcf_filename.c_str());

	structure_data *data = new structure_data;
	data->nloci = nloci;
	data->nind = nrow;
	data->data = new short*[nrow];
	for (int i = 0; i < nrow; i++) data->data[i] = new short[nloci];
	data->ind_names = new string[nrow];

	for (int i = 0; i < numComments; i++) getline(fin, line);
	
	string gt, a1, a2;
	for (int locus = 0; locus < nloci; locus++){
		getline(fin, line);

		ss.clear();
		ss.str(line);
		for (int field = 0; field < 9; field++) ss >> gt;
	
		for (int ind = 0; ind < nrow; ind++) {
			ss >> gt;
			data->ind_names[ind] = names[ind];

			if (gt.compare("./.") == 0 || gt.compare(".|.") == 0 || gt.compare(".") == 0) {
				data->data[ind][locus] = -9;
				continue;
			}
			else {
				a1 = gt[0];
				a2 = gt[2];

				if( (a1.compare("0") != 0 && a1.compare("1") != 0) || (a2.compare("0") != 0 && a2.compare("1") != 0) ){
					LOG.err("ERORR: --biallelic flag set, but found more than 2 alleles at locus", locus + 1);
					throw - 1;
				}
			}
			data->data[ind][locus] = atoi(a1.c_str()) + atoi(a2.c_str());
		}
	}

	nrow *= 2;
	fin.close();
	return data;
}

structure_data *readData_vcf2(string vcf_filename, int &nrow, int &nloci) {
	int numComments = 0;
	string line, junk;

	igzstream fin;

	fin.open(vcf_filename.c_str());
	if (fin.fail()){
		LOG.err("ERROR: Coult not open", vcf_filename);
		throw - 1;
	}

	while(getline(fin, line)){
		numComments++;
		if(line[0] == '#' && line[1] == 'C') break;
	}

	stringstream ss;
	ss.str(line);
	for (int i = 0; i < 9; i++) ss >> line;

	vector<string> names;
	while(ss >> line){
		names.push_back(line);
		nrow++;
	}

	nloci = 0;
	string pos, id;
	while(getline(fin, line)){
		nloci++;
	}
	fin.close();
	fin.clear();

	fin.open(vcf_filename.c_str());
	structure_data *data = new structure_data;
	data->nloci = nloci;
	data->nind = nrow;
	data->data = new short*[2*nrow];
	for (int i = 0; i < 2*nrow; i++) data->data[i] = new short[nloci];
	data->ind_names = new string[nrow];

	for (int i = 0; i < numComments; i++) getline(fin, line);

	string gt, a1, a2;

	for (int locus = 0; locus < nloci; locus++)
	{
		getline(fin, line);

		ss.clear();
		ss.str(line);
		for (int field = 0; field < 9; field++) ss >> gt;
		int index = 0;

		for (int ind = 0; ind < nrow; ind++) {
			ss >> gt;

			data->ind_names[ind] = names[ind];

			if (gt.compare("./.") == 0 || gt.compare(".|.") == 0 || gt.compare(".") == 0) {
				data->data[index][locus] = -9;
				index++;
				data->data[index][locus] = -9;
				index++;
			}
			else {
				a1 = gt[0];
				a2 = gt[2];
				data->data[index][locus] = atoi(a1.c_str());
				index++;
				data->data[index][locus] = atoi(a2.c_str());
				index++;
			}
		}
	}

	fin.close();
	nrow *=2;
	return data;
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
	short int *a = new short int[fields];
	for (int i = 0; i < fields; i++){
		fin >> a[i];
	}
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

