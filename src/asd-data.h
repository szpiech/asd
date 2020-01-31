#ifndef __ASD_DATA_H__
#define __ASD_DATA_H__

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <pthread.h>
#include <map>
#include "gzstream.h"
#include "param_t.h"
#include "asd-cli.h"
#include "pbar.h"
#include "errlog.h"

extern const char DEL;
extern const string EMPTY_STRING;

extern double **DIST_MAT;
extern int **NUM_LOCI;
extern int **IBS_0_MAT;
extern int **IBS_1_MAT;
extern int **IBS_2_MAT;
extern pthread_mutex_t mutex_dist_mat;
extern pthread_mutex_t mutex_loci_mat;
extern pthread_mutex_t mutex_ibs_0;
extern pthread_mutex_t mutex_ibs_1;
extern pthread_mutex_t mutex_ibs_2;

typedef struct
{
	short **data;
	string *locus_names;
	string *ind_names;
	int nloci;
	int nind;
	int nindAllocated;
} structure_data;

typedef struct
{
	int first_index;
	int last_index;
	structure_data *stru_data;
	int missing;
	bool CALC_ALL_IBS;
	Bar *bar;
	int threads;
} work_order_t;
/*
typedef struct
{
	bool PRINT_FULL;
	bool PRINT_FULL_LOG;
	ostream *out;
	string *ind_names;
	int ncols;
	string type;
	int nind;
} output_order_t;
*/
/*
typedef struct
{
	double **data;
	string *pop_names;
	int *nind;
	int npop;
	int nloci;
} population_data;
*/

map<string,bool>* readSubsetFile(string infile); 

void combine_partial_files(param_t *params);

/*
structure_data *readData_stru(string infile, int sort, int &nrows, int &ncols, double MAF, string STRU_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS);
structure_data *readData_stru2(string infile, int sort, int &nrows, int &ncols, double MAF, string STRU_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS);
*/

structure_data *readData_stru(string infile, int sort, int &nrows, int &ncols, string STRU_MISSING);
structure_data *readData_stru2(string infile, int sort, int &nrows, int &ncols, string STRU_MISSING);

vector<double>* multiallelic_maf_stru(string infile, int sort, int &nrow, int &ncols, int &ndrows, int &ndcols, int &nkeep, double MAF, string STRU_MISSING);

/*
structure_data *readData_tped_tfam(string tped_filename, string tfam_filename, int &nrow, int &nloci, double MAF, string TPED_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS);
structure_data *readData_tped_tfam2(string tped_filename, string tfam_filename, int &nrow, int &nloci, double MAF, string TPED_MISSING,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS);
*/
structure_data *readData_tped_tfam(string tped_filename, string tfam_filename, int &nrow, int &nloci, string TPED_MISSING);
structure_data *readData_tped_tfam2(string tped_filename, string tfam_filename, int &nrow, int &nloci, string TPED_MISSING);

vector<double>* biallelic_maf_tped(string tped_filename, string tfam_filename, int &nrow, int &nloci, int &nkeep, double MAF, string TPED_MISSING);
vector<double>* multiallelic_maf_tped(string tped_filename, string tfam_filename, int &nrow, int &nloci, int &nkeep, double MAF, string TPED_MISSING);

/*
structure_data *readData_vcf(string vcf_filename, int &nrow, int &nloci, double MAF,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS);
structure_data *readData_vcf2(string vcf_filename, int &nrow, int &nloci, double MAF,
	map<string, bool> *keepIND, map<string, bool> *keepSNP, map<string, bool> *keepPOS, bool KEEP_IND, bool KEEP_SNP, bool KEEP_POS);
*/

structure_data *readData_vcf(string vcf_filename, int &nrow, int &nloci);
structure_data *readData_vcf2(string vcf_filename, int &nrow, int &nloci);


vector<double>* biallelic_maf_vcf(string vcf_filename, int &nrow, int &nloci, int &nkeep, int &numComments, double MAF);
vector<double>* multiallelic_maf_vcf(string vcf_filename, int &nrow, int &nloci, int &nkeep, int &numComments, double MAF);

double biallelic_ehom(double MAF);
/*
void readData_ind_asd(igzstream &fin, structure_data &data,
                      int sort, int ndcols, int ndrows, int nrows, int ncols, string STRU_MISSING);
*/

/*
void readData_ind_asd_tped_tfam(string tped_filename, string tfam_filename, structure_data &data,
                                int &nrow, int &nloci, string TPED_MISSING);
*/
int countFields(string junk);
void readData_pop_freq(igzstream &fin, structure_data &data,
                       int sort, int ndcols, int ndrows,
                       int nrows, int ncols);

short int *split_int(igzstream &fin, int fields);
string *split_str_str(int &size, const char *s, char c);

void output(void *order);
int search(string *s, int size, string key);
int put(string *s, int size, string key);

bool init_storage(int nind, bool CALC_ALL_IBS);
bool finalize_calculations(int nind, int ncols, bool CALC_ALL_IBS, bool GRM);

void write_ibs_matrices(string outfile, int nind, int ncols, string *ind_names, bool PRINT_PARTIAL, bool PRINT_LONG_IBS);
void write_dist_matrix(string outfile, int nind, int ncols, string *ind_names, bool PRINT_PARTIAL, bool PRINT_LOG, bool PRINT_LONG, bool GRM);

#endif