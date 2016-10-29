#ifndef __ASD_CLI_H__
#define __ASD_CLI_H__

#include <iostream>
#include <string>
#include "param_t.h"
#include "errlog.h"

using namespace std;

extern const string VERSION;

extern const string ARG_CALC_IBS;
extern const bool DEFAULT_CALC_IBS;
extern const string HELP_CALC_IBS;

extern const string ARG_OUTFILE;
extern const string DEFAULT_OUTFILE;
extern const string HELP_OUTFILE;

extern const string ARG_FILENAME;
extern const string DEFAULT_FILENAME;
extern const string HELP_FILENAME;

extern const string ARG_TPED_FILENAME;
extern const string DEFAULT_TPED_FILENAME;
extern const string HELP_TPED_FILENAME;

extern const string ARG_TFAM_FILENAME;
extern const string DEFAULT_TFAM_FILENAME;
extern const string HELP_TFAM_FILENAME;

extern const string ARG_VCF_FILENAME;
extern const string DEFAULT_VCF_FILENAME;
extern const string HELP_VCF_FILENAME;

extern const string ARG_MAF;
extern const double DEFAULT_MAF;
extern const string HELP_MAF;

extern const string ARG_SORT;
extern const int DEFAULT_SORT;
extern const string HELP_SORT;

extern const string ARG_THREAD;
extern const int DEFAULT_THREAD;
extern const string HELP_THREAD;

extern const string ARG_PARTIAL;
extern const bool DEFAULT_PARTIAL;
extern const string HELP_PARTIAL;

extern const string ARG_LOG;
extern const bool DEFAULT_LOG;
extern const string HELP_LOG;

extern const string ARG_STRU_MISSING;
extern const string DEFAULT_STRU_MISSING;
extern const string HELP_STRU_MISSING;

extern const string ARG_TPED_MISSING;
extern const string DEFAULT_TPED_MISSING;
extern const string HELP_TPED_MISSING;

extern const string ARG_BIALLELIC;
extern const bool DEFAULT_BIALLELIC;
extern const string HELP_BIALLELIC;

extern const string ARG_COMBINE;
extern const string DEFAULT_COMBINE;
extern const string HELP_COMBINE;

extern const string ARG_LONG_FORMAT;
extern const bool DEFAULT_LONG_FORMAT;
extern const string HELP_LONG_FORMAT;

extern const string ARG_IBS_LONG;
extern const bool DEFAULT_IBS_LONG;
extern const string HELP_IBS_LONG;

param_t *getCLI(int argc, char *argv[]);
string getCommandLineString(int argc, char *argv[]);
bool check_maf(double maf);
bool check_int_gt_0(int n);
bool check_int_ge_0(int n);
bool check_sort(int sort, int ndcols);
bool check_file_type(bool STRU, bool TPED, bool TFAM, bool VCF);

#endif