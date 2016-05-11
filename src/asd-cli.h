#ifndef __ASD_CLI_H__
#define __ASD_CLI_H__

#include <iostream>
#include <string>
#include "param_t.h"
#include "errlog.h"

using namespace std;

extern const string VERSION;

extern const string ARG_CHECK_FILE;
extern const bool DEFAULT_CHECK_FILE;
extern const string HELP_CHECK_FILE;

extern const string ARG_CHECK_FILE_DEEP;
extern const bool DEFAULT_CHECK_FILE_DEEP;
extern const string HELP_CHECK_FILE_DEEP;

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

extern const string ARG_NROWS;
extern const int DEFAULT_NROWS;
extern const string HELP_NROWS;

extern const string ARG_NCOLS;
extern const int DEFAULT_NCOLS;
extern const string HELP_NCOLS;

extern const string ARG_SORT;
extern const int DEFAULT_SORT;
extern const string HELP_SORT;

extern const string ARG_NDCOLS;
extern const int DEFAULT_NDCOLS;
extern const string HELP_NDCOLS;

extern const string ARG_NDROWS;
extern const int DEFAULT_NDROWS;
extern const string HELP_NDROWS;

extern const string ARG_THREAD;
extern const int DEFAULT_THREAD;
extern const string HELP_THREAD;

extern const string ARG_FULL;
extern const bool DEFAULT_FULL;
extern const string HELP_FULL;

extern const string ARG_FULL_LOG;
extern const bool DEFAULT_FULL_LOG;
extern const string HELP_FULL_LOG;

extern const string ARG_STRU_MISSING;
extern const int DEFAULT_STRU_MISSING;
extern const string HELP_STRU_MISSING;

extern const string ARG_TPED_MISSING;
extern const string DEFAULT_TPED_MISSING;
extern const string HELP_TPED_MISSING;

/*
extern const string ARG_CALC_ASD;
extern const bool DEFAULT_CALC_ASD;
extern const string HELP_CALC_ASD;

extern const string ARG_CALC_FST;
extern const bool DEFAULT_CALC_FST;
extern const string HELP_CALC_FST;
*/

param_t *getCLI(int argc, char *argv[]);
string getCommandLineString(int argc, char *argv[]);
bool check_int_gt_0(int n);
bool check_int_ge_0(int n);
bool check_sort_ge_ndcols(int sort, int ndcols);
bool check_print_full(bool PRINT_FULL, bool PRINT_FULL_LOG);
bool check_file_type(bool STRU, bool TPED, bool TFAM);
bool check_file_check(bool CHECK_FILE, bool CHECK_FILE_DEEP, bool STRU);
#endif