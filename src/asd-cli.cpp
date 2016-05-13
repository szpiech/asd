#include "asd-cli.h"

const string VERSION = "1.0.0";

const string ARG_CALC_IBS = "--ibs-all";
const bool DEFAULT_CALC_IBS = false;
const string HELP_CALC_IBS = "Set to output all IBS sharing matrices.";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "Basename for output files.";

const string ARG_FILENAME = "--stru";
const string DEFAULT_FILENAME = "__none";
const string HELP_FILENAME = "The input data filename (stru format), data must be biallelic and coded 0/1 (-9 for missing data or change with --stru-mis).";

const string ARG_TPED_FILENAME = "--tped";
const string DEFAULT_TPED_FILENAME = "__none";
const string HELP_TPED_FILENAME = "The input data filename (tped format), also requires a .tfam file (--tfam).  Data must be biallelic, but alleles can be arbitrarily coded.  Missing data is coded 0 (or change with --tped-mis).";

const string ARG_TFAM_FILENAME = "--tfam";
const string DEFAULT_TFAM_FILENAME = "__none";
const string HELP_TFAM_FILENAME = "The input data filename (tfam format), also requires a .tped file (--tped).";

const string ARG_SORT = "--id-col";
const int DEFAULT_SORT = 1;
const string HELP_SORT = "Column where individual IDs are. (stru files only)";

const string ARG_THREAD = "--threads";
const int DEFAULT_THREAD = 1;
const string HELP_THREAD = "Number of threads to spawn for faster calculation.";

const string ARG_PARTIAL = "--partial";
const bool DEFAULT_PARTIAL = false;
const string HELP_PARTIAL = "If set, outputs two r/2 x r/2 matrices. The first is the number of loci used for each pairwise comparison, and the second is the untransformed IBS matrix. In order to recover the allele sharing distance matrix divide the untransformed IBS counts by the number of loci in the appropriate cell, and then transform by 1-(ps).";

const string ARG_LOG = "--log";
const bool DEFAULT_LOG = false;
const string HELP_LOG = "Transforms allele sharing distance by -ln(ps) instead of 1-(ps).";

const string ARG_STRU_MISSING = "--missing-stru";
const string DEFAULT_STRU_MISSING = "-9";
const string HELP_STRU_MISSING = "For stru files, set the missing data value.";

const string ARG_TPED_MISSING = "--missing-tped";
const string DEFAULT_TPED_MISSING = "0";
const string HELP_TPED_MISSING = "For stru files, set the missing data value.";

const string ARG_BIALLELIC = "--biallelic";
const bool DEFAULT_BIALLELIC = false;
const string HELP_BIALLELIC = "Set for more efficient computations when all loci are biallelic.";

string getCommandLineString(int argc, char *argv[])
{
	string str = argv[0];
	for (int i = 1; i < argc; i++) {
		str += " " + string(argv[i]);
	}
	return str;
}

param_t *getCLI(int argc, char *argv[])
{
	param_t *params = new param_t;
	params->addFlag(ARG_THREAD, DEFAULT_THREAD, "", HELP_THREAD);
	params->addFlag(ARG_CALC_IBS, DEFAULT_CALC_IBS, "", HELP_CALC_IBS);
	params->addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
	params->addFlag(ARG_FILENAME, DEFAULT_FILENAME, "", HELP_FILENAME);
	params->addFlag(ARG_TPED_FILENAME, DEFAULT_TPED_FILENAME, "", HELP_TPED_FILENAME);
	params->addFlag(ARG_TFAM_FILENAME, DEFAULT_TFAM_FILENAME, "", HELP_TFAM_FILENAME);
	params->addFlag(ARG_SORT, DEFAULT_SORT, "", HELP_SORT);
	params->addFlag(ARG_PARTIAL, DEFAULT_PARTIAL, "", HELP_PARTIAL);
	params->addFlag(ARG_LOG, DEFAULT_LOG, "", HELP_LOG);
	params->addFlag(ARG_STRU_MISSING, DEFAULT_STRU_MISSING, "", HELP_STRU_MISSING);
	params->addFlag(ARG_TPED_MISSING, DEFAULT_TPED_MISSING, "", HELP_TPED_MISSING);
	params->addFlag(ARG_BIALLELIC, DEFAULT_BIALLELIC, "", HELP_BIALLELIC);
	
	params->setPreamble("asd v" + VERSION);

	if (!params->parseCommandLine(argc, argv))
	{
		delete params;
		return NULL;
	}
	return params;
}


bool check_int_gt_0(int n) {
	return (n > 0);
}
bool check_int_ge_0(int n) {
	return (n >= 0);
}
bool check_sort_ge_ndcols(int sort, int ndcols) {
	return (sort >= ndcols);
}
bool check_file_type(bool STRU, bool TPED, bool TFAM) {
	return ((STRU && !TPED && !TFAM) || (!STRU && TFAM && TPED));
}
