#include "asd-cli.h"

const string VERSION = "1.1.0a";

const string ARG_CALC_IBS = "--ibs";
const bool DEFAULT_CALC_IBS = false;
const string HELP_CALC_IBS = "Set to output all IBS sharing matrices too.";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "Basename for output files.";

const string ARG_FILENAME = "--stru";
const string DEFAULT_FILENAME = "__none";
const string HELP_FILENAME = "The input data filename (stru format).\
	\n\tChange missing data code with --missing-stru.\
	\n\tRequires two header rows: locus names, and locus positions.";

const string ARG_TPED_FILENAME = "--tped";
const string DEFAULT_TPED_FILENAME = "__none";
const string HELP_TPED_FILENAME = "The input data filename (tped format).\
	\n\tRequires a .tfam file (--tfam).\
	\n\tChange missing data code with --missing-tped.";

const string ARG_TFAM_FILENAME = "--tfam";
const string DEFAULT_TFAM_FILENAME = "__none";
const string HELP_TFAM_FILENAME = "The input data filename (tfam format).\
	\n\tRequires a .tped file (--tped).";

const string ARG_VCF_FILENAME = "--vcf";
const string DEFAULT_VCF_FILENAME = "__none";
const string HELP_VCF_FILENAME = "The input data filename (vcf format).";

/*
const string ARG_MAF = "--maf";
const double DEFAULT_MAF = 0.0;
const string HELP_MAF = "Filter sites with a MAF below this. This part may be slow for very large files.";
*/

const string ARG_SORT = "--id";
const int DEFAULT_SORT = 1;
const string HELP_SORT = "Column where individual IDs are. (stru files only)";

const string ARG_THREAD = "--threads";
const int DEFAULT_THREAD = 1;
const string HELP_THREAD = "Number of threads to spawn for faster calculation.";

const string ARG_PARTIAL = "--partial";
const bool DEFAULT_PARTIAL = false;
const string HELP_PARTIAL =
    "If set, outputs two nind x nind matrices.\
    \n\tThe first is the number of loci used for each pairwise\
    \n\tcomparison, and the second is the untransformed dist/IBS\
	\n\tmatrix. Useful for splitting up very large jobs.\
	\n\tCan combine with --combine.";

const string ARG_LOG = "--log";
const bool DEFAULT_LOG = false;
const string HELP_LOG = "Transforms allele sharing distance by -ln(ps) instead of 1-(ps).";

const string ARG_STRU_MISSING = "--missing-stru";
const string DEFAULT_STRU_MISSING = "-9";
const string HELP_STRU_MISSING = "For stru files, set the missing data value.";

const string ARG_TPED_MISSING = "--missing-tped";
const string DEFAULT_TPED_MISSING = "0";
const string HELP_TPED_MISSING = "For stru files, set the missing data value.";

const string ARG_MULTIALLELIC = "--multiallelic";
const bool DEFAULT_MULTIALLELIC = false;
const string HELP_MULTIALLELIC = "Set if there are more than two alleles at at least one locus.\
	\n\tComputations are less efficient with this flag set.";

const string ARG_COMBINE = "--combine";
const string DEFAULT_COMBINE = "__none";
const string HELP_COMBINE = "Combine several files\
	\n\tgenerated with --partial.";

const string ARG_LONG_FORMAT = "--long";
const bool DEFAULT_LONG_FORMAT = false;
const string HELP_LONG_FORMAT = "Instead of printing a matrix, print allele sharing distances one\
	\n\tper row. Formatted <ID1> <ID2> <allele sharing distance>.\
	\n\tNot compatible with --partial.";

const string ARG_IBS_LONG = "--long-ibs";
const bool DEFAULT_IBS_LONG = false;
const string HELP_IBS_LONG = "Instead of printing a matrix, print IBS calculations one\
	\n\tper row. Formatted <ID1> <ID2> <IBS0/1/2>.\
	\n\tNot compatible with --partial.";

const string ARG_GRM = "--grm";
const bool DEFAULT_GRM = false;
const string HELP_GRM = "Calculate the genomic relationship matrix.";

extern const string ARG_WEIGHTED_ASD = "--weighted";
extern const bool DEFAULT_WEIGHTED_ASD = false;
extern const string HELP_WEIGHTED_ASD = "Calculate weighted allele sharing similarity scores from Greenbaum et al. 2019.";


/*
const string ARG_KEEP_SITES_ID = "--keep-sites-id";
const string DEFAULT_KEEP_SITES_ID = "__none";
const string HELP_KEEP_SITES_ID = "A file containing a list of site IDs,\
\n\tone per line, to keep for computations.";

const string ARG_KEEP_SITES_POS = "--keep-sites-pos";
const string DEFAULT_KEEP_SITES_POS = "__none";
const string HELP_KEEP_SITES_POS = "A file containing a list of site positions,\
\n\tone per line, to keep for computations.";

const string ARG_KEEP_IND = "--keep-ind";
const string DEFAULT_KEEP_IND = "__none";
const string HELP_KEEP_IND = "A file containing a list of individual IDs,\
\n\tone per line, to keep for computations.";
*/

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
	params->addFlag(ARG_VCF_FILENAME, DEFAULT_VCF_FILENAME, "", HELP_VCF_FILENAME);
	//params->addFlag(ARG_MAF, DEFAULT_MAF, "", HELP_MAF);
	params->addFlag(ARG_SORT, DEFAULT_SORT, "", HELP_SORT);
	params->addFlag(ARG_PARTIAL, DEFAULT_PARTIAL, "", HELP_PARTIAL);
	params->addFlag(ARG_LOG, DEFAULT_LOG, "", HELP_LOG);
	params->addFlag(ARG_STRU_MISSING, DEFAULT_STRU_MISSING, "", HELP_STRU_MISSING);
	params->addFlag(ARG_TPED_MISSING, DEFAULT_TPED_MISSING, "", HELP_TPED_MISSING);
	params->addFlag(ARG_MULTIALLELIC, DEFAULT_MULTIALLELIC, "", HELP_MULTIALLELIC);
	params->addListFlag(ARG_COMBINE, DEFAULT_COMBINE, "", HELP_COMBINE);
	params->addFlag(ARG_LONG_FORMAT, DEFAULT_LONG_FORMAT, "", HELP_LONG_FORMAT);
	params->addFlag(ARG_IBS_LONG, DEFAULT_IBS_LONG, "", HELP_IBS_LONG);
	params->addFlag(ARG_GRM, DEFAULT_GRM, "", HELP_GRM);
	params->addFlag(ARG_WEIGHTED_ASD, DEFAULT_WEIGHTED_ASD, "", HELP_WEIGHTED_ASD);
	/*
	params->addFlag(ARG_KEEP_SITES_ID, DEFAULT_KEEP_SITES_ID, "", HELP_KEEP_SITES_ID);
	params->addFlag(ARG_KEEP_SITES_POS, DEFAULT_KEEP_SITES_POS, "", HELP_KEEP_SITES_POS);
	params->addFlag(ARG_KEEP_IND, DEFAULT_KEEP_IND, "", HELP_KEEP_IND);
	*/
	params->setPreamble("asd v" + VERSION);

	if (!params->parseCommandLine(argc, argv))
	{
		delete params;
		return NULL;
	}
	return params;
}

bool check_maf(double maf){
	return (maf >= 0 && maf <= 1);
}

bool check_int_gt_0(int n) {
	return (n > 0);
}
bool check_int_ge_0(int n) {
	return (n >= 0);
}
bool check_sort(int sort, int ndcols) {
	return (sort <= ndcols && sort >= 1);
}
bool check_file_type(bool STRU, bool TPED, bool TFAM, bool VCF) {
	return ((STRU && !TPED && !TFAM && !VCF) || (!STRU && TFAM && TPED && !VCF) || (!STRU && !TFAM && !TPED && VCF));
}
