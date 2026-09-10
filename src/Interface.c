
// File: Interface.c
// Date: 18 February 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines LODESTAR configuration from command line interface and
//              output formatting of LODESTAR results.

#include "Interface.h"
#include "ketopt.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h> // Check if file exits.

// Make sure user defined arguments are valid.
// Accepts:
//  LodestarConfig_t* lodestarConfig -> The configured parameters.
// Returns: int, 0 if all parameters are valid. -1 if user supplied an invalid value.
int check_configuration(LodestarConfig_t* lodestarConfig) {
    // Check num reps and sample size for bootstrap.
    //  0 is legal for both and means "no bootstrap" / "use every block", so the old
    //  message ("must be an integer greater than 0") described the wrong condition.
    if (lodestarConfig -> numReps < 0) {
        fprintf(stderr, "-r %d must be INT >= 0 (0 for no bootstrap).\n", lodestarConfig -> numReps); 
        return -1;
    }
    if (lodestarConfig -> sampleSize < 0) {
        fprintf(stderr, "-s %d must be INT >= 0 (0 to use every block).\n", lodestarConfig -> sampleSize); 
        return -1;
    }
    // Check the projection dimension. The upper bound is tested in main(), once the
    //  number of samples in the VCF is known.
    if (lodestarConfig -> k < 1) {
        fprintf(stderr, "-k %d must be INT >= 1.\n", lodestarConfig -> k);
        return -1;
    }
    // Check maf.
    if (lodestarConfig -> maf < 0 || lodestarConfig -> maf > 1) { 
        fprintf(stderr, "--maf %lf must be in [0, 1].\n", lodestarConfig -> maf); 
        return -1;
    }
    // Check afMissing.
    if (lodestarConfig -> afMissing < 0 || lodestarConfig -> afMissing > 1) { 
        fprintf(stderr, "--afMissing %lf must be in [0, 1].\n", lodestarConfig -> afMissing); 
        return -1;
    }
    if (access(lodestarConfig -> inputFileName, F_OK) != 0) {
        fprintf(stderr, "Input file %s does not exist.\n", lodestarConfig -> inputFileName);
        return -1;
    }
    // Check target file if exists.
    if (lodestarConfig -> targetFileName != NULL && access(lodestarConfig -> targetFileName, F_OK) != 0) {
        fprintf(stderr, "-y %s does not exist.\n", lodestarConfig -> targetFileName); 
        return -1;
    }
    // Check haplotype and block size are valid.
    if (lodestarConfig -> BLOCK_SIZE <= 0) { 
        fprintf(stderr, "-b %d must be INT > 0.\n", lodestarConfig -> BLOCK_SIZE); 
        return -1;
    }
    // HAP_SIZE is a number of loci and BLOCK_SIZE is a number of base pairs, so the old
    //  `HAP_SIZE > BLOCK_SIZE` comparison related two different units and meant nothing.
    if (lodestarConfig -> HAP_SIZE <= 0) { 
        fprintf(stderr, "-H %d must be INT > 0.\n", lodestarConfig -> HAP_SIZE); 
        return -1;
    }
    // Check number of threads. 0 means "use every core we can see".
    if (lodestarConfig -> threads < 0) { 
        fprintf(stderr, "--threads %d must be INT >= 0 (0 for all available cores).\n", lodestarConfig -> threads); 
        return -1;
    }
    // --geo only means anything alongside -y.
    if (lodestarConfig -> geo && lodestarConfig -> targetFileName == NULL) {
        fprintf(stderr, "--geo has no effect without -y. Exiting!\n");
        return -1;
    }
    // Check drop threshhold.
    if (lodestarConfig -> dropThreshold < 1) {
        fprintf(stderr, "--drop %d must be INT > 0.\n", lodestarConfig -> dropThreshold); 
        return -1;
    }
    // All parameters are valid, return success.
    return 0;
}

// Print help menu for LODESTAR.
// Accepts: void.
// Returns: void.
void print_help() {
    fprintf(stderr, "\n");
    fprintf(stderr, "LODESTAR v%s\n", LODESTAR_VERSION);
    fprintf(stderr, "----------------------\n\n");
    fprintf(stderr, "Written by T. Quinn Smith\n");
    fprintf(stderr, "Principal Investigator: Zachary A. Szpiech\n");
    fprintf(stderr, "The Pennsylvania State University\n\n");
    fprintf(stderr, "Usage: lodestar [options] <input>.vcf.gz\n");
    fprintf(stderr, "Writes <basename>.tsv (per-block summary) and <basename>.json (coordinates).\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -h, --help              Print this menu and exit.\n");
    fprintf(stderr, "   -v, --version           Print the version and exit.\n");
    fprintf(stderr, "   -o, --out STR           Output file basename.\n");
    fprintf(stderr, "                               Default <input> with the .vcf(.gz) suffix removed.\n");
    fprintf(stderr, "   -H, --hap-size INT      Number of loci in a haplotype.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   -b, --block-size INT    The block size in base pairs.\n");
    fprintf(stderr, "                               Default 2 Mb.\n");
    fprintf(stderr, "   -k, --dimension INT     Dimension to project samples into.\n");
    fprintf(stderr, "                               Default 2. Must be less than number of samples in VCF.\n");
    fprintf(stderr, "   -y, --target file.tsv   A n-by-k tsv file containing user defined coordinates.\n");
    fprintf(stderr, "                               The set of points to use in Procrustes analysis.\n");
    fprintf(stderr, "                               Rows must be in the same order as the VCF samples.\n");
    fprintf(stderr, "   --geo                   Requires -y. Assume coordinates are (long, lat) in degrees and\n");
    fprintf(stderr, "                               convert to rectangular coordinates, where R = 6371000, using\n");
    fprintf(stderr, "                               (R * pi * sqrt(2) * long / 360, R * sqrt(2) * sin(lat))\n");
    fprintf(stderr, "   -t, --threads INT       Number of threads to use in computation.\n");
    fprintf(stderr, "                               Default 1. 0 uses all available cores.\n");
    fprintf(stderr, "   -d, --drop INT          Block is dropped if less than INT number of haplotypes are within block.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   --maf DOUBLE            Drops biallelic VCF records with a MAF less than threshold.\n");
    fprintf(stderr, "                               Default 0.05.\n");
    fprintf(stderr, "   --afMissing DOUBLE      Drops VCF records with fraction of genotypes missing greater than threshold.\n");
    fprintf(stderr, "                               Default 0.1\n");
    fprintf(stderr, "   -r, --replicates INT    Number of replicates for bootstrap. 0 for no bootstrap.\n");
    fprintf(stderr, "                               Default 1000.\n");
    fprintf(stderr, "   -s, --boot-blocks INT   Number of blocks resampled per bootstrap replicate.\n");
    fprintf(stderr, "                               Default 0, meaning the number of blocks kept in the input.\n");
    fprintf(stderr, "   --seed INT              Base RNG seed for the bootstrap. Thread i uses seed + i, so a\n");
    fprintf(stderr, "                               given (--seed, -t) pair reproduces exactly.\n");
    fprintf(stderr, "                               Default 42. The effective seed is echoed to the output.\n");
    fprintf(stderr, "   -q, --quiet             Suppress the per-block and per-replicate progress lines.\n");
    fprintf(stderr, "\n");
}

// Long options used for LODESTAR.
//  NOTE on renames from v1.0:
//   -h is now --help (it was the haplotype size, which is now -H/--hap-size), and
//   --samples is now --boot-blocks, because in bcftools/plink --samples selects samples
//   and it selected the bootstrap resample size here. The old long spellings
//   (--haplotypeSize, --samples) are kept as hidden aliases.
static ko_longopt_t long_options[] = {
    {"help",            ko_no_argument,         'h'},
    {"version",         ko_no_argument,         'v'},
    {"out",             ko_required_argument,   'o'},
    {"threads",         ko_required_argument,   't'},
    {"maf",             ko_required_argument,   309},
    {"afMissing",       ko_required_argument,   310},
    {"geo",             ko_no_argument,         311},
    {"seed",            ko_required_argument,   312},
    {"quiet",           ko_no_argument,         'q'},
    {"target",          ko_required_argument,   'y'},
    {"dimension",       ko_required_argument,   'k'},
    {"hap-size",        ko_required_argument,   'H'},
    {"haplotypeSize",   ko_required_argument,   'H'},
    {"block-size",      ko_required_argument,   'b'},
    {"blockSize",       ko_required_argument,   'b'},
    {"drop",            ko_required_argument,   'd'},
    {"replicates",      ko_required_argument,   'r'},
    {"boot-blocks",     ko_required_argument,   's'},
    {"samples",         ko_required_argument,   's'},
    {0, 0, 0}
};

LodestarConfig_t* init_lodestar_config(int argc, char *argv[]) {

    // Print help menu when no options/arguments were given.
    if (argc == 1) {
        print_help();
        exit(0);
    } 

    // Parse user options.
    const char *opt_str = "hvqH:b:k:t:y:d:o:s:r:";
    ketopt_t options = KETOPT_INIT;
    int c;

    // Pass through options. Check for options requiring an argument that were not given one.
    //  Also, check if any options are supplied that are not defined. If user supplies the help
    //  argument, print help menu and exit, likewise for version.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            // help and version are a successful request for information, so they exit 0
            //  here rather than returning NULL -- which main() now treats as an error.
            case 'h': print_help(); exit(0);
            case 'v': fprintf(stderr, "LODESTAR v%s\n", LODESTAR_VERSION); exit(0);
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return NULL;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return NULL;
        }
	}

    // Set defaults for LODESTAR configuration.
    LodestarConfig_t* lodestarConfig = calloc(1, sizeof(LodestarConfig_t));
    lodestarConfig -> inputFileName = NULL;
    lodestarConfig -> outputBasename = NULL;
    lodestarConfig -> HAP_SIZE = 1;
    lodestarConfig -> BLOCK_SIZE = 2000000;
    lodestarConfig -> k = 2;
    lodestarConfig -> threads = 1;
    lodestarConfig -> targetFileName = NULL;
    lodestarConfig -> dropThreshold = 1;
    lodestarConfig -> maf = 0.05;
    lodestarConfig -> afMissing = 0.1;
    lodestarConfig -> numReps = 1000;
    lodestarConfig -> sampleSize = 0;
    lodestarConfig -> geo = false;
    lodestarConfig -> seed = 42;
    lodestarConfig -> quiet = false;

    // Parse command line arguments.
    options = KETOPT_INIT;
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'o': lodestarConfig -> outputBasename = strdup(options.arg); break;
            case 'H': lodestarConfig -> HAP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'b': lodestarConfig -> BLOCK_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'k': lodestarConfig -> k = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 't': lodestarConfig -> threads = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'y': lodestarConfig -> targetFileName = strdup(options.arg); break;
            case 'd': lodestarConfig -> dropThreshold = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'r': lodestarConfig -> numReps = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 's': lodestarConfig -> sampleSize = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'q': lodestarConfig -> quiet = true; break;
            case 309: lodestarConfig -> maf = strtod(options.arg, (char**) NULL); break;
            case 310: lodestarConfig -> afMissing = strtod(options.arg, (char**) NULL); break;
            case 311: lodestarConfig -> geo = true; break;
            case 312: lodestarConfig -> seed = strtoul(options.arg, (char**) NULL, 10); break;
        }
	}

    // --threads 0 means "use every core we can see".
    if (lodestarConfig -> threads == 0) {
        long cores = sysconf(_SC_NPROCESSORS_ONLN);
        lodestarConfig -> threads = cores > 0 ? (int) cores : 1;
    }

    // If the input file was not give, exit.
    if (argc - options.ind != 1) {
        fprintf(stderr, "No input file given. Exiting!\n");
        destroy_lodestar_config(lodestarConfig);
        return NULL;
    }

    // Set the input filename.
    lodestarConfig -> inputFileName = strdup(argv[options.ind]);

    // Check to that user input is all valid.
    if (check_configuration(lodestarConfig) != 0) {
        destroy_lodestar_config(lodestarConfig);
        return NULL;
    }

    // If output base name was not specified, then use the input file basename.
    //  The old code used strstr(), which finds the FIRST ".vcf" -- so a path like
    //  /data/run.vcf.d/chr1.vcf.gz produced the basename "/data/run" -- and did not
    //  handle strstr() returning NULL, in which case it passed a large negative value
    //  as the size_t argument of strndup().
    if (lodestarConfig -> outputBasename == NULL) {
        const char* name = lodestarConfig -> inputFileName;
        // Only look for the suffix in the final path component.
        const char* slash = strrchr(name, '/');
        const char* base = slash != NULL ? slash + 1 : name;
        const char* suffix = strstr(base, ".vcf");
        if (suffix != NULL)
            lodestarConfig -> outputBasename = strndup(name, (size_t) (suffix - name));
        else
            lodestarConfig -> outputBasename = strdup(name);
    }

    // Copy command to echo in output files.
    kstring_t* cmd = calloc(1, sizeof(kstring_t));
    for (int i = 0; i < argc; i++)
        ksprintf(cmd, "%s ", argv[i]);
    lodestarConfig -> cmd = cmd -> s;
    free(cmd);

    return lodestarConfig;

}

void destroy_lodestar_config(LodestarConfig_t* lodestarConfig) {
    if (lodestarConfig == NULL)
        return;
    if (lodestarConfig -> targetFileName != NULL)
        free(lodestarConfig -> targetFileName);
    if (lodestarConfig -> inputFileName != NULL)
        free(lodestarConfig -> inputFileName);
    if (lodestarConfig -> outputBasename != NULL)
        free(lodestarConfig -> outputBasename);
    if (lodestarConfig -> cmd != NULL)
        free(lodestarConfig -> cmd);
    free(lodestarConfig);
}