
// File: Main.c
// Date: 
// Version 1: 
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Run LODESTAR analysis.

#include "Interface.h"
#include "LODESTAR.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_json_matrix(FILE* out, double** X, int N, int K) {
    fprintf(out, "[\n");
    for (int i = 0; i < N; i++) {
        fprintf(out, "[ ");
        for (int j = 0; j < K - 1; j++) {
            fprintf(out, "%lf, ", X[i][j]);
        }
        fprintf(out, "%lf ]", X[i][K - 1]);
        if (i != N - 1)
            fprintf(out, ",\n");
        else
            fprintf(out, "\n");
    }
    fprintf(out, "]");
}

// Write a string as a JSON string body, escaping the characters JSON requires.
//  The command line is interpolated into the output, so a filename containing a quote
//  or a backslash used to produce a file that no JSON parser would accept.
void print_json_escaped(FILE* out, const char* s) {
    for (; s != NULL && *s != '\0'; s++) {
        switch (*s) {
            case '"':  fprintf(out, "\\\""); break;
            case '\\': fprintf(out, "\\\\"); break;
            case '\n': fprintf(out, "\\n"); break;
            case '\r': fprintf(out, "\\r"); break;
            case '\t': fprintf(out, "\\t"); break;
            default:
                if ((unsigned char) *s < 0x20)
                    fprintf(out, "\\u%04x", (unsigned char) *s);
                else
                    fputc(*s, out);
        }
    }
}

void print_json(LodestarConfig_t* lodestarConfig, BlockList_t* globalList, char** sampleNames, double** y, double* y0) {
    kstring_t* outName = calloc(1, sizeof(kstring_t));
    ksprintf(outName, "%s.json", lodestarConfig -> outputBasename);
    FILE* out = fopen(outName -> s, "w");
    if (out == NULL) {
        fprintf(stderr, "Could not write %s.\n", outName -> s);
        free(outName -> s); free(outName);
        return;
    }
    free(outName -> s); free(outName);

    fprintf(out, "{\n");

    // Print out all the global information.
    fprintf(out, "\"Version\": \"%s\",\n", LODESTAR_VERSION);
    fprintf(out, "\"Command\": \"");
    print_json_escaped(out, lodestarConfig -> cmd);
    fprintf(out, "\",\n");
    fprintf(out, "\"Seed\": %lu,\n", lodestarConfig -> seed);
    fprintf(out, "\"BootstrapReplicates\": %d,\n", globalList -> numReps);
    fprintf(out, "\"NumberOfBlocks\": %d,\n", globalList -> numBlocks);
    fprintf(out, "\"NumberOfBlocksDropped\": %d,\n", globalList -> numBlocks - globalList -> numKept);
    // Sample identifiers, in the row order of every coordinate matrix below. Without
    //  these every X matrix was an unlabelled block in VCF column order.
    fprintf(out, "\"Samples\": [ ");
    for (int i = 0; i < globalList -> numSamples; i++) {
        fprintf(out, "\"");
        print_json_escaped(out, sampleNames[i]);
        fprintf(out, "\"%s", i + 1 < globalList -> numSamples ? ", " : " ");
    }
    fprintf(out, "],\n");
    fprintf(out, "\"GlobalNumberOfLoci\": %d,\n", globalList -> numLoci);
    fprintf(out, "\"GlobalNumberOfHaplotypes\": %d,\n", globalList -> numHaps);
    // Spelling corrected from "GlobalVarainceCaptured"; the misspelling is also emitted
    //  for one release so existing parsers keep working.
    fprintf(out, "\"GlobalVarianceCaptured\": %lf,\n", globalList -> varCapt);
    fprintf(out, "\"GlobalVarainceCaptured\": %lf,\n", globalList -> varCapt);
    if (lodestarConfig -> targetFileName == NULL) {
        fprintf(out, "\"GlobalProcrustesStatistic\": null,\n");
        fprintf(out, "\"GlobalProcrustesStatisticPvalue\": null,\n");
    } else {
        fprintf(out, "\"GlobalProcrustesStatistic\": %lf,\n", globalList -> procrustesT);
        fprintf(out, "\"GlobalProcrustesStatisticPvalue\": %lf,\n", globalList -> pvalue);
    }
    fprintf(out, "\"GlobalX\":");
    print_json_matrix(out, globalList -> X, globalList -> numSamples, lodestarConfig -> k);
    fprintf(out, ",\n");
    if (lodestarConfig -> targetFileName == NULL) {
        fprintf(out, "\"Y\": null,\n");
        fprintf(out, "\"y0\": null,\n");
    } else {
        fprintf(out, "\"Y\":");
        print_json_matrix(out, y, globalList -> numSamples, lodestarConfig -> k);
        fprintf(out, ",\n");
        fprintf(out, "\"y0\": [ ");
        for (int i = 0; i < lodestarConfig -> k - 1; i++)
            fprintf(out, "%lf, ", y0[i]);
        fprintf(out, "%lf ],\n", y0[lodestarConfig -> k - 1]);
    }
    
    // Print out all the blocks.
    //  The separator is driven by whether we have already emitted a block, not by
    //  `temp -> next != NULL`. With the old test, a dropped block at the end of the list
    //  left a trailing comma before the closing bracket and the file was not valid JSON
    //  (jsonlite::fromJSON, and so lodestarPlots.R, could not read it).
    fprintf(out, "\"Blocks\": [");
    bool wroteBlock = false;
    for (Block_t* temp = globalList -> head; temp != NULL; temp = temp -> next) {
        if (temp -> isDropped)
            continue;
        if (wroteBlock)
            fprintf(out, ",\n");
        wroteBlock = true;
        fprintf(out, "\t{\n");
        fprintf(out, "\t\"BlockNumber\": %d,\n", temp -> blockNum);
        fprintf(out, "\t\"BlockNumberOnChromosome\": %d,\n", temp -> blockNumOnChrom);
        fprintf(out, "\t\"Chromosome\": \"%s\",\n", temp -> chrom);
        fprintf(out, "\t\"StartCoordinate\": %d,\n", temp -> startCoordinate);
        fprintf(out, "\t\"EndCoordinate\": %d,\n", temp -> endCoordinate);
        fprintf(out, "\t\"NumberOfLoci\": %d,\n", temp -> numLoci);
        fprintf(out, "\t\"NumberOfHaplotypes\": %d,\n", temp -> numHaps);
        // The isDropped branch that used to live here was unreachable -- dropped blocks
        //  are skipped above. A block whose cMDS did not converge has X == NULL, and its
        //  statistics are reported as null rather than as the -1 sentinel.
        if (temp -> X == NULL) {
            fprintf(out, "\t\"VarianceCaptured\": null,\n");
            fprintf(out, "\t\"ProcrustesStatistic\": null,\n");
            fprintf(out, "\t\"ProcrustesPValue\": null,\n");
        } else {
            fprintf(out, "\t\"VarianceCaptured\": %lf,\n", temp -> varCapt);
            fprintf(out, "\t\"ProcrustesStatistic\": %lf,\n", temp -> procrustesT);
            if (globalList -> numReps > 0)
                fprintf(out, "\t\"ProcrustesPValue\": %lf,\n", temp -> pvalue);
            else
                fprintf(out, "\t\"ProcrustesPValue\": null,\n");
        }
        fprintf(out, "\t\"X\": ");
        if (temp -> X != NULL)
            print_json_matrix(out, temp -> X, globalList -> numSamples, lodestarConfig -> k);
        else 
            fprintf(out, "null");
        fprintf(out, "}");
    }
    fprintf(out, "\n\t]\n");
    fprintf(out, "}\n");

    fclose(out);
}

// Write one double, or NA when it was not computed.
//  The old output used -1 as the sentinel for "not computed", which reads as data to
//  R's read.delim and pandas' read_csv. NA is what both expect.
static void print_tsv_double(FILE* out, double value, bool computed, const char* terminator) {
    if (computed)
        fprintf(out, "%lf%s", value, terminator);
    else
        fprintf(out, "NA%s", terminator);
}

void print_summary(LodestarConfig_t* lodestarConfig, BlockList_t* globalList) {
    kstring_t* outName = calloc(1, sizeof(kstring_t));
    ksprintf(outName, "%s.tsv", lodestarConfig -> outputBasename);
    FILE* out = fopen(outName -> s, "w");
    if (out == NULL) {
        fprintf(stderr, "Could not write %s.\n", outName -> s);
        free(outName -> s); free(outName);
        return;
    }
    free(outName -> s); free(outName);

    // Echo the command and the run parameters needed to reproduce it.
    fprintf(out, "#LODESTAR v%s\n", LODESTAR_VERSION);
    fprintf(out, "#%s\n", lodestarConfig -> cmd);
    fprintf(out, "#seed=%lu threads=%d replicates=%d blocks=%d blocksDropped=%d\n",
            lodestarConfig -> seed, lodestarConfig -> threads, globalList -> numReps,
            globalList -> numBlocks, globalList -> numBlocks - globalList -> numKept);

    // Print header.
    fprintf(out, "BlockNum\tBlockNumOnChr\tChr\tStart\tEnd\tNumLoci\tNumHaps\tVarianceCaptured\tProcrustesStatistic\tP-Value\n");

    // Print out each block.
    for(Block_t* temp = globalList -> head; temp != NULL; temp = temp -> next) {
        if (temp -> isDropped)
            continue;
        fprintf(out, "%d\t%d\t%s\t%d\t%d\t%d\t%d\t", temp -> blockNum, temp -> blockNumOnChrom, temp -> chrom, temp -> startCoordinate, temp -> endCoordinate, temp -> numLoci, temp -> numHaps);
        bool converged = temp -> X != NULL;
        print_tsv_double(out, temp -> varCapt, converged, "\t");
        print_tsv_double(out, temp -> procrustesT, converged, "\t");
        print_tsv_double(out, temp -> pvalue, converged && globalList -> numReps > 0, "\n");
    }
    fprintf(out, "0\t0\tGLOBAL\t0\t0\t%d\t%d\t", globalList -> numLoci, globalList -> numHaps);
    print_tsv_double(out, globalList -> varCapt, true, "\t");
    bool haveGlobalT = globalList -> procrustesT != -1;
    print_tsv_double(out, globalList -> procrustesT, haveGlobalT, "\t");
    print_tsv_double(out, globalList -> pvalue, haveGlobalT && globalList -> numReps > 0, "\n");

    fclose(out);
}

double** open_target_file(char* targetFileName, int N, int K) {
    FILE* targetFile = fopen(targetFileName, "r");
    if (targetFile == NULL) {
        fprintf(stderr, "Could not open -y %s.\n", targetFileName);
        return NULL;
    }

    // Create target matrix.
    double** targetPoints = init_matrix(N, K);
    if (targetPoints == NULL) {
        fclose(targetFile);
        return NULL;
    }

    // Read in values from target.
    int numLines = 0;
    size_t length = 0;
    char* line = NULL;
    bool bad = false;
    while (getline(&line, &length, targetFile) != -1) {
        // Skip blank lines so a trailing newline does not look like a short row.
        if (line[0] == '\n' || line[0] == '\r' || line[0] == '\0')
            continue;
        // NOTE: the bound must be tested BEFORE the write. The old loop wrote
        //  targetPoints[numLines][i] and only then checked numLines == N, so a file with
        //  more rows than samples overflowed the row-pointer array
        //  (AddressSanitizer: heap-buffer-overflow).
        if (numLines >= N) {
            fprintf(stderr, "-y %s has more than %d rows.\n", targetFileName, N);
            bad = true;
            break;
        }
        int i = 0;
        char* tok = strtok(line, "\t\n\r");
        for (; tok != NULL && i < K; i++) {
            // atof() silently returns 0 for a non-numeric token; strtod lets us tell.
            char* end = NULL;
            targetPoints[numLines][i] = strtod(tok, &end);
            if (end == tok) {
                fprintf(stderr, "-y %s row %d column %d (\"%s\") is not a number.\n", targetFileName, numLines + 1, i + 1, tok);
                bad = true;
                break;
            }
            tok = strtok(NULL, "\t\n\r");
        }
        // Reject both short rows and rows with trailing extra columns.
        if (bad)
            break;
        if (i != K || tok != NULL) {
            fprintf(stderr, "-y %s row %d does not have exactly %d columns.\n", targetFileName, numLines + 1, K);
            bad = true;
            break;
        }
        numLines++;
    }
    fclose(targetFile);
    free(line);

    if (bad || numLines != N) {
        if (!bad)
            fprintf(stderr, "-y %s has %d rows; expected %d (one per VCF sample).\n", targetFileName, numLines, N);
        destroy_matrix(targetPoints, N);
        return NULL;
    }
    return targetPoints;
}

void convertToRectangular(double** y, int N) {
    for (int i = 0; i < N; i++) {
        y[i][0] = 6371000 * M_PI * sqrt(2) * y[i][0] / 360.0;
        y[i][1] = 6371000 * sqrt(2) * sin(y[i][1]);
    }
}

int main (int argc, char *argv[]) {

    // Get configuration.
    LodestarConfig_t* lodestarConfig = init_lodestar_config(argc, argv);
    // A NULL config now means the options were invalid. --help, --version and a bare
    //  invocation exit(0) inside init_lodestar_config, so reaching here with NULL is an
    //  error and must not report success to the shell.
    if (lodestarConfig == NULL)
        return -1;
   
    // We create the VCF parser and haplotype encoder.
    //  init_vcf_locus_parser() returns NULL for an unreadable or non-VCF input; the old
    //  code dereferenced it immediately.
    VCFLocusParser_t* parser = init_vcf_locus_parser(lodestarConfig -> inputFileName, lodestarConfig -> maf, lodestarConfig -> afMissing, true);
    if (parser == NULL) {
        fprintf(stderr, "Could not read %s as a gzipped VCF. Exiting!\n", lodestarConfig -> inputFileName);
        destroy_lodestar_config(lodestarConfig);
        return -1;
    }
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(parser -> numSamples);
    
    // Check to make sure k is valid.
    if (lodestarConfig -> k >= parser -> numSamples) {
        fprintf(stderr, "k = %d >= numSamples %d. Exiting!\n", lodestarConfig -> k, parser -> numSamples);
        destroy_lodestar_config(lodestarConfig);
        destroy_haplotype_encoder(encoder);
        destroy_vcf_locus_parser(parser);
        return -1;
    }

    // User defined points to perform Procrustes against.
    double** y = NULL;
    double* y0 = NULL;

    // If target file supplied, make sure it is valid, at least numSamples-by-k
    if (lodestarConfig -> targetFileName != NULL) {
        y = open_target_file(lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
        if (y == NULL) {
            fprintf(stderr, "-y %s is not an %d-by-%d matrix. Exiting!\n", lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
            destroy_lodestar_config(lodestarConfig);
            destroy_haplotype_encoder(encoder);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
        if (lodestarConfig -> geo && lodestarConfig -> k != 2) {
            fprintf(stderr, "--geo requires two dimensions for the matrix given with -y. Exiting!\n");
            destroy_matrix(y, parser -> numSamples);
            destroy_lodestar_config(lodestarConfig);
            destroy_haplotype_encoder(encoder);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
        if (lodestarConfig -> geo)
            convertToRectangular(y, parser -> numSamples);
        // Center and normalize.
        y0 = calloc(lodestarConfig -> k, sizeof(double));
        center_matrix(y, y0, parser -> numSamples, lodestarConfig -> k);
        normalize_matrix(y, parser -> numSamples, lodestarConfig -> k);
    }
    
    // Partition genome into blocks and calculate IBS within the blocks.
    if (!lodestarConfig -> quiet)
        fprintf(stderr, "Beginning blocking algorithm ...\n");
    BlockList_t* globalList = block_allele_sharing(parser, encoder, encoder -> numSamples, lodestarConfig -> BLOCK_SIZE, lodestarConfig -> HAP_SIZE, lodestarConfig -> dropThreshold, lodestarConfig -> threads, lodestarConfig -> quiet);
    if (globalList == NULL) {
        destroy_matrix(y, parser -> numSamples);
        if (y0 != NULL) free(y0);
        destroy_lodestar_config(lodestarConfig);
        destroy_haplotype_encoder(encoder);
        destroy_vcf_locus_parser(parser);
        return -1;
    }
    if (globalList -> numBlocks == 0) {
        fprintf(stderr, "No blocks were formed from %s. Is the VCF empty after filtering? Exiting!\n", lodestarConfig -> inputFileName);
        destroy_matrix(y, parser -> numSamples);
        if (y0 != NULL) free(y0);
        destroy_block_list(globalList);
        destroy_lodestar_config(lodestarConfig);
        destroy_haplotype_encoder(encoder);
        destroy_vcf_locus_parser(parser);
        return -1;
    }

    // Convert IBS to ASD and compute MDS on global.
    if (!lodestarConfig -> quiet)
        fprintf(stderr, "\nFinished blocking algorithm. Starting genome-wide MDS calculations ...\n");
    double* asd = calloc(PACKED_SIZE(encoder -> numSamples), sizeof(double));
    RealSymEigen_t* eigen = init_real_sym_eigen(encoder -> numSamples);
    for (int i = 0; i < encoder -> numSamples; i++)
        for (int j = i + 1; j < encoder -> numSamples; j++)
            asd[PACKED_INDEX(i, j)] = ibs_to_asd(globalList -> alleleCounts[PACKED_INDEX(i, j)]);
    globalList -> X = init_matrix(encoder -> numSamples, lodestarConfig -> k);
    globalList -> varCapt = compute_classical_mds(eigen, asd, lodestarConfig -> k, globalList -> X);
    // cMDS results are already centered. Normalize X for symmetric Procrsutes statistic.
    normalize_matrix(globalList -> X, parser -> numSamples, lodestarConfig -> k);
    if (y != NULL)
        globalList -> procrustesT = procrustes_statistic(globalList -> X, NULL, y, y0, eigen, eigen -> N, lodestarConfig -> k, true);
    destroy_real_sym_eigen(eigen);
    free(asd);

    // Convert IBS to ASD and calculate jackknifed procrustes statistic.
    if (!lodestarConfig -> quiet)
        fprintf(stderr, "Finished genome-wide MDS calulations. Starting Procrustes ...\n\n");
    // Default resample size is the number of blocks that were KEPT, since dropped blocks
    //  cannot be drawn.
    if (lodestarConfig -> sampleSize == 0)
        lodestarConfig -> sampleSize = globalList -> numKept;
    procrustes(globalList, y, y0, lodestarConfig -> k, lodestarConfig -> threads, lodestarConfig -> numReps, lodestarConfig -> sampleSize, lodestarConfig -> seed, lodestarConfig -> quiet);
    if (lodestarConfig -> numReps == 0)
        fprintf(stderr, "Writing results to output files...\n");
    else
        fprintf(stderr, "\nFinished Bootstrap. Writing results to output files...\n");

    // Print summary and JSON file.
    print_summary(lodestarConfig, globalList);
    print_json(lodestarConfig, globalList, parser -> sampleNames, y, y0);

    // Free used memory.
    if (y0 != NULL) free(y0);
    destroy_matrix(y, encoder -> numSamples);
    destroy_block_list(globalList);
    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);
    destroy_lodestar_config(lodestarConfig);
    fprintf(stderr, "Done!\n");
    return 0;
}