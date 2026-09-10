// File: VCFLocusParser.h
// Date: 6 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parse a VCF file for samples' genotypes at each record.

#ifndef _VCF_LOCUS_PARSER_H_
#define _VCF_LOCUS_PARSER_H_

#include <zlib.h>
#include "kseq.h"
#include <stdbool.h>

// We use kseq as a stream to read in GZ files.
#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

typedef unsigned int Locus;
#define LEFT_ALLELE(a) (a >> 16)
#define RIGHT_ALLELE(a) (a & 0x0000FFFF)
#define MAX_NUM_ALLELES 1024

// Our structure that represents a VCFLocusParser.
typedef struct {
    // The name of the VCF file.
    char* fileName;
    // The GZ file. Note: Seamlessly works with uncompressed files too.
    gzFile file;
    // The stream to read in the GZ file.
    kstream_t* stream;
    // A buffer to hold a line of the VCF parser from the stream.
    kstring_t* buffer;

    // Number of samples in the VCF file.
    int numSamples;
    // An array to hold the names of each sample.
    char** sampleNames;

    // Flag to drop monomorphic sites.
    bool dropMonomorphicSites;
    // Minor-allele-frequency threshold.
    double maf;
    // Missing genotype threshold.
    double afMissing;
    // An array to hold the number of each allele at a locus.
    //  Used in maf calculation. Has MAX_NUM_ALLELES + 1 slots, since a missing
    //  genotype is encoded as allele index numAlleles.
    int* alleleCounts;

    // One-shot flags so a malformed input warns once rather than once per record.
    bool warnedFormat;
    bool warnedExtraColumns;

    // For convenience, we implement a priming read/peak operation.
    //  This allows us to easily test if the next record belongs to a different chromosome.
    char* nextChrom;
    unsigned int nextCoord;
    int nextNumAlleles;
    // Array that holds the genotypes for each of the samples.
    Locus* nextLocus;
} VCFLocusParser_t;

// Creates a VCFLocusParser_t structure.
// Accepts:
//  char* fileName -> The name of the VCF file.
//  double maf -> Records with a minor-allele-frequency < maf are dropped.
//  double afMissing -> Record with a proportion of missing genotypes >= afMissing are dropped.
//  bool dropMonomorphicSites -> If set, records with monomorphic genotypes are dropped.
// Returns: VCFLocusParser_t*, If fileName does not correspond to a file or regions could not be parsed, return NULL.
//              Otherwise, return the created structure.
VCFLocusParser_t* init_vcf_locus_parser(char* fileName, double maf, double afMissing, bool dropMonomorphicSites);

// Get the next locus from a parser.
// Accepts:
//  VCFLocusParser_t* parser -> The parser structure to read the VCF file.
//  char** chrom -> Sets the chromosome of the read-in record.
//  unsigned int* coord -> Sets the position of the read-in record.
//  int* numOfAlleles -> Sets the number of alleles at the read-in record.
//  Locus** genos -> Sets the array of the samples' genotypes at the read-in record.
// Returns: bool, True is another record exists and False otherwise.
bool get_next_locus(VCFLocusParser_t* parser, char** chrom, int* coord, int* numOfAlleles, Locus** genos);

// Checks if stream has reached end of the VCF file.
bool isEOF(VCFLocusParser_t* parser);

// Free all the memory allocated to a VCFLocusParser_t.
// Accepts:
//  VCFLocusParser_t* parser -> The parser to free.
// Returns: void.
void destroy_vcf_locus_parser(VCFLocusParser_t* parser);

// Encode a string genotype from the VCF file into a Locus.
//  This function will be called a lot.
// Accepts:
//  char* start -> The genotype to parse.
//  int numAlleles -> The number alleles at the given record.
// Returns: Locus, The encoded genotype.
static inline Locus parse_locus(char* start, int numAlleles) {
    // By default, both genotypes at a locus are missing.
    Locus locus = (numAlleles << 16) | numAlleles;
    char* next = start + 1;
    // If the left allele is not missing, then parse integer and set left genotype.
    if (start[0] != '.')
        locus = (strtol(start, &next, 10) << 16) | numAlleles;
    // If there is a second, non-missing genotype, parse and set right genotype.
    if ((next[0] == '|' || next[0] == '/') && next[1] != '.')
        locus = (locus & 0xFFFF0000) | strtol(next + 1, (char**) NULL, 10);
    return locus;
}

#endif