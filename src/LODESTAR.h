// File: LODESTAR.h
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main LODESTAR operations.

#ifndef _LODESTAR_H_
#define _LODESTAR_H_

#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "BlockList.h"
#include "MatrixOperations.h"
#include <math.h>

// The number of elements in the upper triangle of a symmetric matrix.
#define PACKED_SIZE(N) ((N * (N + 1)) / 2)

// We store the upper triangle of a symmetric matrix in a one-dimensional array,
//  This is known as packed storage. Element Aij -> a_{i + j * (j + 1) / 2}
#define PACKED_INDEX(i, j) (i + j * (j + 1) / 2)

// IBS_t stores the three counts as an array (see BlockList.h), so the count for
//  numShared alleles is indexed directly. This replaces a three-way switch on
//  unpredictable data in the hottest loop of the program, without the aliasing
//  problems of casting the struct to unsigned int* and doing pointer arithmetic.
static inline void increment_ibs_value(IBS_t* ibs, int numShared) {
    ibs -> c[numShared]++;
}

// Calculate the number of shared alleles between two genotypes.
// Accepts:
//  Genotype_t s1 -> The genotype of the first sample.
//  Genotype_t s2 -> The genotype of the second sample.
// Returns: int, The number of shared alleles.
// The genotypes are unordered, so the number of shared alleles is the better of the
//  two pairings. Written branchlessly: the comparisons become flag registers rather
//  than jumps, which matters because this is called O(numHaplotypes * N^2) times.
//  Equivalent to the previous nested-if version for every input.
static inline int num_shared_alleles(Genotype_t s1, Genotype_t s2) {
    int direct = (s1.left == s2.left) + (s1.right == s2.right);
    int crossed = (s1.left == s2.right) + (s1.right == s2.left);
    return direct > crossed ? direct : crossed;
}


// Converts IBS counts to asd.
// Accepts:
//  IBS_t ibs -> The IBS counts.
// Returns: double, the transformed distance.
// L is the number of loci at which BOTH samples of this pair were called, which
//  differs between pairs whenever missingness does. Dividing by L puts every pair on
//  the same per-locus scale; without it a pair compared at fewer loci gets a
//  systematically smaller distance and the sample is pulled toward the centroid.
static inline double ibs_to_asd(IBS_t ibs) {
    double L = (double) (ibs.c[0] + ibs.c[1] + ibs.c[2]);
    if (L <= 0)
        return 0;
    return sqrt((4 * L - 2 * (ibs.c[1] + 2 * (double) ibs.c[2])) / L);
}

// Adds counts from right to left.
// Accepts:
//  IBS_t* left -> The accumulating counts.
//  IBS_t* right -> The increment counts.
// Returns: void.
static inline void add_ibs(IBS_t* left, IBS_t* right) {
    left -> c[0] += right -> c[0];
    left -> c[1] += right -> c[1];
    left -> c[2] += right -> c[2];
}

// Subtracts counts from right to left.
// Accepts:
//  IBS_t* left -> The accumulating counts.
//  IBS_t* right -> The decrementing counts.
// Returns: void.
static inline void sub_ibs(IBS_t* left, IBS_t* right) {
    left -> c[0] -= right -> c[0];
    left -> c[1] -= right -> c[1];
    left -> c[2] -= right -> c[2];
}


// Calculate IBS in blocks along the genome.
// Accepts:
//  VCFLocusParser_t* vcfFile -> The VCF file we are reading in.
//  HaplotypeEncoder_t* encoder -> The encoder for haplotypes.
//  int numSamples -> The number of samples.
//  int blockSize -> The size of the genome block in base-pairs.
//  int haplotypeSize -> The size of the haplotypes in number of loci.
//  int dropThreshold -> If numHaps within block is less than this value, we drop it.
//  int NUM_THREADS -> The number of threads to use in the computation.
//  bool quiet -> If set, suppress the per-block progress lines on stderr.
// Returns: BlockList_t*, the list of resulting blocks, or NULL on allocation failure.
BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int dropThreshold, int NUM_THREADS, bool quiet);

// Convert IBS blocks to ASD and perform Procrustes analysis
// Accepts:
//  BlockList_t* globalList -> Blocks along the genome with precomputed pairwise IBS.
//  double** y -> Centered/normalized target matrix to perform Procrustes against. If NULL, use genome-wide matrix and do not caluclate jackknife.
//  double* y0 -> Centroid of user defined points.
//  int k -> The dimension to reduce to.
//  int NUM_THREADS -> The number of threads to use in the computation.
//  int numReps -> The number of replicates for the bootstrap.
//  int sampleSize -> The number of blocks to sample for the bootstrap.
//  unsigned long seed -> Base RNG seed; thread i uses seed + i.
//  bool quiet -> If set, suppress the per-block and per-replicate progress lines.
// Returns: void.
void procrustes(BlockList_t* globalList, double** y, double* y0, int k, int NUM_THREADS, int numReps, int sampleSize, unsigned long seed, bool quiet);

#endif