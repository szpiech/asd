// File: BlockList.h
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#ifndef _BLOCK_LIST_H_
#define _BLOCK_LIST_H_

#include <stdbool.h>

// Our structure that tracks the counts of loci
//  with 0, 1, and 2 alleles IBS between samples.
//  c[n] is the number of loci at which the pair shared n alleles, so the count can be
//  indexed directly by the result of num_shared_alleles() (see LODESTAR.h).
typedef struct {
    unsigned int c[3];
} IBS_t;

// A node in the BlockList.
typedef struct Block {
    // Allele counts and lower dimensional representation for the block.
    IBS_t* alleleCounts;
    double** X;

    // Procrustes t-statistic, t-statistic if bootstrap is executed, and effective rank of points.
    double procrustesT;
    double sampleT;
    double varCapt;
    double pvalue;

    // Block attributes.
    int blockNum;
    int blockNumOnChrom;
    char* chrom;
    int startCoordinate;
    int endCoordinate;
    int numHaps;
    int numLoci;
    bool isDropped;
    int numSamples;
    // The next block in the list.
    struct Block* next;
} Block_t;

// A list of blocks.
typedef struct BlockList {
    // Global counts.
    IBS_t* alleleCounts;
    double** X;
    double varCapt;

    // If bootstrap is computed.
    double procrustesT;
    double pvalue;
    double* samplingDistribution;
    // The number of bootstrap replicates that actually completed (an MDS that failed to
    //  converge is not counted). This is what the p-value denominators use, and it is
    //  echoed into the output so a reader can tell a short bootstrap from a full one.
    int numReps;

    // Global attributes.
    int numSamples;
    int numBlocks;
    int numHaps;
    int numLoci;
    Block_t* head;
    Block_t* tail;

    // Indexable views over the list, built once by finalize_block_list().
    //  blocks holds every block in blockNum order; kept holds only the blocks that
    //  were not dropped. The bootstrap samples uniformly from kept in O(1), which
    //  replaces an O(numBlocks) linked-list walk per draw (and a rejection loop that
    //  could run off the end of the list).
    Block_t** blocks;
    Block_t** kept;
    int numKept;
} BlockList_t;

// Creates a list of blocks.
// Accepts:
//  int numSamples -> The number of samples.
// Returns: An empty block list.
BlockList_t* init_block_list(int numSamples);

// Creates a block.
// Acccepts:
//  char* chrom -> The chromosome the block is on.
//  int startCoordinate -> The coordinate of the first record in the block.
//  int numSamples -> The number of samples.
// Returns: Block_t*, the new block.
Block_t* init_block(char* chrom, int startCoordinate, int numSamples);

// Adds a block to the block list.
// Accepts:
//  BlockList_t* blockList -> The list to add the block to.
//  Block_t* block -> Adds the block to the end of the list.
// Returns: void.
void append_block(BlockList_t* blockList, Block_t* block);

// Sorts the list into blockNum order and builds the blocks/kept index arrays.
//  Must be called after all blocks have been appended and before the Procrustes
//  and bootstrap passes.
// Accepts:
//  BlockList_t* blockList -> The list to finalize.
// Returns: int, 0 on success, -1 if allocation failed.
int finalize_block_list(BlockList_t* blockList);

// Frees the memory occupied by the block list.
// Accepts:
//  BlockList_t* blockList -> The block list to destroy.
// Returns: void.
void destroy_block_list(BlockList_t* blockList);

#endif