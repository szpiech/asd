// File: BlockList.c
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#include "BlockList.h"
#include <string.h>
#include <stdlib.h>

Block_t* init_block(char* chrom, int startCoordinate, int numSamples) {
    Block_t* block = calloc(1, sizeof(Block_t));
    block -> isDropped = false;
    block -> chrom = strdup(chrom);
    block -> startCoordinate = startCoordinate;
    block -> numHaps = 0;
    block -> numLoci = 0;
    block -> numSamples = numSamples;
    block -> procrustesT = -1;
    block -> pvalue = -1;
    block -> alleleCounts = NULL;
    block -> X = NULL;
    block -> next = NULL;
    return block;
}

BlockList_t* init_block_list(int numSamples) {
    BlockList_t* blockList = calloc(1, sizeof(BlockList_t));
    blockList -> numSamples = numSamples;
    blockList -> head = NULL;
    blockList -> tail = NULL;
    blockList -> numHaps = 0;
    blockList -> numLoci = 0;
    blockList -> procrustesT = -1;
    blockList -> pvalue = -1;
    blockList -> alleleCounts = NULL;
    blockList -> X = NULL;
    return blockList;
}

void append_block(BlockList_t* blockList, Block_t* block) {
    if (blockList -> numBlocks == 0) {
        blockList -> head = block;
        blockList -> tail = block;
    } else {
        blockList -> tail -> next = block;
        blockList -> tail = block;
    }
    blockList -> numBlocks++;
}

// Order two blocks by their global block number.
static int compare_block_num(const void* a, const void* b) {
    int x = (*(Block_t* const*) a) -> blockNum;
    int y = (*(Block_t* const*) b) -> blockNum;
    return (x > y) - (x < y);
}

int finalize_block_list(BlockList_t* blockList) {
    if (blockList -> numBlocks == 0)
        return 0;

    // Collect the list into an array. Threads append in completion order, so this is
    //  where the ordering is restored -- with qsort rather than the O(numBlocks^2)
    //  selection sort the threaded path used to run on the linked list.
    blockList -> blocks = calloc(blockList -> numBlocks, sizeof(Block_t*));
    if (blockList -> blocks == NULL)
        return -1;
    int n = 0;
    for (Block_t* temp = blockList -> head; temp != NULL; temp = temp -> next)
        blockList -> blocks[n++] = temp;
    qsort(blockList -> blocks, blockList -> numBlocks, sizeof(Block_t*), compare_block_num);

    // Relink in sorted order so anything that still walks head->next sees block order.
    for (int i = 0; i + 1 < blockList -> numBlocks; i++)
        blockList -> blocks[i] -> next = blockList -> blocks[i + 1];
    blockList -> blocks[blockList -> numBlocks - 1] -> next = NULL;
    blockList -> head = blockList -> blocks[0];
    blockList -> tail = blockList -> blocks[blockList -> numBlocks - 1];

    // Index of the blocks that survived the drop threshold.
    blockList -> kept = calloc(blockList -> numBlocks, sizeof(Block_t*));
    if (blockList -> kept == NULL)
        return -1;
    blockList -> numKept = 0;
    for (int i = 0; i < blockList -> numBlocks; i++)
        if (!blockList -> blocks[i] -> isDropped)
            blockList -> kept[blockList -> numKept++] = blockList -> blocks[i];

    return 0;
}

void destroy_block(Block_t* block) {
    if (block -> X != NULL) {
        for (int i = 0; i < block -> numSamples; i++) {
            free(block -> X[i]);
        }
        free(block -> X);
    }
    if (block -> alleleCounts != NULL)
        free(block -> alleleCounts);
    free(block -> chrom);
    free(block);
}

void destroy_block_list(BlockList_t* blockList) {
    if (blockList == NULL)
        return;
    if (blockList -> blocks != NULL)
        free(blockList -> blocks);
    if (blockList -> kept != NULL)
        free(blockList -> kept);
    if (blockList -> samplingDistribution != NULL)
        free(blockList -> samplingDistribution);
    if (blockList -> X != NULL) {
        for (int i = 0; i < blockList -> numSamples; i++) {
            free(blockList -> X[i]);
        }
        free(blockList -> X);
    }
    if (blockList -> alleleCounts != NULL)
        free(blockList -> alleleCounts);
    Block_t* temp = NULL;
    for (int i = 0; i < blockList -> numBlocks; i++) {
        temp = blockList -> head;
        blockList -> head = blockList -> head -> next;
        destroy_block(temp);
    }
    free(blockList);
}