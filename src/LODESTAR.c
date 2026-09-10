// File: LODESTAR.c
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main LODESTAR operations.

#include "LODESTAR.h"
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

// Used to mutually exclude threads' access to the VCF file
//  and overlapping genotypes for the current window.
pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;

// Used to mutually exclude access to the global list.
pthread_mutex_t listLock = PTHREAD_MUTEX_INITIALIZER;

// Used to mutually exclude access to the global counts.
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

// Saves a few lines to swap values.
#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

// Set by block_allele_sharing()/procrustes() from the --quiet option. When set, the
//  per-block and per-replicate progress lines are suppressed; warnings and errors
//  still go to stderr. A genome-scale run otherwise writes one unbuffered line per
//  block and per bootstrap replicate into the job log.
static bool beQuiet = false;

#define PROGRESS(...) do { if (!beQuiet) fprintf(stderr, __VA_ARGS__); } while (0)

// A list of haplotypes within a block.
typedef struct BlockGenotypes {
    Genotype_t* genotypes;
    struct BlockGenotypes* next;
} BlockGenotypes_t;

// What is given to each thread to block the genome.
typedef struct BlockCounts {

    // A simple linked list to add haplotypes to the block.
    //  We save a little time by not reallocating memory for each block.
    BlockGenotypes_t* head;
    BlockGenotypes_t* tail;
    int numHapsInGenotypeList;
    int numHaps;

    // Protected by mutex for reading in.
    VCFLocusParser_t* vcfFile;
    HaplotypeEncoder_t* encoder;
    int numSamples;
    int blockSize;
    int haplotypeSize;
    int dropThreshold;

    // Protected by mutex for adding values to global.
    BlockList_t* globalList;
    int* blockNum;

    // Thread-private scratch: the indices of the samples with a fully called genotype
    //  at the current haplotype. Sized numSamples, reused for every haplotype.
    int* callable;

} BlockCounts_t;

// A simple method to add a new haplotype to the linked list within a block.
void add_haplotype_to_block(BlockCounts_t* blockCounts) {
    BlockGenotypes_t* temp = calloc(1, sizeof(BlockGenotypes_t));
    temp -> genotypes = calloc(blockCounts -> numSamples, sizeof(Genotype_t));
    temp -> next = NULL;
    if (blockCounts -> numHapsInGenotypeList == 0) {
        blockCounts -> head = temp;
        blockCounts -> tail = temp;
    } else {
        blockCounts -> tail -> next = temp;
        blockCounts -> tail = temp;
    }
    blockCounts -> numHapsInGenotypeList++;
}

// Frees the linked list of haplotypes within a block.
void destroy_block_counts(BlockCounts_t* blockCounts) {
    if (blockCounts == NULL)
        return;
    if (blockCounts -> callable != NULL)
        free(blockCounts -> callable);
    if (blockCounts -> head != NULL) {
        BlockGenotypes_t* temp = NULL;
        for (int i = 0; i < blockCounts -> numHapsInGenotypeList; i++) {
            temp = blockCounts -> head;
            blockCounts -> head = blockCounts -> head -> next;
            free(temp -> genotypes);
            free(temp);
        }
    }
    if (blockCounts != NULL)
        free(blockCounts);
}

// For each thread, we read in a block of haplotypes and compute ASD.
void* partition(void* arg) {
    BlockCounts_t* blockCounts = (BlockCounts_t*) arg;

    Genotype_t* genoTemp;
    BlockGenotypes_t* blockTemp;

    // Points to the new block we add to the global list.
    Block_t* block;

    // Holds global counts.
    IBS_t* globalCounts = calloc(PACKED_SIZE(blockCounts -> numSamples), sizeof(IBS_t));

    while (true) {

        pthread_mutex_lock(&genomeLock);
        if (isEOF(blockCounts -> vcfFile)) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        
        // Reset everything for next block.
        block = init_block(blockCounts -> vcfFile -> nextChrom, blockCounts -> vcfFile -> nextCoord, blockCounts -> numSamples);
        blockCounts -> numHaps = 0;
        blockTemp = blockCounts -> head;
        bool isOnSameChrom = true;
        int endOfBlock = ((int) ((blockCounts -> vcfFile -> nextCoord - 1) / (double) blockCounts -> blockSize) + 1) * blockCounts -> blockSize;
        
        // Read in the next block.
        while (isOnSameChrom && blockCounts -> vcfFile -> nextCoord <= endOfBlock) {
            isOnSameChrom = get_next_haplotype(blockCounts -> vcfFile, blockCounts -> encoder, blockCounts -> haplotypeSize);
            // If we need more haplotypes, then we allocate the necessary memory.
            if (blockCounts -> numHaps == blockCounts -> numHapsInGenotypeList) {
                add_haplotype_to_block(blockCounts);
                blockTemp = blockCounts -> tail;
            }
            // Insert the newly read in haplotype.
            SWAP(blockTemp -> genotypes, blockCounts -> encoder -> genotypes, genoTemp);
            blockTemp = blockTemp -> next;
            blockCounts -> numHaps++;
            block -> numHaps++;
            block -> numLoci += blockCounts -> encoder -> numLoci;
        }
        // Finish the block.
        block -> endCoordinate = blockCounts -> encoder -> endCoord;
        block -> blockNum = *(blockCounts -> blockNum);
        *(blockCounts -> blockNum) = *(blockCounts -> blockNum) + 1;
        if (block -> numHaps < blockCounts -> dropThreshold)
            block -> isDropped = true;
        pthread_mutex_unlock(&genomeLock);

        // Calculate IBS for all haplotypes within the block and accumulate global if block is not dropped.
        if (!block -> isDropped) {
            block -> alleleCounts = calloc(PACKED_SIZE(blockCounts -> numSamples), sizeof(IBS_t));
            IBS_t* blockCountsArray = block -> alleleCounts;
            int N = blockCounts -> numSamples;
            int* callable = blockCounts -> callable;
            BlockGenotypes_t* locus = blockCounts -> head;
            for (int l = 0; l < blockCounts -> numHaps; l++) {
                Genotype_t* g = locus -> genotypes;
                // A pair contributes at this haplotype only when all four haplotypes are
                //  called. The old test was `g[i].left != MISSING || g[j].right != MISSING`
                //  -- an OR, comparing i's LEFT against j's RIGHT -- so a pair involving a
                //  missing sample was still counted, num_shared_alleles(MISSING, x)
                //  returned 0, and every missing genotype was scored as maximal
                //  dissimilarity instead of being skipped.
                //  .left and .right are collapsed to MISSING independently by add_locus(),
                //  so both must be tested. Doing it once per sample here costs N tests
                //  rather than N^2 / 2 inside the pair loop.
                int numCallable = 0;
                for (int s = 0; s < N; s++)
                    if (g[s].left != MISSING && g[s].right != MISSING)
                        callable[numCallable++] = s;
                // PACKED_INDEX(i, j) = i + j * (j + 1) / 2, so holding j fixed and walking
                //  i gives (near-)unit-stride access into the packed array. The previous
                //  i-outer/j-inner order strode by j+1 and thrashed cache for any N whose
                //  packed array exceeds L2 (6 MB at N=1000).
                for (int cj = 1; cj < numCallable; cj++) {
                    int j = callable[cj];
                    IBS_t* column = blockCountsArray + (j * (j + 1) / 2);
                    Genotype_t gj = g[j];
                    for (int ci = 0; ci < cj; ci++)
                        increment_ibs_value(&(column[callable[ci]]),
                                            num_shared_alleles(g[callable[ci]], gj));
                }
                locus = locus -> next;
            }
            // The genome-wide counts are by construction the sum of the per-block counts,
            //  so accumulate once per block rather than once per (locus, pair). This is
            //  O(numBlocks * N^2) instead of O(numHaplotypes * N^2) and halves the stores
            //  in the hottest loop above.
            for (int p = 0; p < PACKED_SIZE(N); p++)
                add_ibs(&(globalCounts[p]), &(blockCountsArray[p]));
        }

        // Append the block to the global list.
        pthread_mutex_lock(&listLock);
        if (block -> isDropped)
            PROGRESS("Block on %s from %d to %d contains %d haplotypes. Dropped block.\n", block -> chrom, block -> startCoordinate, block -> endCoordinate, block -> numHaps);
        else
            PROGRESS("Finished IBS for block number %d on %s from %d to %d.\n", block -> blockNum, block -> chrom, block -> startCoordinate, block -> endCoordinate);
        append_block(blockCounts -> globalList, block);
        pthread_mutex_unlock(&listLock);
    }

    // Accumulate global counts.
    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < blockCounts -> numSamples; i++)
        for (int j = i + 1; j < blockCounts -> numSamples; j++)
            add_ibs(&(blockCounts -> globalList -> alleleCounts[PACKED_INDEX(i, j)]), &(globalCounts[PACKED_INDEX(i, j)]));
    pthread_mutex_unlock(&globalLock);

    destroy_block_counts(blockCounts);
    free(globalCounts);

    return NULL;
}

// Allocate and populate the per-thread argument block for partition().
//  All threads share the parser, encoder, and global list; everything else is private.
static BlockCounts_t* init_block_counts(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder,
                                        int numSamples, int blockSize, int haplotypeSize,
                                        int dropThreshold, BlockList_t* globalList, int* blockNum) {
    BlockCounts_t* blockCounts = calloc(1, sizeof(BlockCounts_t));
    blockCounts -> vcfFile = vcfFile;
    blockCounts -> encoder = encoder;
    blockCounts -> blockSize = blockSize;
    blockCounts -> numSamples = numSamples;
    blockCounts -> haplotypeSize = haplotypeSize;
    blockCounts -> dropThreshold = dropThreshold;
    blockCounts -> globalList = globalList;
    blockCounts -> blockNum = blockNum;
    blockCounts -> callable = calloc(numSamples, sizeof(int));
    return blockCounts;
}

BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int dropThreshold, int NUM_THREADS, bool quiet) {

    beQuiet = quiet;

    BlockList_t* globalList = init_block_list(numSamples);
    globalList -> alleleCounts = calloc(PACKED_SIZE(numSamples), sizeof(IBS_t));
    int* blockNum = calloc(1, sizeof(int));
    *blockNum = 1;

    if (NUM_THREADS == 1) {
        partition((void*) init_block_counts(vcfFile, encoder, numSamples, blockSize,
                                            haplotypeSize, dropThreshold, globalList, blockNum));
    } else {
        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, partition,
                           (void*) init_block_counts(vcfFile, encoder, numSamples, blockSize,
                                                     haplotypeSize, dropThreshold, globalList, blockNum));
        partition((void*) init_block_counts(vcfFile, encoder, numSamples, blockSize,
                                            haplotypeSize, dropThreshold, globalList, blockNum));
        // Wait for all threads to finish.
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);
    }
    free(blockNum);

    // Restore block order and build the indexable views. This replaces the
    //  O(numBlocks^2) selection sort that used to run over the linked list, and gives
    //  the bootstrap O(1) uniform sampling over the non-dropped blocks.
    if (finalize_block_list(globalList) != 0) {
        fprintf(stderr, "Could not allocate the block index. Exiting!\n");
        destroy_block_list(globalList);
        return NULL;
    }

    // Assign blockNumOnChrom and count the global number of haplotypes and loci.
    //  NOTE: the old loop condition was `temp -> next != NULL`, so the final block was
    //  never counted and GlobalNumberOfLoci / GlobalNumberOfHaplotypes were short by it.
    int blockNumOnChrom = 1;
    globalList -> numHaps = 0;
    globalList -> numLoci = 0;
    for (int i = 0; i < globalList -> numBlocks; i++) {
        Block_t* temp = globalList -> blocks[i];
        temp -> blockNumOnChrom = blockNumOnChrom;
        if (temp -> next != NULL && strcmp(temp -> chrom, temp -> next -> chrom) != 0)
            blockNumOnChrom = 1;
        else
            blockNumOnChrom++;
        if (!temp -> isDropped) {
            globalList -> numHaps += temp -> numHaps;
            globalList -> numLoci += temp -> numLoci;
        }
    }

    return globalList;
}

typedef struct BlockProcrustes {
    BlockList_t* globalList;
    double** y;
    double* y0;
    int k;
    // Points to the current block in globalList for the next thread to operate on.
    Block_t** current;
    // For the bootstrap.
    int numReps;
    int sampleSize;
    int* currentReplicate;
    // Deterministic per-thread RNG seed, derived in procrustes() from the user's --seed.
    unsigned long seed;
} BlockProcrustes_t;

// Randomly sample a non-dropped block, uniformly and in O(1).
//  The previous version walked the linked list by index, and -- because `temp` was
//  initialised outside the rejection loop -- resumed the walk from wherever the last
//  rejected draw stopped. That made the draw non-uniform and ran off the end of the
//  list (dereferencing NULL) whenever any block was dropped.
static Block_t* get_random_block(gsl_rng* r, BlockList_t* globalList) {
    int i = (int) (globalList -> numKept * gsl_rng_uniform(r));
    // gsl_rng_uniform() is in [0, 1), but guard the boundary anyway.
    if (i >= globalList -> numKept)
        i = globalList -> numKept - 1;
    return globalList -> kept[i];
}

void* procrustes_bootstrap(void* arg) {
    BlockProcrustes_t* blockProcrustes = (BlockProcrustes_t*) arg;

    // Allocate all required memory.
    double* asdBlock = calloc(PACKED_SIZE(blockProcrustes -> globalList -> numSamples), sizeof(double));
    IBS_t* ibsBlock = calloc(PACKED_SIZE(blockProcrustes -> globalList -> numSamples), sizeof(IBS_t));
    double** bootX = init_matrix(blockProcrustes -> globalList -> numSamples, blockProcrustes -> k);
    RealSymEigen_t* eigen = init_real_sym_eigen(blockProcrustes -> globalList -> numSamples);

    // The current block to operate on.
    Block_t* current = NULL;

    while (true) {

        // Reuse list lock to get the next block.
        pthread_mutex_lock(&listLock);
        if (*(blockProcrustes -> current) == NULL) {
            pthread_mutex_unlock(&listLock);
            break;
        }
        current = *(blockProcrustes -> current);
        *(blockProcrustes -> current) = (*(blockProcrustes -> current)) -> next;
        if (current -> isDropped)
            PROGRESS("Block on %s from %d to %d was dropped. Skipping.\n", current -> chrom, current -> startCoordinate, current -> endCoordinate);
        else
            PROGRESS("Performing MDS and Procrustes for block number %d on %s from %d to %d.\n", current -> blockNum, current -> chrom, current -> startCoordinate, current -> endCoordinate);
        // If we reached the end of the list.
        if (current -> next == NULL && blockProcrustes -> numReps == 0)
            PROGRESS("\nProcrustes finished ...\n\n");
        else if (current -> next == NULL)
            PROGRESS("\nProcrustes finished. Starting bootstrap ...\n\n");
        pthread_mutex_unlock(&listLock);

        if (!current -> isDropped) {
            // Convert to ASD first.
            for (int i = 0; i < blockProcrustes -> globalList -> numSamples; i++)
                for (int j = i + 1; j < blockProcrustes -> globalList -> numSamples; j++)
                    asdBlock[PACKED_INDEX(i, j)] = ibs_to_asd(current -> alleleCounts[PACKED_INDEX(i, j)]);

            // Perform MDS and Procrustes.
            double** X = init_matrix(eigen -> N, blockProcrustes -> k);
            current -> varCapt = compute_classical_mds(eigen, asdBlock, blockProcrustes -> k, X);
            // I am doing this to be safe. Results of cMDS are already centered.
            normalize_matrix(X, blockProcrustes -> globalList -> numSamples, blockProcrustes -> k);

            // In case MDS did not converge, treat the block as having no effect.
            if (current -> varCapt == -1) {
                current -> X = NULL;
                current -> procrustesT = 0;
                destroy_matrix(X, eigen -> N);
            // If we are comparing against the global.
            } else if (blockProcrustes -> y == NULL) {
                current -> X = X;
                current -> procrustesT = procrustes_statistic(current -> X, NULL, blockProcrustes -> globalList -> X, NULL, eigen, eigen -> N, blockProcrustes -> k, true);
            // If we are comparing against user defined coordinates.
            } else {
                current -> X = X;
                current -> procrustesT = procrustes_statistic(current -> X, NULL, blockProcrustes -> y, blockProcrustes -> y0, eigen, eigen -> N, blockProcrustes -> k, true);
            }
        }
    }

    // Create RNG. The seed is derived deterministically from the user's --seed in
    //  procrustes(), so a given (seed, thread count) reproduces exactly.
    //  The old code used time(NULL) + (long) arg -- wall clock plus the numeric value of
    //  a heap pointer -- which made results irreproducible run to run.
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, blockProcrustes -> seed);

    int N = blockProcrustes -> globalList -> numSamples;
    int packedSize = PACKED_SIZE(N);

    // Cap the number of consecutive non-converging replicates so a pathological input
    //  cannot spin here forever (the old `while (numReps != 0)` loop had no such bound).
    int consecutiveFailures = 0;
    const int maxConsecutiveFailures = 1000;

    // Start the bootstrap.
    while (blockProcrustes -> numReps != 0) {
        double bootStrappedT;

        // Claim a replicate slot BEFORE doing the work, so we do not compute up to
        //  NUM_THREADS - 1 replicates past the budget and throw them away.
        pthread_mutex_lock(&genomeLock);
        if (*(blockProcrustes -> currentReplicate) >= blockProcrustes -> numReps) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        pthread_mutex_unlock(&genomeLock);

        // Build the replicate: sum the IBS counts of sampleSize blocks drawn with
        //  replacement, then convert the SUM to ASD.
        //  The old code accumulated into ibsBlock and then never read it -- asdBlock was
        //  built from `temp -> alleleCounts`, the single block drawn on the last
        //  iteration. The null distribution was therefore the distribution of
        //  single-block statistics, and -s / --samples had no effect at all.
        memset(ibsBlock, 0, packedSize * sizeof(IBS_t));
        for (int s = 0; s < blockProcrustes -> sampleSize; s++) {
            Block_t* temp = get_random_block(r, blockProcrustes -> globalList);
            for (int p = 0; p < packedSize; p++)
                add_ibs(&(ibsBlock[p]), &(temp -> alleleCounts[p]));
        }
        for (int p = 0; p < packedSize; p++)
            asdBlock[p] = ibs_to_asd(ibsBlock[p]);

        // Perfrom MDS.
        double effectiveRank = compute_classical_mds(eigen, asdBlock, blockProcrustes -> k, bootX);

        // If MDS does not converge, then we do not count the sample.
        if (effectiveRank == -1) {
            if (++consecutiveFailures >= maxConsecutiveFailures) {
                fprintf(stderr, "MDS failed to converge for %d consecutive bootstrap replicates. Giving up on the bootstrap.\n", maxConsecutiveFailures);
                break;
            }
            continue;
        }
        consecutiveFailures = 0;

        normalize_matrix(bootX, N, blockProcrustes -> k);

        // Calculate our Procrustes statistic.
        if (blockProcrustes -> y == NULL) 
            bootStrappedT = procrustes_statistic(bootX, NULL, blockProcrustes -> globalList -> X, NULL, eigen, eigen -> N, blockProcrustes -> k, false);
        else 
            bootStrappedT = procrustes_statistic(bootX, NULL, blockProcrustes -> y, blockProcrustes -> y0, eigen, eigen -> N, blockProcrustes -> k, false);
        

        pthread_mutex_lock(&genomeLock);
        if (*(blockProcrustes -> currentReplicate) >= blockProcrustes -> numReps) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        blockProcrustes -> globalList -> samplingDistribution[*(blockProcrustes -> currentReplicate)] = bootStrappedT;
        *(blockProcrustes -> currentReplicate) += 1;
        PROGRESS("Completed Replicate %d of the Bootstrap.\n", *(blockProcrustes -> currentReplicate));
        pthread_mutex_unlock(&genomeLock);
    }

    // NOTE: bootX must be freed before eigen -- the old order read eigen -> N after
    //  destroy_real_sym_eigen() had freed it (AddressSanitizer: heap-use-after-free).
    destroy_matrix(bootX, eigen -> N);
    destroy_real_sym_eigen(eigen);
    free(asdBlock);
    free(ibsBlock);
    free(blockProcrustes);
    gsl_rng_free(r);

    return NULL;
}

// Allocate and populate the per-thread argument block for procrustes_bootstrap().
static BlockProcrustes_t* init_block_procrustes(BlockList_t* globalList, double** y, double* y0,
                                                int k, Block_t** current, int numReps,
                                                int sampleSize, int* currentReplicate,
                                                unsigned long seed) {
    BlockProcrustes_t* blockProcrustes = calloc(1, sizeof(BlockProcrustes_t));
    blockProcrustes -> globalList = globalList;
    blockProcrustes -> y = y;
    blockProcrustes -> y0 = y0;
    blockProcrustes -> k = k;
    blockProcrustes -> current = current;
    blockProcrustes -> numReps = numReps;
    blockProcrustes -> sampleSize = sampleSize;
    blockProcrustes -> currentReplicate = currentReplicate;
    blockProcrustes -> seed = seed;
    return blockProcrustes;
}

void procrustes(BlockList_t* globalList, double** y, double* y0, int k, int NUM_THREADS, int numReps, int sampleSize, unsigned long seed, bool quiet) {

    beQuiet = quiet;

    // The bootstrap needs at least one block it is allowed to sample.
    if (numReps > 0 && globalList -> numKept == 0) {
        fprintf(stderr, "Every block was dropped, so no bootstrap can be performed. Continuing without p-values.\n");
        numReps = 0;
    }
    if (numReps > 0 && sampleSize > 0 && globalList -> numKept < globalList -> numBlocks)
        PROGRESS("Bootstrap will resample %d blocks from the %d blocks that were kept.\n", sampleSize, globalList -> numKept);

    int* currentReplicate = calloc(1, sizeof(int));
    *currentReplicate = 0;
    // calloc(0, ...) may return NULL, and the p-value pass must not touch it in that case.
    globalList -> samplingDistribution = numReps > 0 ? calloc(numReps, sizeof(double)) : NULL;
    // NOTE: the old code did `Block_t* current = calloc(1, sizeof(Block_t*));` and then
    //  immediately overwrote the pointer (leaking the allocation), and later called
    //  free(current) on what by then was a Block_t from the list.
    Block_t* current = globalList -> head;

    // Compute Procrustes statistic for each block and bootstrap.
    //  Thread i draws from a stream seeded seed + i, so a given (--seed, -t) is exactly
    //  reproducible. Different -t gives a different partition of the work and hence a
    //  different (but equally valid) set of replicates.
    if (NUM_THREADS == 1) {
        procrustes_bootstrap((void*) init_block_procrustes(globalList, y, y0, k, &current,
                                                          numReps, sampleSize, currentReplicate, seed));
    } else {
        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, procrustes_bootstrap,
                           (void*) init_block_procrustes(globalList, y, y0, k, &current, numReps,
                                                         sampleSize, currentReplicate, seed + i + 1));
        procrustes_bootstrap((void*) init_block_procrustes(globalList, y, y0, k, &current, numReps,
                                                          sampleSize, currentReplicate, seed));
        // Wait for all threads to finish.
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);
    }

    // Use the number of replicates that actually completed, not the number requested;
    //  a replicate whose MDS did not converge was never written.
    int completedReps = *currentReplicate;
    globalList -> numReps = completedReps;
    free(currentReplicate);

    // Compute p-values. We are doing this the lazy way.
    if (completedReps > 0) {
        for (Block_t* temp = globalList -> head; temp != NULL; temp = temp -> next) {
            // Skip if window dropped.
            if (temp -> isDropped || temp -> X == NULL)
                continue;
            int numGreater = 1;
            for (int j = 0; j < completedReps; j++) {
                if (temp -> procrustesT <= globalList -> samplingDistribution[j])
                    numGreater++;
            }
            temp -> pvalue = numGreater / (double) (completedReps + 1);
        }
        // Global p-value.
        if (y != NULL) {
            int numGreater = 1;
            for (int j = 0; j < completedReps; j++) {
                if (globalList -> procrustesT <= globalList -> samplingDistribution[j])
                    numGreater++;
            }
            globalList -> pvalue = numGreater / (double) (completedReps + 1);
        }
    }
}