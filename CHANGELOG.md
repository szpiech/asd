# Changelog

## v1.0.1

A code review of v1.0 found defects in three classes: results that were wrong, memory
errors, and interface problems. This release fixes them. **Results from v1.0 should be
regenerated** — items 1 and 2 change reported statistics.

### Results

1. **The block bootstrap now pools blocks.** `procrustes_bootstrap()` accumulated the IBS
   counts of `sampleSize` resampled blocks into `ibsBlock` and then never read it; the
   matrix passed to cMDS was built from the single block drawn on the last iteration. The
   null distribution was therefore the distribution of single-block Procrustes
   statistics — much wider than intended, giving conservative p-values — and `-s` had no
   effect. `tests/run_tests.sh` now asserts that `-s` changes the null.

2. **Missing genotypes are handled correctly.** The pair guard in the innermost IBS loop
   was `g[i].left != MISSING || g[j].right != MISSING` — an OR, comparing sample *i*'s
   left haplotype against sample *j*'s right. Pairs involving a missing call were still
   counted and `num_shared_alleles(MISSING, x)` returned 0, so missingness was scored as
   maximal dissimilarity. Separately, `ibs_to_asd()` computed `L`, the number of loci at
   which the pair was jointly called, but did not divide by it. Both are fixed: a pair
   contributes only when all four haplotypes are called, and the distance is per-locus.
   With one sample missing at 50% of sites, its distance from the MDS centroid relative to
   the largest among the other samples was 2.93 before and 0.76 after (0.72 with no
   missing data at all). On data with no missing genotypes the global configuration is
   unchanged to machine precision.

3. **`GlobalNumberOfLoci` and `GlobalNumberOfHaplotypes` no longer undercount.** The loop
   that assigned `blockNumOnChrom` stopped at `temp -> next != NULL`, so the last block was
   never added to the genome-wide totals.

4. **`relabel_haplotypes()` relabels both haplotypes.** It skipped a sample entirely when
   `.left == MISSING`, but `add_locus()` sets `.left` and `.right` independently, so a
   half-missing genotype kept a stale `.right` label while every other sample was
   renumbered. Only reachable with `-H > 1`.

5. **`compute_classical_mds()` returns `-1` on a negative input distance.** It returned
   `1`, which every caller read as "converged, 100% variance captured" over a matrix of
   zeros.

6. **p-value denominators use the number of replicates that completed**, not the number
   requested, and a replicate whose MDS fails to converge no longer counts.

### Crashes and memory errors

7. **No more segfault when a block is dropped and a bootstrap is requested.**
   `get_random_block()` initialised its cursor outside the rejection loop, so a rejected
   draw resumed the walk from wherever the previous one stopped and eventually ran off the
   end of the list. Blocks are now held in an array and drawn in O(1); the draw is also
   uniform, which it previously was not. All blocks dropped is reported and handled.

8. **Heap use-after-free in `procrustes_bootstrap()`** — `destroy_matrix(bootX, eigen -> N)`
   ran after `destroy_real_sym_eigen(eigen)`. Fired on every bootstrap run.

9. **Heap buffer overflow in `open_target_file()`** — row `N` of a `-y` file was written
   before the `numLines == N` bound was tested. `-y` input is now validated for row count,
   column count, and numeric content, and `fopen()` is checked.

10. **Parser hardening** — `gzopen()` returning `NULL` was passed to `gzerror()`; the
    metadata loop never terminated on a file with no `#CHROM` line; a record with more
    genotype columns than the header declared wrote past `nextLocus`; `alleleCounts` was
    one slot short of the index used for the missing-allele count, and the
    `numAlleles > MAX_NUM_ALLELES` test ran after the counting loop. A warning is issued
    once if `FORMAT` does not begin with `GT`, which the parser assumes.

11. **`init_matrix()` allocated `n * sizeof(double)` for a `double**` array** — correct only
    because the two are equal on LP64.

12. Leaks fixed: the `strdup`'d VCF header, the `cmd` string, a stray `calloc` in
    `procrustes()` that was overwritten on the next line (and later `free`d as if it were a
    live `Block_t`), and error paths that called `free(config)` instead of
    `destroy_lodestar_config()`. `make asan` (`-fsanitize=address,undefined`) is clean
    across a plain run, a multithreaded run, dropped blocks, all-blocks-dropped, a
    malformed `-y`, a well-formed `-y`, and `-H 5`.

### Interface

13. **`-h` and `--help` print the help menu**, and `-v`/`--version` prints the version.
    Haplotype size moved from `-h` to **`-H` / `--hap-size`**; `--haplotypeSize` still works.
14. **`--seed INT` (default 42).** The RNG was seeded `time(NULL) + (long) arg` — wall clock
    plus the numeric value of a heap pointer — so no run was reproducible. Thread *i* now
    draws from `seed + i`, and the effective seed is echoed into both output files.
    A given `(--seed, -t)` pair reproduces exactly.
15. **`-q` / `--quiet`** suppresses the per-block and per-replicate progress lines, which
    were unconditional and unbuffered (2,500 lines for a 1,500-block, 1,000-replicate run).
16. **`--boot-blocks`** is the new name for `--samples`, which selected the bootstrap
    resample size while meaning "select samples" in bcftools and plink. `--samples` is kept
    as an alias. `-o`/`--out`, `-b`/`--block-size`, `-k`/`--dimension`, `-y`/`--target`,
    `-r`/`--replicates`, `-d`/`--drop` all have long forms now.
17. **`--threads 0`** uses every available core.
18. **Sample identifiers appear in the JSON** as `"Samples"`, in the row order of every
    coordinate matrix. They were parsed and freed but never emitted, so every matrix was an
    unlabelled block in VCF column order and a mis-ordered `-y` file could not be detected.
19. **Invalid options exit non-zero.** A configuration error returned 0 to the shell.
20. **Output basename derivation** uses the last `.vcf` in the final path component, and
    handles an input with no `.vcf` at all. `strstr()` found the *first* occurrence, so
    `/data/run.vcf.d/chr1.vcf.gz` produced the basename `/data/run`, and a name without
    `.vcf` passed a large negative value as the `size_t` argument of `strndup()`.
21. **Validation messages match the checks.** `-r 0` and `-s 0` are legal and documented as
    such; `-h` is no longer compared against `-b` (loci versus base pairs); `-k >= 1` is
    enforced; `--geo` without `-y` is an error rather than silently ignored; `b:` appeared
    twice in the option string.

### Output format

22. **`NA` replaces the `-1` "not computed" sentinel in the TSV**, and `null` in the JSON.
    `-1` is read as data by `read.delim` and `read_csv`.
23. **The JSON is always valid.** A dropped trailing block left a trailing comma before the
    closing bracket, which `jsonlite::fromJSON` — and therefore `lodestarPlots.R` — could
    not read. The command line is now escaped, so a filename containing `"` or `\` no longer
    breaks the file.
24. The TSV header and the JSON record the version, effective seed, thread count,
    completed replicate count, block count, and number of blocks dropped.
25. `GlobalVarainceCaptured` is now also spelled `GlobalVarianceCaptured`. Both keys are
    emitted this release; the misspelling will be removed in the next.

### Performance

26. **`-O2`.** `CFLAGS` was `-c -Wall -g` with no optimisation flag, so every C source and
    all 99 vendored LAPACK/BLAS Fortran files were compiled unoptimised.
27. **`CC` and `FC` are now `=` rather than `?=`.** GNU make defines `CC = cc` and
    `FC = f77` as built-in defaults, and `?=` does not override a built-in default, so a
    bare `make` invoked `f77` (absent on most systems, failing the build) and compiled the C
    with `cc` — clang on macOS, contradicting the gfortran requirement documented in
    `MatrixOperations.c`.
28. **The redundant `globalCounts` update is gone** from the innermost loop. The genome-wide
    counts are by construction the sum of the per-block counts, so they are accumulated once
    per block: `O(numBlocks * N^2)` rather than `O(numHaplotypes * N^2)`.
29. **The IBS loop is unit-stride.** `PACKED_INDEX(i, j) = i + j(j+1)/2`, so the previous
    i-outer/j-inner order strode by `j+1` through an array that exceeds L2 at moderate `N`.
30. **The innermost work is branchless** — `increment_ibs_value()` indexes an array instead
    of switching on unpredictable data, and `num_shared_alleles()` is
    `max((l1==l2)+(r1==r2), (l1==r2)+(r1==l2))`.
31. **Blocks are stored in an array.** This removes the `O(numBlocks)` linked-list walk per
    bootstrap draw and the `O(numBlocks^2)` selection sort that restored block order after
    the threaded pass.
32. **Bootstrap workers claim a replicate slot before computing one**, so up to
    `threads - 1` replicates are no longer computed and discarded.

End-to-end, against v1.0 built with its own makefile:

| case | v1.0 | v1.0.1 | speedup |
|---|---|---|---|
| N=100, no bootstrap | 3.85 s | 0.90 s | 4.3x |
| N=100, 1000 replicates | 6.62 s | 1.37 s | 4.9x |
| N=100, 1000 replicates, 8 threads | 1.40 s | 0.43 s | 3.3x |
| N=300, no bootstrap | 31.90 s | 6.51 s | 4.9x |
| N=300, 1000 replicates | 77.42 s | 15.33 s | 5.1x |
| N=300, 1000 replicates, 8 threads | 17.86 s | 3.40 s | 5.3x |

100,000 records over 50 Mb, 50 blocks, on an 8-core machine. The fixed bootstrap does
strictly more work than v1.0's (it actually pools blocks) and is still ~5x faster.

### Build and tests

33. **`tests/run_tests.sh`** (`make test`) — 24 regression tests, one per defect above,
    depending only on a POSIX shell, `awk` and `gzip`. Fixtures are generated by a
    deterministic LCG so they are byte-identical on every platform. Run against v1.0 the
    suite reports 16 failures.
34. **`make asan`** builds with the address and undefined-behaviour sanitizers;
    **`make install`** installs to `$(PREFIX)/bin`.
35. `make clean` no longer fails on an already-clean tree; `lib/kstring.o` is built to its
    own target name rather than `src/kstring.o` (so it stopped rebuilding every time, and
    is now named explicitly on the link line); the LAPACK rule compiles each source to an
    explicit output instead of building in the repository root and `mv`-ing, which raced
    under `make -j`; object rules list their own sources and headers as prerequisites, so
    editing a `.c` or `.h` actually triggers a rebuild.
36. One local patch to vendored `lib/kstring.c`: `s->s + s->l` is undefined behaviour when
    `s->s` is `NULL`, even at zero offset.

### R plotting script

37. Blocks whose statistics are `NA` are filtered (the `-1` test is kept so v1.0 JSON still
    loads).
38. `as.numeric(gsub("chr", "", CHR))` turned every non-numeric chromosome name —
    scaffolds, `X`, `Y`, Drosophila-style `2L`/`3R` — into `NA`, and those blocks vanished
    from the plot with no warning. Names are now mapped to positions in order of first
    appearance and the original name is used for the axis label.
39. The two-colour chromosome palette is sized to the data rather than hard-coded to 44.
40. `Varaince` corrected in a plot title.
