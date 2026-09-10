
# LODESTAR v1.0.1

LODESTAR partitions a genome into blocks, computes a multidimensional scaling of allele
sharing within each block, and compares each block's configuration to the genome-wide one
(or to user-supplied coordinates) with a Procrustes statistic. Significance comes from a
block bootstrap.

## Manual

The manual can be found at [https://github.com/TQ-Smith/LODESTAR/wiki/LODESTAR-Manual]

## TL;DR

```
conda install bioconda::lodestar
lodestar --help
```

## Building from source

Requires a C compiler, `gfortran` (the vendored netlib LAPACK is Fortran), and zlib.

```
make                 # builds bin/lodestar with -O2
make test            # runs the regression suite (needs only sh, awk, gzip)
make asan            # builds bin/lodestar_asan with -fsanitize=address,undefined
make install         # installs to $(PREFIX)/bin, PREFIX defaults to /usr/local
make CC=clang        # override the compiler as usual
```

## Worked example

```
cd examples
../bin/lodestar -b 1000000 -r 1000 --seed 1 -t 4 -o AR_ST AR_ST.vcf.gz
```

writes two files:

- `AR_ST.tsv` — one row per block: coordinates on the genome, number of loci and
  haplotypes, variance captured by the MDS, the Procrustes statistic, and its bootstrap
  p-value. A final `GLOBAL` row carries the genome-wide values. `NA` marks a statistic
  that was not computed. The header comment records the version, the effective seed, the
  thread count, and how many blocks were dropped.
- `AR_ST.json` — the same summary plus every block's `N`-by-`k` coordinate matrix, the
  genome-wide matrix `GlobalX`, and the sample identifiers in `Samples` (the row order of
  every matrix).

To compare against known coordinates, pass an `N`-by-`k` TSV with `-y`, one row per VCF
sample **in VCF column order**:

```
../bin/lodestar -b 1000000 -r 1000 --seed 1 -y EUR_40_target.tsv EUR_40_chr6.vcf.gz
```

Plots are produced from the JSON by `lodestarPlots.R`:

```
Rscript lodestarPlots.R var    AR_ST.json     # variance captured along the genome
Rscript lodestarPlots.R pvals  AR_ST.json     # -log10(p) along the genome
Rscript lodestarPlots.R tvals  AR_ST.json     # Procrustes statistic along the genome
Rscript lodestarPlots.R mds    AR_ST.json     # the MDS of a chosen block
```

## Reproducibility

The bootstrap is seeded with `--seed` (default 42) and thread *i* draws from `seed + i`, so
a given `(--seed, -t)` pair reproduces exactly. Changing `-t` repartitions the work and
gives a different — equally valid — set of replicates. The effective seed is written into
both output files.

## Input assumptions

- gzipped VCF, **phased** genotypes (`|`).
- `GT` is the first subfield of `FORMAT`. A warning is issued once if it is not.
- Sites are filtered by `--maf` (default 0.05) and `--afMissing` (default 0.1). Note that
  `--afMissing` filters per *site*: a consistently poorly-covered *sample* passes every
  site filter, so check per-sample call rates yourself.

## Bundled dependencies

`lib/` vendors third-party sources so the build is self-contained:

| path | upstream | license |
|---|---|---|
| `lib/lapack` | netlib LAPACK/BLAS (reference) | BSD 3-clause |
| `lib/gsl` | GNU Scientific Library, RNG subset (`mt19937`) | GPL-3.0 |
| `lib/khash.h`, `lib/kseq.h`, `lib/kstring.*`, `lib/ketopt.h` | klib | MIT |

The GSL RNG subset is GPL-3.0, which is more restrictive than the MIT license below; a
distribution that bundles it must comply with the GPL.

## License 

This project is licensed under the [MIT License](./LICENSE), excepting the bundled
third-party sources listed above, which carry their own licenses.

## Changelog

See [CHANGELOG.md](./CHANGELOG.md).

## Please Cite