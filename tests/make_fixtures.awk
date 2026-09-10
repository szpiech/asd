# File: make_fixtures.awk
# Purpose: Emit a small deterministic phased VCF to stdout for the regression tests.
#
# Deterministic by construction: genotypes come from a linear congruential generator
#  seeded from -v seed=, so every platform and every awk implementation produces
#  byte-identical output. Do not substitute rand()/srand() here; they differ between
#  awk implementations.
#
# Variables:
#   n        number of samples (default 20)
#   seed     LCG seed (default 1)
#   layout   "uniform" | "sparse" | "missing"
#   nblocks  number of 1 Mb blocks to emit (default 6)
#   dense    records per dense block (default 400)
#   sparse   records per sparse block (default 5)
#   missfrac fraction of sites at which sample 0 is set to .|. (layout=missing)

function u() {
    state = (state * 1103515245 + 12345) % 2147483648;
    return state / 2147483648.0;
}

BEGIN {
    if (n == 0) n = 20;
    if (seed == 0) seed = 1;
    if (layout == "") layout = "uniform";
    if (nblocks == 0) nblocks = 6;
    if (dense == 0) dense = 400;
    if (sparse == 0) sparse = 5;
    state = seed;

    printf "##fileformat=VCFv4.2\n";
    printf "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (s = 0; s < n; s++) printf "\tS%d", s;
    printf "\n";

    for (b = 0; b < nblocks; b++) {
        # In the sparse layout every other block gets few records, INCLUDING the last
        #  one -- a dropped trailing block is what used to produce invalid JSON.
        nrec = dense;
        if (layout == "sparse" && b % 2 == 1) nrec = sparse;

        step = int(1000000 / (nrec + 1));
        for (r = 1; r <= nrec; r++) {
            pos = b * 1000000 + r * step;
            # Allele frequency per site, kept away from the MAF filter boundary.
            p = 0.2 + 0.6 * u();
            printf "1\t%d\t.\tA\tT\t.\tPASS\t.\tGT", pos;
            for (s = 0; s < n; s++) {
                if (layout == "missing" && s == 0 && u() < missfrac) {
                    printf "\t.|.";
                    # Keep the stream aligned with the non-missing case.
                    u(); u();
                    continue;
                }
                left  = (u() < p) ? 1 : 0;
                right = (u() < p) ? 1 : 0;
                printf "\t%d|%d", left, right;
            }
            printf "\n";
        }
    }
}
