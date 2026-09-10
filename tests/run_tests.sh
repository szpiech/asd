#!/bin/sh
# File: run_tests.sh
# Purpose: Regression tests for LODESTAR. Run with `make test` from the repository root.
#
# Every test here corresponds to a defect that was found in v1.0. Dependencies are awk,
# gzip and a POSIX shell only -- no python, no jq -- so this runs anywhere the tool does.

set -u

LODESTAR=${LODESTAR:-./bin/lodestar}
DIR=$(dirname "$0")
WORK=$(mktemp -d "${TMPDIR:-/tmp}/lodestar_tests.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

PASS=0
FAIL=0

ok()   { PASS=$((PASS + 1)); printf '  ok   %s\n' "$1"; }
bad()  { FAIL=$((FAIL + 1)); printf '  FAIL %s\n' "$1"; }
check() { if [ "$1" = "0" ]; then ok "$2"; else bad "$2"; fi; }

mkfix() {
    # mkfix <output> <awk -v assignments...>
    out=$1; shift
    awk "$@" -f "$DIR/make_fixtures.awk" | gzip -c > "$WORK/$out"
}

printf '\nLODESTAR regression tests\n'
printf -- '-------------------------\n'

if [ ! -x "$LODESTAR" ]; then
    printf 'Cannot find an executable at %s. Run make first.\n' "$LODESTAR"
    exit 1
fi

printf '\nBuilding fixtures ...\n'
mkfix uniform.vcf.gz -v n=30 -v seed=7  -v layout=uniform -v nblocks=6 -v dense=400
mkfix sparse.vcf.gz  -v n=20 -v seed=11 -v layout=sparse  -v nblocks=6 -v dense=400 -v sparse=5
mkfix nomiss.vcf.gz  -v n=20 -v seed=13 -v layout=uniform -v nblocks=2 -v dense=1500
mkfix miss.vcf.gz    -v n=20 -v seed=13 -v layout=missing -v nblocks=2 -v dense=1500 -v missfrac=0.5

# ---------------------------------------------------------------- interface
printf '\nInterface\n'

"$LODESTAR" --help > "$WORK/help.txt" 2>&1
if grep -q 'Usage: lodestar' "$WORK/help.txt"; then ok "--help prints usage"; else bad "--help prints usage"; fi

"$LODESTAR" -h > "$WORK/help2.txt" 2>&1
if grep -q 'Usage: lodestar' "$WORK/help2.txt"; then ok "-h prints usage"; else bad "-h prints usage"; fi

"$LODESTAR" --version 2>&1 | grep -q 'LODESTAR v'
check $? "--version prints a version"

"$LODESTAR" -q -r 0 -o "$WORK/nosuch" "$WORK/does_not_exist.vcf.gz" > /dev/null 2>&1
if [ $? -ne 0 ]; then ok "missing input exits nonzero"; else bad "missing input exits nonzero"; fi

# ---------------------------------------------------------------- basic run
printf '\nBasic run\n'

"$LODESTAR" -q -b 1000000 -r 0 -o "$WORK/basic" "$WORK/uniform.vcf.gz" > /dev/null 2>&1
check $? "runs without a bootstrap"

[ -f "$WORK/basic.tsv" ] && [ -f "$WORK/basic.json" ]
check $? "writes both .tsv and .json"

# The GLOBAL row's NumLoci must equal the sum over the blocks. The old code's
#  blockNumOnChrom loop stopped at the tail, so the last block was never counted.
awk -F'\t' '
    /^#/ || /^BlockNum/ { next }
    $3 == "GLOBAL" { global = $6; next }
    { sum += $6 }
    END { exit (sum == global) ? 0 : 1 }
' "$WORK/basic.tsv"
check $? "GLOBAL NumLoci equals the sum over blocks"

# Non-computed statistics must be NA, not the -1 sentinel that reads as data.
awk -F'\t' '/^[0-9]/ && $10 == "-1.000000" { found = 1 } END { exit found ? 1 : 0 }' "$WORK/basic.tsv"
check $? "no -1 sentinels in the TSV"

grep -q '"Samples": \[ "S0"' "$WORK/basic.json"
check $? "sample identifiers appear in the JSON"

grep -q '"Seed":' "$WORK/basic.json"
check $? "effective seed is echoed into the JSON"

# ---------------------------------------------------------------- reproducibility
printf '\nReproducibility\n'

"$LODESTAR" -q -b 1000000 -r 200 --seed 99 -t 1 -o "$WORK/seedA" "$WORK/uniform.vcf.gz" > /dev/null 2>&1
"$LODESTAR" -q -b 1000000 -r 200 --seed 99 -t 1 -o "$WORK/seedB" "$WORK/uniform.vcf.gz" > /dev/null 2>&1
grep -v '^#' "$WORK/seedA.tsv" > "$WORK/a.txt"
grep -v '^#' "$WORK/seedB.tsv" > "$WORK/b.txt"
cmp -s "$WORK/a.txt" "$WORK/b.txt"
check $? "same seed gives identical p-values"

"$LODESTAR" -q -b 1000000 -r 200 --seed 12345 -t 1 -o "$WORK/seedC" "$WORK/uniform.vcf.gz" > /dev/null 2>&1
grep -v '^#' "$WORK/seedC.tsv" > "$WORK/c.txt"
if cmp -s "$WORK/a.txt" "$WORK/c.txt"; then
    bad "a different seed gives a different bootstrap"
else
    ok "a different seed gives a different bootstrap"
fi

# ---------------------------------------------------------------- bootstrap pooling
printf '\nBootstrap\n'

# -s controls how many blocks are pooled per replicate, so the width of the null must
#  depend on it. In v1.0 the pooled counts were discarded and -s had no effect at all.
spread() {
    awk -F'\t' '
        /^#/ || /^BlockNum/ || $3 == "GLOBAL" { next }
        $9 != "NA" { s += $9; ss += $9 * $9; n++ }
        END { if (n < 2) { print 0; exit } m = s / n; print (ss / n - m * m) }
    ' "$1"
}
"$LODESTAR" -q -b 1000000 -r 400 -s 1 --seed 5 -t 1 -o "$WORK/s1" "$WORK/uniform.vcf.gz" > /dev/null 2>&1
"$LODESTAR" -q -b 1000000 -r 400 -s 6 --seed 5 -t 1 -o "$WORK/s6" "$WORK/uniform.vcf.gz" > /dev/null 2>&1
awk -F'\t' 'FNR==NR { if ($3 != "GLOBAL" && $10 != "NA" && /^[0-9]/) a[++i] = $10; next }
            { if ($3 != "GLOBAL" && $10 != "NA" && /^[0-9]/) b[++j] = $10 }
            END { for (k = 1; k <= i; k++) if (a[k] != b[k]) { differ = 1 } exit differ ? 0 : 1 }' \
    "$WORK/s1.tsv" "$WORK/s6.tsv"
check $? "-s changes the null distribution (blocks are pooled)"

# ---------------------------------------------------------------- dropped blocks
printf '\nDropped blocks\n'

# v1.0 segfaulted here: get_random_block() never reset its cursor on a rejected draw.
"$LODESTAR" -q -b 1000000 -r 300 -d 100 --seed 3 -t 1 -o "$WORK/drop" "$WORK/sparse.vcf.gz" > /dev/null 2>&1
check $? "bootstrap with dropped blocks does not crash"

"$LODESTAR" -q -b 1000000 -r 300 -d 100 --seed 3 -t 4 -o "$WORK/dropt" "$WORK/sparse.vcf.gz" > /dev/null 2>&1
check $? "the same, multithreaded"

# The sparse fixture drops its LAST block, which used to leave a trailing comma.
if tr -d ' \t\n' < "$WORK/drop.json" | grep -q -e ',]' -e ',}'; then
    bad "JSON has no trailing comma when the last block is dropped"
else
    ok "JSON has no trailing comma when the last block is dropped"
fi

grep -q '"NumberOfBlocksDropped":' "$WORK/drop.json"
check $? "dropped-block count is reported"

# Every block dropped: report and continue rather than spin or crash.
"$LODESTAR" -q -b 1000000 -r 100 -d 100000 --seed 3 -o "$WORK/alldrop" "$WORK/sparse.vcf.gz" > /dev/null 2>&1
check $? "all blocks dropped does not crash"

# ---------------------------------------------------------------- target file
printf '\nTarget coordinates (-y)\n'

# 20 samples in nomiss.vcf.gz; hand it 25 rows.
awk 'BEGIN { for (i = 0; i < 25; i++) print "0.1\t0.2" }' > "$WORK/too_many.tsv"
"$LODESTAR" -q -r 0 -y "$WORK/too_many.tsv" -o "$WORK/ty" "$WORK/nomiss.vcf.gz" > /dev/null 2>&1
if [ $? -ne 0 ]; then ok "-y with too many rows is rejected"; else bad "-y with too many rows is rejected"; fi

awk 'BEGIN { for (i = 0; i < 20; i++) print "0.1" }' > "$WORK/too_narrow.tsv"
"$LODESTAR" -q -r 0 -y "$WORK/too_narrow.tsv" -o "$WORK/tn" "$WORK/nomiss.vcf.gz" > /dev/null 2>&1
if [ $? -ne 0 ]; then ok "-y with too few columns is rejected"; else bad "-y with too few columns is rejected"; fi

awk 'BEGIN { for (i = 0; i < 20; i++) print "not_a_number\t0.2" }' > "$WORK/nonnumeric.tsv"
"$LODESTAR" -q -r 0 -y "$WORK/nonnumeric.tsv" -o "$WORK/tnn" "$WORK/nomiss.vcf.gz" > /dev/null 2>&1
if [ $? -ne 0 ]; then ok "-y with a non-numeric entry is rejected"; else bad "-y with a non-numeric entry is rejected"; fi

awk 'BEGIN { for (i = 0; i < 20; i++) printf "%f\t%f\n", i / 20.0, (i % 5) / 5.0 }' > "$WORK/good.tsv"
"$LODESTAR" -q -b 1000000 -r 100 --seed 1 -y "$WORK/good.tsv" -o "$WORK/tg" "$WORK/nomiss.vcf.gz" > /dev/null 2>&1
check $? "-y with a well-formed matrix is accepted"

"$LODESTAR" -q -r 0 --geo -o "$WORK/geo" "$WORK/nomiss.vcf.gz" > /dev/null 2>&1
if [ $? -ne 0 ]; then ok "--geo without -y is rejected"; else bad "--geo without -y is rejected"; fi

# ---------------------------------------------------------------- missing genotypes
printf '\nMissing genotypes\n'

# In v1.0 a sample missing at half its sites became the dominant MDS outlier, because
#  the pair guard was an OR and ASD was not divided by the loci actually compared.
#  Test: sample S0's distance from the centroid must stay within the range set by the
#  other samples, and must not blow up relative to the no-missing-data run.
outlier_ratio() {
    awk '
        /"GlobalX"/ { inX = 1; next }
        inX && /^\]/ { inX = 0; next }
        inX {
            gsub(/[\[\],]/, " ");
            n++; x[n] = $1; y[n] = $2;
            sx += $1; sy += $2;
        }
        END {
            cx = sx / n; cy = sy / n;
            for (i = 1; i <= n; i++) {
                d = sqrt((x[i] - cx) ^ 2 + (y[i] - cy) ^ 2);
                if (i == 1) d0 = d; else if (d > dmax) dmax = d;
            }
            printf "%.4f\n", (dmax > 0) ? d0 / dmax : 0;
        }
    ' "$1"
}
"$LODESTAR" -q -b 10000000 -r 0 --afMissing 0.9 -o "$WORK/nomiss" "$WORK/nomiss.vcf.gz" > /dev/null 2>&1
"$LODESTAR" -q -b 10000000 -r 0 --afMissing 0.9 -o "$WORK/miss"   "$WORK/miss.vcf.gz"   > /dev/null 2>&1
R_NOMISS=$(outlier_ratio "$WORK/nomiss.json")
R_MISS=$(outlier_ratio "$WORK/miss.json")
printf '       S0 distance / max other: no missing = %s, 50%% missing = %s\n' "$R_NOMISS" "$R_MISS"
awk -v r="$R_MISS" 'BEGIN { exit (r <= 1.2) ? 0 : 1 }'
check $? "a half-missing sample is not the dominant outlier"

# ---------------------------------------------------------------- summary
printf -- '\n-------------------------\n'
printf '%d passed, %d failed\n\n' "$PASS" "$FAIL"
[ "$FAIL" -eq 0 ]
