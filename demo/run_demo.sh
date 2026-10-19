#!/usr/bin/env bash
# Run the filtering pipeline on the small synthetic dataset in demo/data/.
# Outputs are written to demo/output/ and compared against demo/expected_output/.
#
# Usage:  bash demo/run_demo.sh            (uses `python3` on PATH)
#         PYTHON=/path/to/python bash demo/run_demo.sh
set -euo pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Scripts live in the repository root (GitHub) or in scripts/ (submission zip).
if [[ -d "$DEMO_DIR/../scripts" ]]; then
    SCRIPTS="$DEMO_DIR/../scripts"
else
    SCRIPTS="$DEMO_DIR/.."
fi
DATA="$DEMO_DIR/data"
OUT="$DEMO_DIR/output"
EXPECTED="$DEMO_DIR/expected_output"
PYTHON="${PYTHON:-python3}"
# remove_clusters.py writes cluster members from a Python set; fix the hash seed so
# the member order (and therefore the output) is reproducible.
export PYTHONHASHSEED=0

rm -rf "$OUT"
mkdir -p "$OUT"
cp "$DATA/demo_dsa.fasta" "$OUT/demo_dsa.fasta"   # keep the .fai index out of demo/data

step() { echo; echo "=== $*"; }

# ---------------------------------------------------------------------------
step "Stage A: remove clustered variants (>=4 SNVs within 500 bp)"
"$PYTHON" "$SCRIPTS/inter_snp_distance.py" \
    -i "$DATA/candidates.in_good_regions.tsv" \
    -o "$OUT/A.intersnp_annotated.tsv"
"$PYTHON" "$SCRIPTS/get_clusters.py" \
    -i "$OUT/A.intersnp_annotated.tsv" \
    -o "$OUT/A.clusters.tsv"
"$PYTHON" "$SCRIPTS/remove_clusters.py" \
    -v "$OUT/A.intersnp_annotated.tsv" \
    -c "$OUT/A.clusters.tsv" \
    -o "$OUT/A.no_clustered_variants.tsv"

# ---------------------------------------------------------------------------
step "Stage B: remove variants with gnomAD v4.1 AF >= 1e-4"
"$PYTHON" "$SCRIPTS/filter_gnomad.py" \
    -v "$DATA/gnomad_annotated.tsv" \
    -o "$OUT/B.no_gnomad.tsv"

# ---------------------------------------------------------------------------
step "Stage C1: combine HiFi + Illumina read support within each tissue"
for tissue in blood liver lung; do
    "$PYTHON" "$SCRIPTS/single_tissue_combine_platform.py" \
        -v "$DATA/variants.tsv" \
        -hifi "$DATA/$tissue.hifi.read_counts.tsv" \
        -ilm "$DATA/$tissue.illumina.read_counts.tsv" \
        -o "$OUT/C1.$tissue.tsv"
done

step "Stage C2: combine evidence across tissues (phasing, binomial and strand-bias tests)"
"$PYTHON" "$SCRIPTS/combine_tissues.py" \
    -v "$DATA/variants.tsv" \
    -t blood "$OUT/C1.blood.tsv" \
    -t liver "$OUT/C1.liver.tsv" \
    -t lung "$OUT/C1.lung.tsv" \
    -o "$OUT/C2.combined_tissues.tsv"

step "Stage C3: drop variants failing the binomial or strand-bias tests (same logic as remove_tests.sh)"
awk -F'\t' '$18 != "False" && $19 != "False"' "$OUT/C2.combined_tissues.tsv" > "$OUT/C3.passed_binomial_test.tsv"
awk -F'\t' '$20 != "FAIL"' "$OUT/C3.passed_binomial_test.tsv" > "$OUT/C3.passed_strand_bias.tsv"

step "Stage C4: remove GA-rich / homopolymer-adjacent variants"
"$PYTHON" "$SCRIPTS/ga_richness_filter.py" \
    -v "$OUT/C3.passed_strand_bias.tsv" \
    -f "$OUT/demo_dsa.fasta" \
    -o "$OUT/C4.no_homopolymer.tsv"

# ---------------------------------------------------------------------------
if command -v bedtools >/dev/null 2>&1 && "$PYTHON" -c "import pybedtools" 2>/dev/null; then
    step "Stage D (optional): remove variants in NucFlag-flagged and distorted regions"
    "$PYTHON" "$SCRIPTS/remove_flagged_regions.py" \
        -i "$DATA/candidates.in_good_regions.tsv" \
        -b "$DATA/nucflag.bed" \
        -o "$OUT/D.no_flagged_regions.tsv"
    "$PYTHON" "$SCRIPTS/remove_distorted_regions.py" \
        -i "$DATA/candidates.in_good_regions.tsv" \
        -b1 "$DATA/distortions_hap1.bed" \
        -b2 "$DATA/distortions_hap2.bed" \
        -o "$OUT/D.no_distortions.tsv"
else
    step "Stage D (optional): skipped (requires bedtools on PATH and pybedtools)"
fi

# ---------------------------------------------------------------------------
step "Comparing outputs with demo/expected_output/"
status=0
for expected in "$EXPECTED"/*.tsv; do
    name="$(basename "$expected")"
    if [[ ! -f "$OUT/$name" ]]; then
        echo "SKIPPED  $name"
    elif diff -q "$expected" "$OUT/$name" >/dev/null; then
        echo "OK       $name"
    else
        echo "MISMATCH $name"
        status=1
    fi
done

echo
if [[ $status -eq 0 ]]; then
    echo "Demo finished successfully. Final call set: $OUT/C4.no_homopolymer.tsv"
else
    echo "Demo finished, but some outputs differ from expected_output/." >&2
fi
exit $status
