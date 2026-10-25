#!/bin/bash
# Example usage and roundtrip test for extraction from FLARE outputs.
set -euo pipefail

CHROM="21"
STEM="agd_250k_chr${CHROM}_flare"

time python -m tractor_prep extract_flare \
    "./test_flare/${STEM}.anc.vcf.gz" \
    --vcf-gzip \
    --out-gzip \
    --outdir="./test_flare/recon_${CHROM}" \
    --no-header-char \
    --num-ancestors=6


#### Test the outputs (the rest of the script is just for testing)
function test_output() {
    TRACT="./test_flare/tract_${CHROM}/${STEM}.anc.anc$1.$2.txt"
    RECON="./test_flare/recon_${CHROM}/${STEM}.anc$1.$2.txt"
    # Need to decompress the output files before they can be tested with diff
    pigz --keep --decompress --force "$TRACT.gz" "$RECON.gz"
    diff -q "$TRACT" "$RECON"
    rm "$TRACT" "$RECON"
}
# Check the tracts for differences (run all in paralell)
for idx in `seq 0 5`; do
    test_output "$idx" "dosage" &
    test_output "$idx" "hapcount" &
done
wait
