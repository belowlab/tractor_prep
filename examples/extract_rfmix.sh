#!/bin/bash
set -euo pipefail
# Example using RFMix outputs (VCFs and MSP files). The source VCFs in this
# example are assumed to have much more data than simply 'GT' format (e.g.,
# "GT:DS:HDS:GP"); the VCF files must be stripped before tract extraction.
# Since only the first five fields are kept by extraction, we can also strip
# QUAL, FILTER, and INFO fields if present to reduce the tractor_prep.py input
# file size.

CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
RSIDS="./pqtl_rsid.txt"
mkdir -p results/{stripped,tracts,unkinked}
mkdir -p timings

function run_extraction() {
    # The VCF files need to have annotations stripped. Here we write to an
    # intermediate file ($VCF); but we could also use an input pipe.
    # tractor_prep.py will read from stdin if the VCF positional argument is
    # "-" or "stdin":
    #
    # bcftools annotate -x FORMAT,INFO "./src_vcf/chr$1.vcf" -Ov \
    #   |\
    # { time python -m tractor_prep extract_rfmix "-" \
    #   # ...etc
    VCF="./results/stripped/chr$1.vcf"
    bcftools annotate -x FORMAT,INFO "./src_vcf/chr$1.vcf" -Ov -o "$VCF"
    { time python -m tractor_prep extract_rfmix "$VCF" \
        --msp="./src_msp/chr$1.msp" \
        --outdir "./results/tracts/" \
        --iids "$RSIDS" \
        --stem "chr$1" \
        --no-vcf-gzip \
        --out-gzip \
        --num-ancestors 4 \
        --unkink \
        --emit-msp \
        --emit-vcf ;
    } 2>"timings/chrom_$1.txt"
    # Move the unkinked files out of the tracts directory
    mv "./results/tracts/chr$1_switches.msp.gz" results/unkinked/
    mv "./results/tracts/chr$1_unkinked.msp.gz" results/unkinked/
    mv "./results/tracts/chr$1_unkinked.vcf.gz" results/unkinked/
}

# Run all 22 in parallel - each gets a single core.
for chrom in "${CHR[@]}"; do
    run_extraction $chrom &
done
# Wait for all to finish running before script exit.
wait
