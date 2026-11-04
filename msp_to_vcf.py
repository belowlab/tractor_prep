"""
Utility script to convert RFMix style outputs to FLARE style outputs for a VCF
input with arbitrary FORMAT. I.e., it reads a file of ancestry windows and a
corresponding VCF, and adds the ancestry calls to the VCF sample format.

Uses PySAM to write arbitrarily formatted VCF files. (tractor_prep *requires*
non-FLARE VCFs to have only `GT` format, and thus does *not* require PySAM).
"""
import argparse
import builtins
import contextlib
from pathlib import Path

import numpy as np
import pysam
from tractor_prep import parse_msp, read_file


def edit_headers(vcf: pysam.VariantFile, subpops: list[(bytes, bytes)]):
    """Add ancestry information to VCF file headers"""
    # Need to add the FORMAT headers for AN1 and AN2 before adding to samples
    vcf.header.formats.add('AN1', 1, 'Integer', 'Ancestry of first haplotype')
    vcf.header.formats.add('AN2', 1, 'Integer', 'Ancestry of second haplotype')
    # Then the (purely informational) line with the list of ancestry codes
    ancs = str(b','.join([b'='.join(sp) for sp in subpops]), 'utf-8')
    vcf.header.add_line(f'##ANCESTRY=<{ancs}>')


# NOTE: all three of vcf, msp, and out may be backed by open files.
def edit_records(vcf, msp, out):
    # Set up the generator over the MSP windows. For compat reasons with pysam
    # we have to cast the numpy int8s to Python ints, and the chrom ID to str
    windows = (
        (
            str(idx[0], 'utf-8'),
            int(idx[1]), int(idx[2]),
            call_a.tolist(), call_b.tolist()
        )
        for idx, (call_a, call_b), *_ in iter(msp)
    )
    msp_chrom, start, stop, call_a, call_b = next(windows)

    for rec in vcf.fetch():
        # Align the VCF position until it is in an MSP window
        if rec.chrom == msp_chrom and rec.pos < start:
            continue
        skip_line = False
        while not (rec.chrom == msp_chrom and (start <= rec.pos < stop)):
            try:
                msp_chrom, start, stop, call_a, call_b = next(windows)
            except StopIteration:
                break
            if rec.chrom == msp_chrom and rec.pos < start:
                skip_line = True
                break
        if skip_line:
            continue

        sorted_samples = sorted(rec.samples.iteritems())
        for (_, sample), an1, an2 in zip(sorted_samples, call_a, call_b):
            sample['AN1'] = an1
            sample['AN2'] = an2

        out.write(rec)


def edit_vcf(out_path, vcf_path, msp_path,
             cohort_path=None,
             msp_open=builtins.open):
    """Rewrite a VCF with Ancestry information from an MSP file"""
    # Read bytes: remember that parse_msp / read_file operate on bytes
    if cohort_path:
        cohort = np.sort(cohort_path.read_bytes().strip().split(b'\n'))
    else:
        cohort = None

    msp = parse_msp(read_file(msp_path, opener=msp_open),
                    cohort=cohort, parse_subpops=True)
    subpops = next(msp)
    msp_ids = next(msp)[1]

    with pysam.VariantFile(vcf_path) as vcf:
        # NOTE: subset_samples *must* be called before reading any records.
        # Also it needs a python list of strings (not bytes).
        vcf.subset_samples(msp_ids.astype(str).tolist())
        edit_headers(vcf, subpops)
        with pysam.VariantFile(out_path, 'w', header=vcf.header) as out:
            edit_records(vcf, msp, out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--msp', type=Path, required=True)
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--cohort', '--iids', type=Path)
    args = parser.parse_args()
    edit_vcf(args.out, args.vcf, args.msp, arg.cohort)
