"""
Convert FLARE output VCFs to Tractor input hapcount / dosage files.
"""
import builtins
import itertools
import os
import sys

import numpy as np
try:
    from isal.igzip_threaded import open as gzip_open
    HAS_ISAL = True
except ImportError:
    from gzip import open as gzip_open
    HAS_ISAL = False


def read_file(filename, opener=builtins.open):
    """Read the contents of a (maybe compressed) file or stdin"""
    # If we're being piped to or otherwise want to read stdin, read stdin
    if filename in ('-' or 'stdin') or not os.isatty(sys.stdin.fileno()):
        yield from sys.stdin.buffer
    else:
        with opener(filename, 'rb') as file:
            yield from file


def parse_flare_vcf(stream):
    ## This is how to make an array for checking the cell format:
    # seps = np.array(list(map(ord, b'|::\t')), dtype=np.int8)
    for i, line in enumerate(stream):
        # Skip the bcftools VCF headers
        if line.startswith(b'##'):
            continue

        # Emit the header with the 5 kept field names and IIDs
        elif line.startswith(b'#'):
            src_iids = np.array(line.strip().split(b'\t')[9:])
            yield (b'CHROM', b'POS', b'ID', b'REF', b'ALT'), src_iids
            continue

        # Do *not* strip the line: the final \n is required for reshaping.
        # Keep the first 5 fields; ditch QUAL, FILTER, INFO, and FORMAT
        *idx, _, _, _, _, calls = line.split(b'\t', 9)
        # `frombuffer` of mv before reshaping etc; `copy` after. The copy
        # *does* make all downstream operations significantly faster.
        calls = memoryview(calls)
        calls = np.frombuffer(calls, dtype=np.int8)
        calls = calls.reshape(-1, 8).T[::2]
        calls = np.copy(calls, order='c')
        np.subtract(calls, 48, out=calls)   # ord(b'0') == 48
        yield idx, calls


def extract_ancestors(stream, ancestors):
    ancestors = np.asarray(ancestors).reshape(-1, 1).astype(np.int8)

    # Emit the header line, and learn the size of the call grid from it
    idx, iids = next(stream)
    yield idx, iids

    call_grid = np.tile(ancestors, len(iids))
    for (idx, calls) in stream:
        # `calls` ought to have four rows (shape [4, n_iids]: see above)
        geno_a, geno_b, call_a, call_b = calls

        # Get the haplotype counts
        call_grid_a = (call_a == call_grid)
        call_grid_b = (call_b == call_grid)
        hap_cts = call_grid_a.view(np.int8) + call_grid_b.view(np.int8)

        # Get the dosage counts
        geno_grid_a = (geno_a & call_grid_a)
        geno_grid_b = (geno_b & call_grid_b)
        dos_cts = geno_grid_a.view(np.int8) + geno_grid_b.view(np.int8)

        yield idx, hap_cts, dos_cts


def write_tracts(stream, hap_files, dos_files):
    # Write the headers
    idx, iids = next(stream)
    header = b'\t'.join(itertools.chain(idx, iids))
    for file in itertools.chain(hap_files, dos_files):
        file.write(header)
        file.write(b'\n')

    # Make a line buffer
    sbuff = np.full(2 * len(iids), b'\t', dtype='S1')
    sbuff[-1] = b'\n'
    vtype = f'S{len(sbuff)}'
    n_anc = len(hap_files)

    for (idx, hap, dos) in stream:
        hap_ords = hap + 48
        dos_ords = dos + 48
        index = b'\t'.join(idx)
        for i in range(n_anc):
            hap_files[i].write(index)
            hap_files[i].write(b'\t')
            sbuff[::2] = hap_ords[i].view('S1')
            hap_files[i].write(sbuff.view(vtype)[0])

            dos_files[i].write(index)
            dos_files[i].write(b'\t')
            sbuff[::2] = dos_ords[i].view('S1')
            dos_files[i].write(sbuff.view(vtype)[0])


def cli_app():
    from argparse import ArgumentParser
    from pathlib import Path
    parser = ArgumentParser()
    f = parser.add_argument

    f('vcf', type=str, help='Path to FLARE ancestry VCF file')
    f('--outdir', type=Path, help='Directory in which to write output files')
    f('--stem', type=str, help='Output filename stem')
    f('--vcf-gzip', action='store_true', help='Input VCF is GZIP compressed')
    f('--out-gzip', action='store_true', help='GZIP compress outputs')
    f('--ancestors', nargs='*', type=int, help='Specific anc pop(s) to extract')
    f('--num-ancestors', type=int, help='Number of anc pop(s) to extract')
    f('--gzip-threads', type=int, default=0,
      help='Attempt to use multithreaded gzip (default: no threads)')

    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # A specific list of ancestors is preferred to a number range
    if args.ancestors:
        ancestors = np.array(args.ancestors)
    elif args.num_ancestors:
        ancestors = np.arange(args.num_ancestors)
    else:
        raise RuntimeError('--ancestors or --num-ancestors is required')

    # Use user-passed stem if given; otherwise strip *all* suffices to get stem
    if args.stem:
        stem = args.stem
    else:
        p0 = args.vcf
        while p0 != p0.with_suffix(''):
            p0 = p0.with_suffix('')
        stem = p0.stem

    # Define a function for opening gzip files at the right level & threading
    def gz_opener(name, mode, threads=args.gzip_threads):
        if HAS_ISAL:
            return gzip_open(name, mode, compresslevel=3, threads=threads)
        else:
            return gzip_open(name, mode, compresslevel=9)

    # Format the output file names and determine output opener
    hap_names = [f'{stem}.anc{i}.hapcount.txt' for i in ancestors]
    dos_names = [f'{stem}.anc{i}.dosage.txt' for i in ancestors]
    if args.out_gzip:
        hap_names = [f'{name}.gz' for name in hap_names]
        dos_names = [f'{name}.gz' for name in dos_names]
        open_output = gz_opener
    else:
        open_output = builtins.open

    # Determine input opener
    if args.vcf_gzip:
        open_input = gz_opener
    else:
        open_input = builtins.open

    #### RUN THE EXTRACTION
    stream = read_file(args.vcf, opener=open_input)
    stream = parse_flare_vcf(stream)
    stream = extract_ancestors(stream, ancestors)
    try:
        hap_files = [open_output(args.outdir/name, 'wb') for name in hap_names]
        dos_files = [open_output(args.outdir/name, 'wb') for name in dos_names]
        write_tracts(stream, hap_files, dos_files)
    finally:
        for fh in itertools.chain(hap_files, dos_files):
            fh.close


if __name__ == '__main__':
    cli_app()
