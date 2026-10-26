"""
Convert RFMix or FLARE output files to ancestry specific hapcount and dosage
files for downstream analysis by Tractor.
"""
import builtins
import contextlib
import itertools
import sys
from argparse import ArgumentParser, BooleanOptionalAction
from pathlib import Path
from pprint import pprint

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
    if (filename == '-') or (filename == 'stdin'):
        yield from sys.stdin.buffer
    else:
        with opener(filename, 'rb') as file:
            yield from file


def bytes_to_ords(bts, *, cols=Ellipsis, rows=4, stride=2):
    """Process a line of formatted chars into shaped single byte integers.

    E.g., if `bts` is something like `b'0|0\t0|1\t...'`
    The output with default kwargs is an ndarray: `[[0, 0, ...], [0, 1, ...]]`
    """
    # `frombuffer` of memoryview before reshaping etc; `copy` after. The copy
    # *does* make all downstream operations significantly faster.
    bts = memoryview(bts)
    bts = np.frombuffer(bts, dtype=np.int8)
    # .squeeze() In case the rows/stride combination selects a single row.
    bts = bts.reshape(-1, rows).T[::stride, cols].squeeze()
    bts = np.copy(bts, order='c')
    # ord(b'0') == 48: so if we subtract 48, then we cast from b'0' to 0, etc.
    np.subtract(bts, 48, bts)
    return bts


def process_iid_header(line, splitfn, cohort=None):
    """Process the header line into fields, iids, and iid cohort col indices"""
    index, iids = splitfn(line.strip().split(b'\t'))
    if cohort is not None:
        _, cols, _ = np.intersect1d(iids, cohort, return_indices=True)
        iids = np.copy(iids[cols], order='c')
    else:
        cols = Ellipsis
    return index, iids, cols


def parse_msp(stream, cohort=None):
    """Parse an iterator of newline-terminated byte as MSP window data"""
    def header_fn(parts):
        iids = np.array([x.strip().split(b'.')[0] for x in parts[6::2]])
        return tuple(parts[:6]), iids

    for line in (next(stream), next(stream)):
        if line.startswith(b'#Subpopulation'):
            pass
        elif line.startswith(b'#chm'):
            index, iids, cols = process_iid_header(line, header_fn, cohort)
            yield index, iids
            break

    for line in stream:
        # Do *not* strip the line: the final \n is required for the reshaping.
        *idx, calls = line.split(b'\t', 6)
        calls = bytes_to_ords(calls, cols=cols, rows=4, stride=2)
        yield idx, calls


def parse_vcf(stream, cohort=None, flare=False):
    """Parse an iterator of newline-terminated bytes as Phased Genotypes"""
    def header_fn(parts):
        return tuple(parts[:5]), np.array(parts[9:])

    for line in stream:
        # Skip the bcftools VCF headers
        if line.startswith(b'##'):
            continue
        # Emit the header with the 5 kept field names and IIDs
        elif line.startswith(b'#'):
            index, iids, cols = process_iid_header(line, header_fn, cohort)
            yield index, iids
            break
    ## This is how to make an array for checking the cell format (not impl.)
    ## XXX: checking would be pricy - double call time - maybe as an assert?
    ## Then could be run with the python "strip asserts" optimization mode.
    # seps = np.array(list(map(ord, b'|::\t')), dtype=np.int8)    # FLARE
    # seps = np.array(list(map(ord, b'|\t')), dtype=np.int8)      # RFMix

    # RFMix: Cell format is GT0|GT1\t (4 bytes). Reshape to 4 rows, then
    #   rows 0 and 2 are GT0 and GT1, while rows 1 and 3 are '|' and
    #   whitespace.
    # FLARE: Cell format is GT0|GT1:ANC0:ANC1\t (8 bytes). Reshape to 8
    #   rows then rows 1, 3, 5, 7 should be identically '|', ':', ':', and
    #   whitespace, while rows 0, 2, 4, and 6 are GT0, GT1, ANC0, and ANC1.
    rows = 8 if flare else 4

    for line in stream:
        # Do *not* strip the line: the final \n is required for the reshaping.
        # Keep the first 5 fields; ditch QUAL, FILTER, INFO, and FORMAT
        *idx, _, _, _, _, calls = line.split(b'\t', 9)
        calls = bytes_to_ords(calls, cols=cols, rows=rows, stride=2)
        yield idx, calls


def parse_anc(stream, cohort=None):
    """Parse an iterator of dosage (or hapcount) byte rows"""
    def header_fn(parts):
        return tuple(parts[:5]), np.array(parts[5:])

    # First line is always the header. NOTE: header not always b'#' prefixed!
    index, iids, cols = process_iid_header(next(stream), header_fn, cohort)
    yield index, iids

    for line in stream:
        *idx, calls = line.split(b'\t', 5)
        calls = bytes_to_ords(calls, cols=cols, rows=2, stride=2)
        yield idx, calls


def unkink_msp(stream):
    """Detect and correct strand flips in an RFMix MSP file.

    The input stream is (index, call data) tuples.
    The output stream is (index, call data, switched flags) tuples.
    """
    idx, iids = next(stream)
    yield idx, iids

    n = len(iids)
    prev_start = np.zeros(n, dtype=float)
    prev_a = np.zeros(n, dtype=np.int8)
    prev_b = np.zeros(n, dtype=np.int8)
    switched = np.zeros(n, dtype=bool)

    for idx, (a, b) in stream:
        start = np.float64(idx[3]) - 1.0
        egpos = np.float64(idx[4])

        update_mask = (a != b)
        change_mask = (a != prev_a) | (b != prev_b)
        starts_mask = (start < prev_start)
        sm = (update_mask & change_mask & starts_mask)
        switched[sm] = np.logical_not(switched[sm])

        # Update state tracking ndarrays
        prev_start[update_mask] = egpos
        prev_a[update_mask] = a[update_mask]
        prev_b[update_mask] = b[update_mask]

        # Perform the switch
        a_out = np.where(switched, b, a)
        b_out = np.where(switched, a, b)
        yield idx, (a_out, b_out), switched


def _check_stream_iids(vcf_stream, msp_stream, check_iids=True):
    """Raise an error if the IIDs do not match"""
    vcf_head, vcf_iids = next(vcf_stream)
    msp_head, msp_iids = next(msp_stream)
    if check_iids and not np.all(vcf_iids == msp_iids):
        raise ValueError('VCF and MSP iids do not match')
    return (vcf_head, vcf_iids), (msp_head, msp_iids)


def unkink_vcf(vcf_stream, msp_stream, check_iids=True):
    """
    Use switches and MSP windows from msp_stream to detect and correct strand
    flips in a Phased Genotypes in vcf_stream (a parsed Phased VCF file).

    Arguments
    ---------
    vcf_stream: iterator of (idx, (geno_a, geno_b)) tuples
    msp_stream: iterator of (idx, (call_a, call_b), switched) tuples
    check_iids: bool

    Yields first (VCF field header tuple, VCF IID ndarray)
    Yields after that (idx, (geno_a, geno_b)
    """
    yield _check_stream_iids(vcf_stream, msp_stream, check_iids=check_iids)[0]

    windows = ((idx[0], int(idx[1]), int(idx[2]), switched)
               for idx, _, switched in iter(msp_stream))
    chrom, start, stop, switched = next(windows)

    for vcf_idx, (a, b) in vcf_stream:
        pos = int(vcf_idx[1])
        while not (vcf_idx[0] == chrom and (start <= pos < stop)):
            try:
                chrom, start, stop, switched = next(windows)
            except StopIteration:
                break
        a_out = np.where(switched, b, a)
        b_out = np.where(switched, a, b)
        yield vcf_idx, (a_out, b_out)


def _hap_and_dos(geno_a, geno_b, call_a, call_b, call_grid):
    # Get the haplotype counts
    call_grid_a = (call_a == call_grid)
    call_grid_b = (call_b == call_grid)
    hap_cts = call_grid_a.view(np.int8) + call_grid_b.view(np.int8)

    # Get the dosage counts
    geno_grid_a = (geno_a.view(bool) & call_grid_a)
    geno_grid_b = (geno_b.view(bool) & call_grid_b)
    dos_cts = geno_grid_a.view(np.int8) + geno_grid_b.view(np.int8)
    return hap_cts, dos_cts


def extract_rfmix(vcf_stream, msp_stream, ancestors, check_iids=True):
    """Extract tracts from RFMix VCF and MSP streams"""
    idx, iids = _check_stream_iids(vcf_stream, msp_stream, check_iids)[0]
    ancestors = np.asarray(ancestors).reshape(-1, 1).astype(np.int8)
    call_grid = np.tile(ancestors, len(iids))
    yield idx, iids
    del idx, iids

    windows = ((idx[0], int(idx[1]), int(idx[2]), call_a, call_b)
               for idx, (call_a, call_b), *_ in iter(msp_stream))
    msp_chrom, start, stop, call_a, call_b = next(windows)

    for idx, (geno_a, geno_b) in vcf_stream:
        # First align the VCF and MSP: the VCF site must be in an MSP window
        vcf_chrom, pos = idx[0], int(idx[1])
        if vcf_chrom == msp_chrom and pos < start:
            continue
        skip_line = False
        while not (vcf_chrom == msp_chrom and (start <= pos < stop)):
            try:
                msp_chrom, start, stop, call_a, call_b = next(windows)
            except StopIteration:
                break
            if vcf_chrom == msp_chrom and pos < start:
                skip_line = True
                break
        if skip_line:
            continue

        hap, dos = _hap_and_dos(geno_a, geno_b, call_a, call_b, call_grid)
        yield idx, hap, dos


def extract_flare(stream, ancestors):
    """Extract tracts from FLARE VCF stream"""
    # Emit the header line, and learn the size of the call grid from it
    idx, iids = next(stream)
    ancestors = np.asarray(ancestors).reshape(-1, 1).astype(np.int8)
    call_grid = np.tile(ancestors, len(iids))
    yield idx, iids
    del idx, iids

    # `calls` ought to have four rows (shape [4, n_iids]: see parse_anc)
    # i.e.: >>> geno_a, geno_b, call_a, call_b = calls
    for (idx, calls) in stream:
        hap, dos = _hap_and_dos(*calls, call_grid)
        yield idx, hap, dos


def write_tracts(stream, hap_files, dos_files, header_char=True):
    """Consume tract data from stream and write to all input files in tandem"""
    idx, iids = next(stream)

    # Handle compatibility with respect to the header char '#'... the original
    # does not use '#' to prefix the header for the dosage / hapcount files.
    if not header_char:
        idx = (idx[0].lstrip(b'#'), *idx[1:])
    elif not idx[0].startswith(b'#'):
        idx = (b'#'+idx[0], *idx[1:])

    header = b'\t'.join(itertools.chain(idx, iids))
    for file in itertools.chain(hap_files, dos_files):
        file.write(header)
        file.write(b'\n')

    # Make a line buffer
    sbuff = np.full(2*len(iids), ord(b'\t'), dtype=np.int8)
    sbuff[-1] = ord(b'\n')
    n_anc = len(hap_files)
    for (idx, hap, dos) in stream:
        hap_ords = hap + 48
        dos_ords = dos + 48
        index = b'\t'.join(idx)
        for i in range(n_anc):
            hap_files[i].write(index)
            hap_files[i].write(b'\t')
            sbuff[::2] = hap_ords[i]
            hap_files[i].write(sbuff)

            dos_files[i].write(index)
            dos_files[i].write(b'\t')
            sbuff[::2] = dos_ords[i]
            dos_files[i].write(sbuff)


def hstack_ancs(streams, sort=True, check_overlap=True, check_index=True):
    """
    Given N input streams with identical index fields and non-overlapping
    sample ids, merge the columns of the input streams.
    """
    streams = tuple(streams)
    if len(streams) < 2:
        raise ValueError('Not enough streams to merge')

    fields, iids = zip(*[next(s) for s in streams])
    if not np.all(fields[0] == np.array(fields[1:])):
        raise ValueError('Header index field mismatch')

    iids = np.concatenate(iids)
    if check_overlap and not np.all(np.sort(iids) == np.unique(iids)):
        raise ValueError('Overlapping sample ids')

    if sort:
        reorder = np.argsort(iids)
        iids = iids[reorder]

    yield fields[0], iids
    while True:
        try:
            rows = [next(s) for s in streams]
        except StopIteration:
            break
        else:
            idx, vals = zip(*rows)
            if check_index:
                idx_arr = np.array(idx)
                if not np.all(idx_arr[0] == idx_arr[1:]):
                    raise ValueError('Index mismatch')
            vals = np.concatenate(vals)
            if sort:
                vals = vals[reorder]
            yield idx[0], vals


def write_msp(stream, msp_name, switch_name, opener=builtins.open):
    """Write unkinked MSP stream to file and forward to consumer function."""
    # NOTE: function is assumed to be called with unkink_msp stream as input
    with opener(msp_name, 'wb') as file, opener(switch_name, 'wb') as sw_file:
        # Write the MSP header, ignoring the Subpopulations header
        fields, iids = next(stream)
        if not fields[0].startswith(b'#'):
            file.write(b'#')
        # MSP-format the IID byte strings
        msp_iids = ((ix + b'.0', ix + b'.1') for ix in iids)
        msp_iids = itertools.chain.from_iterable(msp_iids)
        file.write(b'\t'.join(itertools.chain(fields, msp_iids)))
        file.write(b'\n')
        # Finally remember to pass the header data through
        yield fields, iids

        sbuff = np.full(4*len(iids), ord(b'\t'), dtype=np.int8)
        sbuff[-1] = ord(b'\n')
        for idx, (call_a, call_b), switches in stream:
            idx_bts = b'\t'.join(idx)
            # Write to the unkinked MSP file
            sbuff[0::4] = call_a + 48
            sbuff[2::4] = call_b + 48
            file.write(idx_bts)
            file.write(b'\t')
            file.write(sbuff)
            # Write to the switches file
            sw_file.write(idx_bts)
            sw_file.write(b'\t')
            sw_file.write(b'\t'.join(switches.astype('S5')))
            sw_file.write(b'\n')
            yield idx, (call_a, call_b), switches


def write_vcf(stream, filename, opener=builtins.open):
    """Write unkinked VCF stream to file and forward to consumer function"""
    with opener(filename, 'wb') as file:
        fields, iids = next(stream)
        if not fields[0].startswith(b'#'):
            file.write(b'#')
        file.write(b'\t'.join(itertools.chain(fields, iids)))
        file.write(b'\n')
        yield fields, iids

        sbuff = np.full(4*len(iids), ord(b'\t'), dtype=np.int8)
        sbuff[1::4] = ord(b'|')
        sbuff[-1] = ord(b'\n')
        for idx, (geno_a, geno_b) in stream:
            sbuff[0::4] = geno_a + 48
            sbuff[2::4] = geno_b + 48
            file.write(b'\t'.join(idx))
            # Write The next four header fields: QUAL, FILTER, INFO, FORMAT.
            # These are stripped out when we parse the file (they are not
            # involved in the ANC hapcount / dosage outputs).
            file.write(b'\t.\t.\t.\tGT\t')
            file.write(sbuff)
            yield idx, (geno_a, geno_b)


def write_anc(stream, filename, opener=builtins.open):
    """
    Write a single dosage / hapcount stream to a file. Used to write the
    output of, e.g., hstack_ancs, and *not* the output(s) of extract_ancestors.
    """
    with opener(filename, 'wb') as file:
        fields, iids = next(stream)
        file.write(b'\t'.join(itertools.chain(fields, iids)))
        file.write(b'\n')

        sbuff = np.full(2*len(iids), ord(b'\t'), dtype=np.int8)
        sbuff[-1] = ord(b'\n')
        for idx, vals in stream:
            file.write(b'\t'.join(idx))
            file.write(b'\t')
            sbuff[::2] = vals + 48
            file.write(sbuff)


def cli():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    parser_rfmix = subparsers.add_parser('extract_rfmix')
    parser_flare = subparsers.add_parser('extract_flare')
    parser_merge = subparsers.add_parser('merge_tracts')
    all_parsers = (parser, parser_rfmix, parser_flare, parser_merge)

    # RFMix specific arguments
    f = parser_rfmix.add_argument
    f('vcf', type=str, help='Path to Phased VCF Genotypes file')
    f('--msp', type=Path, help='Path to RFMix MSP file')

    f('--msp-gzip', action=BooleanOptionalAction, default=False,
      help='(default False) Treat MSP input as gzip compressed')

    f('--unkink', action=BooleanOptionalAction, default=True,
      help='(default True) Correct RFMix strand flips in MSP and VCF files')

    f('--emit-msp', action=BooleanOptionalAction, default=True, help="""
        (default True) Write unkinked MSP file and MSP Switches file. Files are
        written to `--outdir` as `{stem}_unkinked.msp[.gz]` and
        `{stem}_switches.msp[.gz]` (see also `--stem` and `--out-gzip`)
    """)

    f('--emit-vcf', action=BooleanOptionalAction, default=True, help="""
        (default True) Write unkinked VCF file. File is written to `--outdir`
        as `{stem}_unkinked.vcf[.gz]` (see also `--stem` and `--out-gzip`
    """)

    # FLARE specific arguments
    f = parser_flare.add_argument
    f('vcf', type=str, help='Path to FLARE formatted VCF file')

    # Dosage / hapcount merging arguments
    f = parser_merge.add_argument
    f('files', type=Path, nargs='+', help='Input files to merge')
    f('-o', '--output', type=Path, help='Output file for merge')

    f('--sort', action=BooleanOptionalAction, default=True,
      help='(default True) Sort merged output by sample IIDs')

    f('--check-overlap', action=BooleanOptionalAction, default=True,
      help='(default True) Raise an error if the input file IIDs overlap')

    f('--check-index', action=BooleanOptionalAction, default=True, help="""
        (default True) Raise an error if input file rows mismatch on CHROM,
        'POS, ID, REF, or ALT fields.
    """)

    # Arguments common to both RFMix and FLARE
    for f in (parser_rfmix.add_argument, parser_flare.add_argument):
        f('--ancestors', nargs='*', type=int, help="""
            Specific ancestral population(s) to extract. E.g., `--ancestors 0
            3` would output dosage and hapcount files for only population codes
            `0` and `3`. See `--num-ancestors` for filename format.
        """)

        f('--num-ancestors', type=int, help="""
            Number (0..N) of ancestral population(s) to extract. E.g.,
            `--num-ancestors=6` will extract codes `[0, 1, 2, 3, 4, 5]. Output
            files are written to `--outdir` with filenames formatted as
            `{stem}.anc{i}.hapcount.txt[.gz]` and
            `{stem}.anc{i}.dosage.txt[.gz]`. See `--stem` and `--out-gzip`.
        """)

        f('--stem', type=str, help="""
            Output filename stem. If not provided, the input VCF filename is
            used, stripped of all suffixes. I.e., the stem taken from
            `phased_chr1.vcf.gz` is `phased_chr`.
        """)

    # CLI options common to all commands must be added to all parsers
    for f in (p.add_argument for p in all_parsers):
        f('--dry-run', type=BooleanOptionalAction,
          help='Print parsed command line arguments and exit')

        f('--outdir', type=Path, required=True,
          help='Directory in which to write all output files')

        f('--vcf-gzip', action=BooleanOptionalAction, default=True,
          help='(default True) Treat input VCF as gzip compressed')

        f('--out-gzip', action=BooleanOptionalAction, default=True, help="""
            (default True) Gzip compress all output files. Uses highest
            compression level. Uses ISAL if installed.
        """)

        f('--gzip-threads', type=int, default=0, help="""
            Use multiple threads for gzip compression. Requires ISAL. Default:
            single-threaded. NOTE: decompression does not benefit from extra
            threads, and thread communication overhead may, depending on input,
            make multiple threads slower than the default.
        """)

        f('--header-char', action=BooleanOptionalAction, default=False,
          help="""
            (Default False) Prefix output dosage and hapcount file header lines
            with a '#' character.
        """)

        f('--iids', type=Path, help="""
            Optional IID file to filter inputs. The file must have one IID per
            line and no header line. IIDs in this file missing from the inputs
            are silently ignored. If provided, the IIDs are leically sorted,
            and the output columns are therefore also sorted.
        """)

    args = parser.parse_args()
    if args.dry_run:
        pprint(vars(args))
        return

    #### Process all the common args first
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Define a function for opening gzip files at the right level & threading
    def gz_opener(name, mode, threads=args.gzip_threads):
        if HAS_ISAL:
            return gzip_open(name, mode, compresslevel=3, threads=threads)
        else:
            return gzip_open(name, mode, compresslevel=9)

    # Determine if the inputs / outputs ought to use the gz_opener function
    if args.vcf_gzip:
        vcf_open = gz_opener
    else:
        vcf_open = builtins.open

    if args.out_gzip:
        suffix = '.gz'
        out_open = gz_opener
    else:
        suffix = ''
        out_open = builtins.open

    # If we must filter the input files to a given set of IIDs, read the
    # cohort. NB: read_bytes() is pathilb.Path function (fun!)
    if args.iids:
        cohort = np.sort(args.iids.read_bytes().strip().split(b'\n'))
    else:
        cohort = None

    #### Run the actual command
    if args.command in ('extract_rfmix', 'extract_flare'):
        # Process args common to both RFMix and Flare extraction first
        # A specific list of ancestors is preferred to a number range
        if args.ancestors:
            ancestors = np.array(args.ancestors)
        elif args.num_ancestors:
            ancestors = np.arange(args.num_ancestors)
        else:
            raise RuntimeError('--ancestors or --num-ancestors is required')

        # Use user stem if given; otherwise strip *all* suffices to get stem
        if args.stem:
            stem = args.stem
        else:
            if args.vcf == '-' or args.vcf == 'stdin':
                raise RuntimeError('Cannot determine output file name stem')
            p0 = Path(args.vcf)
            while p0 != p0.with_suffix(''):
                p0 = p0.with_suffix('')
            stem = p0.stem

        # file names for the extracted hapcount / dosage files
        hap_names = [f'{stem}.anc{i}.hapcount.txt{suffix}' for i in ancestors]
        dos_names = [f'{stem}.anc{i}.dosage.txt{suffix}' for i in ancestors]

        #### RFMix extraction stream: command is extract_rfmix
        if args.command == 'extract_rfmix':
            if args.msp_gzip:
                msp_open = gz_opener
            else:
                msp_open = builtins.open

            vcf_stream = parse_vcf(read_file(args.vcf, opener=vcf_open), cohort)
            msp_stream = parse_msp(read_file(args.msp, opener=msp_open), cohort)

            if args.unkink:
                msp_stream = unkink_msp(msp_stream)
                if args.emit_msp:
                    msp_stream = write_msp(
                        msp_stream,
                        args.outdir/f'{stem}_unkinked.msp{suffix}',
                        args.outdir/f'{stem}_switches.msp{suffix}',
                        opener=out_open
                    )
                # The MSP file is iterated *twice* - once when unkinking and
                # again when doing the extraction. So here we consume the
                # interator so that we have an *iterable* to work with.
                msp_data = list(msp_stream)

                vcf_stream = unkink_vcf(vcf_stream, iter(msp_data))
                if args.emit_vcf:
                    path = args.outdir/f'{stem}_unkinked.vcf{suffix}'
                    vcf_stream = write_vcf(vcf_stream, path, opener=out_open)

            stream = extract_rfmix(vcf_stream, iter(msp_data), ancestors)

        #### FLARE extraction stream: command is extract_flare
        else:
            stream = read_file(args.vcf, opener=vcf_open)
            stream = parse_vcf(stream, cohort=cohort, flare=True)
            stream = extract_flare(stream, ancestors)

        #### Run the extraction (consume the `stream` determined just above)
        try:
            hap_files = [out_open(args.outdir/name, 'wb') for name in hap_names]
            dos_files = [out_open(args.outdir/name, 'wb') for name in dos_names]
            write_tracts(stream, hap_files, dos_files, args.header_char)
        finally:
            for fh in itertools.chain(hap_files, dos_files):
                fh.close

    elif args.command == 'merge_tracts':
        merged = hstack_ancs(
            streams=[
                parse_anc(read_file(f, opener=vcf_open), cohort=cohort)
                for f in args.files
            ],
            sort=args.sort,
            check_overlap=args.check_overlap,
            check_index=args.check_index,
        )
        write_anc(merged, filename=args.outdir/args.output, opener=out_open)

    else:
        raise RuntimeError('Unknown command')


if __name__ == '__main__':
    cli()
