"""
Convert RFMix outputs (VCF and MSP files) to Tractor inputs.
"""
import argparse
import builtins
import contextlib
import copy
import io
import itertools
import logging
import numbers
import operator
import pathlib
import shutil
from argparse import ArgumentParser, BooleanOptionalAction
from pathlib import Path

import numpy as np
try:
    from isal import igzip as gzip
except ImportError:
    import gzip


__all__ = [
    'read_iids',
    'Reader',
    'Writer',
    'BaseVCF',
    'MSP',
    'MSPSwitches',
    'PhasedVCF',
    'AncVCF',
    'MultiAncVCF',
    'unkink_msp',
    'unkink_geno',
    'extract_tracts',
]


def read_iids(fname, header=False, sort=True):
    """Convenience function for loading IID from file (one per line)"""
    with open(fname, 'r') as file:
        iids = [x.strip() for x in file]
    if header:
        iids = iids[1:]
    if sort:
        iids = sorted(iids)
    return np.array(iids)


def zip_on_index(keys, its, index_fn):
    """
    Transforms a dictionary `its` of keys to iterators to a stream of values
    from each iterator keyed by the provided `keys`, with a common index
    provided by `index_fn`. Yields 2-tuples (index, values_dict).
    """
    for row_i in itertools.count(0):
        row_index = None
        row_value = {}
        for k in keys:
            try:
                idx, row = index_fn(next(its[k]))
            except StopIteration:
                return
            else:
                if row_index is None:
                    row_index = idx
                elif row_index != idx:
                    raise RuntimeError(f'Row Index {row_i} Mismatch')
                row_value[k] = row
        yield row_index, row_value


def msp_fmt_iids(iids):
    for ix in iids:
        yield f'{ix}.0'
        yield f'{ix}.1'


class Reader:
    def __init__(self, path, compress=False, buffer=False):
        self.path = Path(path).resolve()
        self.compress = compress
        self.buffer = buffer

    def open(self):
        match (self.buffer, self.compress):
            case (True, True):
                self._buff = io.BytesIO()
                with builtins.open(self.path, 'rb') as src:
                    shutil.copyfileobj(src, self._buff)
                self._buff.seek(0)
                self._file = gzip.open(self._buff, 'rt')
            case (True, False):
                self._file = io.StringIO()
                with builtins.open(self.path, 'rt') as src:
                    shutil.copyfileobj(src, self._file)
                self._file.seek(0)
            case (False, True):
                self._file = gzip.open(self.path, 'rt')
            case (False, False):
                self._file = builtins.open(self.path, 'r')

    def close(self):
        with contextlib.suppress(AttributeError):
            self._file.close()
            self._buff.close()

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, *exc):
        self.close()

    def __iter__(self):
        yield from self._file


class Writer:
    def __init__(self, path, compress=False, buffer=False):
        self.path = path
        self.compress = compress
        self.buffer = buffer

    def open(self):
        match (self.buffer, self.compress):
            case (True, True):
                self._buff = io.BytesIO()
                self._file = gzip.open(self._buff, 'wt')
            case (True, False):
                self._buff = io.BytesIO()
                self._file = io.TextIOWrapper(self._buff)
            case (False, True):
                self._file = gzip.open(self.path, 'wt')
            case (False, False):
                self._file = builtins.open(self.path, 'wt')
        self.write = self._file.write

    def close(self):
        def _write_buff():
            self._buff.seek(0)
            with builtins.open(self.path, 'wb') as disk:
                shutil.copyfileobj(self._buff, disk)
            self._buff.close()

        # Very annoyingly, the order of closing file vs buff is different
        match (self.buffer, self.compress):
            case (True, True):
                self._file.flush()
                self._file.close()
                _write_buff()
            case (True, False):
                self._file.flush()
                _write_buff()
                self._file.close()
            case _:
                self._file.close()

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, *exc):
        self.close()


def parse_reader(reader, char='#', split='\t', n=-1):
    """Extracts headers and rows"""
    reader = iter(reader)
    if isinstance(char, numbers.Integral):
        headers = [next(reader) for _ in range(char)]
        peek = next(reader)
    else:
        headers = []
        while (peek := next(reader)).startswith(char):
            headers.append(peek)
    headers = [h.strip() for h in headers]
    rows = itertools.chain([peek], reader)
    rows = (r.strip().split(split, n) for r in rows)
    return headers, rows


def make_iid_header(obj, iids):
    fields = [f'#{obj.index_fields[0]}']
    fields.extend(obj.index_fields[1:])
    fields.extend(tuple(iids))
    return '\t'.join(fields)


def write_headers(write, headers):
    for h in headers:
        write(h.strip())
        write('\n')

def write_idx_rows(write, idx, rows):
    for ix, rx in zip(idx, rows):
        write('\t'.join(ix))
        write('\t')
        write('\t'.join(rx))
        write('\n')


class BaseVCF:
    fields = ('headers', 'index', 'X0', 'X1')

    def __init__(self, headers, index, X0, X1):
        self.headers = np.asarray(headers)
        self.index = np.asarray(index)
        self.X0 = np.asarray(X0)
        self.X1 = np.asarray(X1)

    @classmethod
    def from_npz(cls, path):
        data = np.load(path, allow_pickle=True)
        return cls(**{k: data[k] for k in cls.fields})

    def to_npz(self, path):
        output = {k: getattr(self, k) for k in self.fields}
        np.savez_compressed(path, **output)

    def __eq__(self, other):
        for field in self.fields:
            if not np.array_equal(getattr(self, field), getattr(other, field)):
                return False
        return True


class MSP(BaseVCF):
    index_fields = ('chm', 'spos', 'epos', 'sgpos', 'egpos', 'n_snps')
    index_dtypes = (str, int, int, float, float, int)

    @property
    def iids(self):
        if not hasattr(self, '_iids'):
            raw = self.headers[-1].strip().split('\t')[6::2]
            raw = [x.split('.')[0] for x in raw]
            self._iids = np.array(raw)
        return self._iids

    @classmethod
    def from_tsv(cls, path, compress=False, buffer=False):
        index, X0, X1 = [], [], []
        with Reader(path, compress, buffer) as reader:
            headers, rows = parse_reader(reader)
            for r in rows:
                index.append(r[:6])
                calls = r[6:]
                X0.append(np.array(calls[0::2], dtype=np.int8))
                X1.append(np.array(calls[1::2], dtype=np.int8))
        return cls(headers=headers, index=index, X0=X0, X1=X1)

    def to_tsv(self, path, compress=False, buffer=False):
        m, n = self.X0.shape
        X = np.empty((m, 2*n), dtype=str)
        X[:, 0::2] = self.X0
        X[:, 1::2] = self.X1
        with Writer(path, compress, buffer) as file:
            write_headers(file.write, self.headers)
            write_idx_rows(file.write, self.index, X)

    def subset(self, iids, inplace=False):
        """Produce a subset of this MSP filtered to specific iids"""
        _, _, locs = np.intersect1d(iids, self.iids, return_indices=True)
        *headers, _ = self.headers.tolist()
        headers.append(make_iid_header(self, msp_fmt_iids(self.iids[locs])))
        if inplace:
            self.headers = headers
            self.X0 = np.copy(self.X0[:, locs])
            self.X1 = np.copy(self.X1[:, locs])
            return self
        else:
            return self.__class__(headers=headers,
                                  index=np.copy(self.index),
                                  X0=np.copy(self.X0[:, locs]),
                                  X1=np.copy(self.X1[:, locs]))

    def sort_by_iids(self, inplace=False):
        return self.subset(self.iids, inplace=inplace)

    def merge(self, other):
        """Produce a single MSP by merging self with other"""
        A_iids = np.copy(self.iids)
        B_iids = np.copy(other.iids)
        iids = np.concatenate((A_iids, B_iids))
        sort = np.argsort(iids)
        iids = iids[sort]

        *headers, _ = self.headers.tolist()
        headers.append(make_iid_header(self, msp_fmt_iids(iids)))

        # XXX: verify the correctness of this algorithm!
        index, X0, X1 = [], [], []
        A_zip = zip(self.index, self.X0, self.X1)
        B_zip = zip(other.index, other.X0, other.X1)
        A_idx, A_x0, A_x1 = next(A_zip, (None, None, None))
        B_idx, B_x0, B_x1 = next(B_zip, (None, None, None))
        while (A_idx is not None) and (B_idx is not None):
            if A_idx[2] == B_idx[2]:
                index.append(np.copy(A_idx))
                X0.append(np.concatenate((A_x0, B_x0))[sort])
                X1.append(np.concatenate((A_x1, B_x1))[sort])
                A_idx, A_x0, A_x1 = next(A_zip, (None, None, None))
                B_idx, B_x0, B_x1 = next(B_zip, (None, None, None))
            elif int(A_idx[2]) > int(B_idx[2]):
                index.append(np.copy(A_idx))
                index[-1][2] = B_idx[2]
                index[-1][4] = B_idx[4]
                X0.append(np.concatenate((A_x0, B_x0))[sort])
                X1.append(np.concatenate((A_x1, B_x1))[sort])
                B_idx, B_x0, B_x1 = next(B_zip, (None, None, None))
            else:
                index.append(np.copy(A_idx))
                X0.append(np.concatenate((A_x0, B_x0))[sort])
                X1.append(np.concatenate((A_x1, B_x1))[sort])
                A_idx, A_x0, A_x1 = next(A_zip, (None, None, None))
        X0 = np.copy(np.asarray(X0), order='C')
        X1 = np.copy(np.asarray(X1), order='C')
        return self.__class__(headers=headers, index=index, X0=X0, X1=X1)


class MSPSwitches(BaseVCF):
    # The original switches files don't have any header at all - but there's no
    # harm in keeping the header with the sample IID information if available.
    fields = ('headers', 'index', 'switches')
    index_fields = tuple(MSP.index_fields)
    index_dtypes = tuple(MSP.index_dtypes)

    def __init__(self, index, switches, headers=None):
        if headers is None:
            headers = []
        self.headers = np.asarray(headers)
        self.index = np.asarray(index)
        self.switches = np.asarray(switches)

    @classmethod
    def from_tsv(cls, path, compress=False, buffer=False):
        index, switches = [], []
        with Reader(path, compress, buffer) as reader:
            headers, rows = parse_reader(reader)
            for r in rows:
                index.append(r[:6])
                switches.append(np.array(r[6:]) == 'True')
        return cls(headers=headers, index=index, switches=switches)

    def to_tsv(self, path, compress=False, buffer=False):
        with Writer(path, compress, buffer) as file:
            write_headers(file.write, self.headers)
            write_idx_rows(file.write, self.index, self.switches.astype(str))


class PhasedVCF(BaseVCF):
    """Wraps conversion and manipulation functions for an Phased VCF file"""
    index_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT')
    index_dtypes = (str, int, str, str, str, str, str, str, str)

    @property
    def iids(self):
        if not hasattr(self, '_iids'):
            self._iids = np.array(self.headers[-1].strip().split('\t')[9:])
        return self._iids

    @classmethod
    def from_tsv(cls, path, compress=False, buffer=False):
        index, X0, X1 = [], [], []
        with Reader(path, compress, buffer) as reader:
            headers, rows = parse_reader(reader, n=9)
            for *idx, calls in rows:
                index.append(tuple(idx))
                calls = np.asarray([calls]).view('U1')
                X0.append(calls[0::4] == '1')
                X1.append(calls[2::4] == '1')
        return cls(headers=headers, index=index, X0=X0, X1=X1)

    def to_tsv(self, path, compress=False, buffer=False):
        X = np.full(4 * len(self.X0[0]), '0')
        X[1::4] = '|'
        X[3::4] = '\t'
        X[-1] = '\n'
        vtype = f'U{len(X)}'
        with Writer(path, compress, buffer) as file:
            write_headers(file.write, self.headers)
            for idx, X0, X1 in zip(self.index, self.X0, self.X1):
                X[0::4] = np.where(X0, '1', '0')
                X[2::4] = np.where(X1, '1', '0')
                file.write('\t'.join(idx))
                file.write('\t')
                file.write(X.view(vtype)[0])

    @classmethod
    def from_tsv_set(cls, paths, iid_path=None, compress=False, buffer=False):
        """
        Filters, merges, and sorts by sample id a set of VCF files with
        identical index elements and non-overlapping sample ids.
        """
        row_iters = {}
        file_iids = {}
        all_iids = set()

        with contextlib.ExitStack() as stack:
            for p in paths:
                s = stack.enter_context(Reader(p, compress, buffer))
                headers, row_iters[p] = parse_reader(s, n=9)
                file_iids[p] = np.array(headers[-1].strip().split('\t')[9:])
                all_iids.update(file_iids[p])

            if iid_path is not None:
                ref_iids = read_iids(iid_path, sort=True)
            else:
                ref_iids = np.array(sorted(all_iids))

            indices = {}
            for p, iids in file_iids.items():
                _, *idx = np.intersect1d(ref_iids, iids, return_indices=True)
                indices[p] = idx

            def index_fn(r):
                # Strip the INFO and FORMAT fields (except 'GT')
                # Equivalent to bcftools annotate -x INFO,FORMAT
                return (*r[:7], '.', 'GT'), r[9]

            index, X0, X1 = [], [], []
            for idx, calld in zip_on_index(paths, row_iters, index_fn):
                index.append(tuple(idx))
                x0 = np.zeros(len(ref_iids), dtype=bool)
                x1 = np.zeros(len(ref_iids), dtype=bool)
                for k, calls in calld.items():
                    ref_locs, iid_locs = indices[k]
                    calls = np.asarray([calls]).view('U1')
                    x0[ref_locs] = (calls[0::4] == '1')[iid_locs]
                    x1[ref_locs] = (calls[2::4] == '1')[iid_locs]
                X0.append(x0)
                X1.append(x1)

        headers = [make_iid_header(cls, ref_iids)]
        return cls(headers=headers, index=index, X0=X0, X1=X1)


class AncVCF(BaseVCF):
    """Used for `*.hapcount.tsv` and `*.dosage.tsv` files"""
    fields = ('headers', 'index', 'X')
    index_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT')
    index_dtypes = (str, int, str, str, str)

    def __init__(self, headers, index, X):
        self.headers = np.asarray(headers)
        self.index = np.asarray(index)
        self.X = np.asarray(X)

    @property
    def iids(self):
        if not hasattr(self, '_iids'):
            self._iids = np.array(self.headers[-1].strip().split('\t')[5:])
        return self._iids

    @classmethod
    def from_tsv(cls, path, compress=False, buffer=False, no_header_char=True):
        index, X = [], []
        with Reader(path, compress, buffer) as reader:
            char = 1 if no_header_char else '#'
            headers, rows = parse_reader(reader, char=char)
            for r in rows:
                index.append(r[:5])
                X.append(np.array(r[5:], dtype=np.int8))
        return cls(headers=headers, index=index, X=X)

    def to_tsv(self, path, compress=False, buffer=False, compat=False):
        with Writer(path, compress, buffer) as file:
            # Write headers (if compat, without the header char!)
            if compat:
                headers = [h.lstrip('#') for h in self.headers]
            else:
                headers = self.headers
            write_headers(file.write, headers)

            buff = np.full(2*self.X.shape[1], '\t')
            buff[-1] = '\n'
            vtype = f'U{len(buff)}'
            for idx, row in zip(self.index, self.X):
                file.write('\t'.join(idx))
                file.write('\t')
                buff[0::2] = row
                file.write(buff.view(vtype)[0])


class MultiAncVCF(BaseVCF):
    """Stores Ancestor VCF hapcount & dosage files with a common index"""
    fields = ('headers', 'index', 'ancestors', 'hapcounts', 'dosages')
    index_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT')
    index_dtypes = (str, int, str, str, str)

    def __init__(self, headers, index, ancestors, hapcounts, dosages):
        self.headers = np.asarray(headers)
        self.index = np.asarray(index)
        self.ancestors = np.asarray(ancestors)
        self.hapcounts = np.asarray(hapcounts)
        self.dosages = np.asarray(dosages)
        self._ancidx = {v: i for i, v in enumerate(ancestors)}

    def get_anc(self, anc, which='hapcount'):
        idx = self._ancidx[anc]
        match which:
            case 'hapcount':
                X = self.hapcounts[:, idx]
            case 'dosage':
                X = self.dosages[:, idx]
            case _:
                raise ValueError(f'{which=} must be "hapcount" or "dosage"')
        return AncVCF(headers=self.headers, index=self.index, X=X)


def unkink_msp(msp, inplace=False, compat=True):
    if not inplace:
        msp = copy.deepcopy(msp)
    assert msp.X0.shape == msp.X1.shape
    assert len(msp.X0) == len(msp.index)
    start = msp.index[:, 3].astype(float) - 1.0
    egpos = msp.index[:, 4].astype(float)
    S = np.zeros(msp.X0.shape, dtype=bool)

    m, n = msp.X0.shape
    prev_start = np.zeros(n, dtype=float)
    prev_a = np.zeros(n, dtype=msp.X0.dtype)
    prev_b = np.zeros(n, dtype=msp.X1.dtype)
    switched = np.zeros(n, dtype=bool)

    for i in range(m):
        # copies have to be made so the switcheroo is not befuddled by mutation
        a = np.copy(msp.X0[i])
        b = np.copy(msp.X1[i])
        update_mask = (a != b)
        change_mask = (a != prev_a) | (b != prev_b)
        starts_mask = (start[i] < prev_start)
        sm = (update_mask & change_mask & starts_mask)
        switched[sm] = np.logical_not(switched[sm])
        # Update state tracking ndarrays
        prev_start[update_mask] = egpos[i]
        prev_a[update_mask] = a[update_mask]
        prev_b[update_mask] = b[update_mask]
        # Perform the switch
        S[i] = switched
        msp.X0[i, switched] = b[switched]
        msp.X1[i, switched] = a[switched]

    # For some reason the original script removes the other header lines
    if compat:
        msp.headers = msp.headers[1:]

    switches = MSPSwitches(index=np.copy(msp.index), switches=S)
    return msp, switches


def unkink_geno(vcf, switches, inplace=False):
    if not inplace:
        vcf = copy.deepcopy(vcf)

    def _iter_windows():
        for idx, S in zip(switches.index, switches.switches):
            yield idx[0], int(idx[1]), int(idx[2]), S
    iter_windows = _iter_windows()
    chrom, start, stop, mask = next(iter_windows)

    for i, idx in enumerate(vcf.index):
        pos = int(idx[1])
        while not (idx[0] == chrom and (start <= pos < stop)):
            try:
                chrom, start, stop, mask = next(iter_windows)
            except StopIteration:
                break
        a = np.copy(vcf.X0[i, mask])
        b = np.copy(vcf.X1[i, mask])
        vcf.X0[i, mask] = b
        vcf.X1[i, mask] = a
    return vcf


def extract_tracts(vcf, msp, ancestors, compat=True, compress=True):
    # TODO: choose method for implementing `compress`
    # sanity check the data shapes
    assert vcf.X0.shape == vcf.X1.shape
    assert msp.X0.shape == msp.X1.shape
    assert vcf.X0.shape[1] == msp.X0.shape[1]
    assert len(vcf.index) == vcf.X0.shape[0]
    assert len(msp.index) == msp.X0.shape[0]
    n_iids = vcf.X0.shape[1]

    # Get the output header for the dosage / hapcount files
    header = vcf.headers[-1].strip()
    # For some reason the original strips out the '#' header marker char
    if compat:
        header = header.lstrip('#')
    header = header.split('\t', 9)
    header = np.array(['\t'.join(header[:5] + header[9:])])

    def _iter_windows():
        for idx, x0, x1 in zip(msp.index, msp.X0, msp.X1):
            yield idx[0], int(idx[1]), int(idx[2]), x0, x1
    iter_windows = _iter_windows()
    msp_chrom, start, stop, call_a, call_b = next(iter_windows)

    # XXX: "output ancestry-specific VCF" - should we bother keeping?
    vcf_0 = []
    vcf_1 = []
    vcf_index = []

    # TODO: use gzip with BytesIO to incrementally compress the objects in mem?
    hap = []
    dos = []
    index = []

    # numbers.Integral will catch python int but also np.int<8|16|etc> types
    if isinstance(ancestors, numbers.Integral):
        ancestors = np.arange(ancestors)
    ancestors = np.asarray(ancestors)
    call_grid = np.tile(ancestors.reshape(-1, 1), n_iids)

    for idx, geno_a, geno_b in zip(vcf.index, vcf.X0, vcf.X1):
        # First align the VCF and MSP: the VCF site must be in an MSP window
        vcf_chrom, pos = str(idx[0]), int(idx[1])
        if vcf_chrom == msp_chrom and pos < start:
            continue
        skip_line = False
        while not (vcf_chrom == msp_chrom and (start <= pos < stop)):
            try:
                msp_chrom, start, stop, call_a, call_b = next(iter_windows)
            except StopIteration:
                break
            if vcf_chrom == msp_chrom and pos < start:
                skip_line = True
                break
        if skip_line:
            continue

        # If we're processing the row, save the right index parts
        index.append(tuple(idx[:5]))

        # Get the haplotype counts
        call_grid_a = (call_a == call_grid)
        call_grid_b = (call_b == call_grid)
        hap_cts = call_grid_a.astype(np.int8) + call_grid_b.astype(np.int8)
        hap.append(hap_cts)

        # Get the dosage counts
        geno_grid_a = (geno_a & call_grid_a)
        geno_grid_b = (geno_b & call_grid_b)
        dos_cts = geno_grid_a.astype(np.int8) + geno_grid_b.astype(np.int8)
        dos.append(dos_cts)

#        vcf_0.append(np.where(call_grid_a, geno_a, '.'))
#        vcf_1.append(np.where(call_grid_b, geno_b, '.'))
    return header, index, ancestors, hap, dos


def run_cli():
    parser = ArgumentParser()
    f = parser.add_argument

    f('--vcf', required=True, action='append', type=Path,
      help='Path to phased genotypes VCF file')

    f('--msp', required=True, action='append', type=Path,
      help='Path to ancestry calls MSP file')

    f('--iids', type=Path,
      help='IID file (one IID per line, no header)')

#    # XXX: implement this? The results are different (annoyingly)
#    f('--pre-merge', action=BooleanOptionalAction, default=False,
#      help='Unkink MSP and VCF files *before* merging (if multiple given)')

    f('--unkink-msp', action=BooleanOptionalAction, default=True,
      help='Detect and correct strand flips in the RFMix MSP output')

    f('--unkink-geno', action=BooleanOptionalAction, default=True,
      help='Detect and correct strand flips in the RFMix VCF output')

    f('--extract', action=BooleanOptionalAction, default=True,
      help='Extract tracts (requires --num-ancestors or --ancestors)')

    f('--msp-compress', action=BooleanOptionalAction, default=False,
      help='Assume the input MSP files are gzip compressed')

    f('--vcf-compress', action=BooleanOptionalAction, default=True,
      help='Assume the input VCF files are gzip compressed')

    f('--buffer-input', action=BooleanOptionalAction, default=False,
      help='Read the input files into memory before processing')

    f('--buffer-output', action=BooleanOptionalAction, default=False,
      help='Write the output files to a memory buffer before writing to disk')

    f('--outdir', type=Path, default=Path().resolve(),
      help='Directory in which to write outputs (created if non-existent)')

    f('--ancestors', nargs='*', type=int,
      help='Ancestral population(s) to extract from VCFs')

    f('--num-ancestors', type=int,
      help='Ancestral population(s) to extract from VCFs')

    f('--write-msp-npz', action=BooleanOptionalAction, default=True,
      help='Write Unkinked MSP and MSPSwitches as .npz files')

    f('--write-vcf-npz', action=BooleanOptionalAction, default=True,
      help='Write Unkinked Phased VCF as an .npz file')

    f('--write-multi-anc-npz', action=BooleanOptionalAction, default=True,
      help='Write extracted Ancestral tracts as a single .npz file')

    f('--write-anc-tsv', action=BooleanOptionalAction, default=False,
      help='Write each extracted Ancestral tracts as a .tsv file')

    f('--write-anc-npz', action=BooleanOptionalAction, default=False,
      help='Write each extracted Ancestral tracts as an .npz file')

    f('--compress-ancs', action=BooleanOptionalAction, default=False,
      help='Compress output ancestor *.tsv files (if writing any)')

    f('--log', action='store_true', help='Write logs to stdout')
    f('--log-file', help='Write logs to the given filename')
    f('--dry-run', action='store_true', help='Log arguments and exit')

    # Process the arguments
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.unkink_geno:
        args.unkink_msp = True

    if args.num_ancestors:
        args.ancestors = list(range(args.num_ancestors))

    if args.log or args.log_file:
        if args.log_file:
            filename = args.outdir/args.log_file
        else:
            # Logs to stream instead of file
            filename = None
        logging.basicConfig(filename=filename,
                            level=logging.INFO,
                            format='%(asctime)s %(name)s %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    # Pretty-print the arguments
    nsdict = vars(args)
    fmt = f'>{max(len(k) for k in nsdict)}s'
    for key, val in sorted(nsdict.items()):
        logging.info('%s: %s', format(key, fmt), val)

    if args.dry_run:
        return

    #### Run the CLI App
    if args.iids is not None:
        logging.info(f'Reading IID file: {args.iids}')
        iids = read_iids(args.iids)
    else:
        logging.info(f'Using IIDs in files')
        iids = None

    logging.info('Loading VCF file(s)')
    kws = {'compress': args.vcf_compress, 'buffer': args.buffer_input}
    if len(args.vcf) > 1:
        vcf = PhasedVCF.from_tsv_set(args.vcf, args.iids, **kws)
    else:
        vcf = PhasedVCF.from_tsv(args.vcf[0], **kws)

    logging.info('Loading MSP file(s)')
    kws = {'compress': args.msp_compress, 'buffer': args.buffer_input}
    msp = MSP.from_tsv(args.msp[0], **kws)
    if iids is not None:
        msp = msp.subset(iids)
    for more_msp in args.msp[1:]:
        new = MSP.from_tsv(more_msp, **kws)
        if iids is not None:
            new = new.subset(iids)
        msp = msp.merge(new)

    logging.info('Checking that IIDs match between MSP and VCF')
    if not np.all(msp.iids == vcf.iids):
        raise RuntimeError('IID mismatch between MSP and VCF files')

    if args.unkink_msp:
        logging.info('Unkinking MSP')
        msp, switches = unkink_msp(msp)
        if args.write_msp_npz:
            logging.info('Writing Unkinked MSP and Switches .npz files')
            msp.to_npz(args.outdir/'unkinked_msp.npz')
            switches.to_npz(args.outdir/'unkinked_msp_switches.npz')

    if args.unkink_geno:
        logging.info('Unkinking Genotypes')
        vcf = unkink_geno(vcf, switches)
        if args.write_vcf_npz:
            logging.info('Writing Unkinked Genotypes .npz file')
            vcf.to_npz(args.outdir/'unkinked_vcf.npz')

    if not args.extract:
        logging.info('DONE')
        return

    logging.info('Extracting tracts')
    tract_data = extract_tracts(vcf, msp, args.ancestors, compat=False)
    multi_ancs = MultiAncVCF(*tract_data)

    if args.write_multi_anc_npz:
        logging.info('Persisting tracts npz')
        multi_ancs.to_npz(args.outdir/'all_ancestor_tracts.npz')

    if args.write_anc_npz or args.write_anc_tsv:
        tsv_kws = {'compress': args.compress_ancs,
                   'buffer': args.buffer_output}
        for anc in args.ancestors:
            hap = multi_ancs.get_anc(anc, 'hapcount')
            dos = multi_ancs.get_anc(anc, 'dosage')
            if args.write_anc_npz:
                logging.info(f'Writing anc {anc} hapcount npz')
                hap.to_npz(args.outdir/f'anc{anc}_hapcount.npz')

                logging.info(f'Writing anc {anc} dosage npz')
                dos.to_npz(args.outdir/f'anc{anc}_dosage.npz')

            if args.write_anc_tsv:
                logging.info(f'Writing anc {anc} hapcount tsv')
                hap.to_tsv(args.outdir/f'anc{anc}_hapcount.tsv', **tsv_kws)

                logging.info(f'Writing anc {anc} dosage tsv')
                dos.to_tsv(args.outdir/f'anc{anc}_dosage.tsv', **tsv_kws)
    logging.info('DONE')


if __name__ == '__main__':
    run_cli()
