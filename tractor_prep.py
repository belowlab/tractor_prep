"""
Convert RFMix outputs (VCF and MSP files) to Tractor inputs.
"""
import copy
import numbers

import numpy as np
try:
    from isal import igzip as gzip
except ImportError:
    import gzip


__all__ = [
    'get_file_headers',
    'parse_file',
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


def get_file_headers(file, char='#'):
    """Extract the headers and move filestream to immediately after"""
    tell = file.seek(0)
    headers = []
    while (line := file.readline()).startswith(char):
        headers.append(line.strip())
        tell = file.tell()
    file.seek(tell)
    return headers


def parse_file(file, char='#', split='\t', n=-1):
    """Extract file headers and split remaining rows on tabs"""
    headers = get_file_headers(file, char=char)
    rows = (r.strip().split(split, n) for r in file)
    return headers, rows


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

    def write_headers(self, file):
        for line in self.headers:
            file.write(line)
            file.write('\n')

    @staticmethod
    def write_row(file, idx, row):
        file.write('\t'.join(idx))
        file.write('\t')
        file.write('\t'.join(row))
        file.write('\n')

    @staticmethod
    def write_rows(file, row_iter):
        for idx, row in row_iter:
            BaseVCF.write_row(file, idx, row)


class MSP(BaseVCF):
    index_fields = ('chm', 'spos', 'epos', 'sgpos', 'egpos', 'n_snps')
    index_dtypes = (str, int, int, float, float, int)

    @classmethod
    def from_tsv(cls, path):
        index, X0, X1 = [], [], []
        with open(path, 'r') as file:
            headers, rows = parse_file(file)
            for r in rows:
                index.append(r[:6])
                calls = r[6:]
                X0.append(np.array(calls[0::2], dtype=np.int8))
                X1.append(np.array(calls[1::2], dtype=np.int8))
        return cls(headers=headers, index=index, X0=X0, X1=X1)

    def to_tsv(self, path):
        m, n = self.X0.shape
        X = np.empty((m, 2*n), dtype=str)
        X[:, 0::2] = self.X0
        X[:, 1::2] = self.X1
        with open(path, 'w') as file:
            self.write_headers(file)
            self.write_rows(file, zip(self.index, X))


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
    def from_tsv(cls, path):
        index, switches = [], []
        with open(path, 'r') as file:
            headers, rows = parse_file(file)
            for r in rows:
                index.append(r[:6])
                switches.append(np.array(r[6:]) == 'True')
        return cls(headers=headers, index=index, switches=switches)

    def to_tsv(self, path):
        with open(path, 'w') as file:
            self.write_rows(file, zip(self.index, self.switches.astype(str)))


class PhasedVCF(BaseVCF):
    """Wraps conversion and manipulation functions for an Phased VCF file"""
    index_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT')
    index_dtypes = (str, int, str, str, str, str, str, str, str)

    @classmethod
    def from_tsv(cls, path):
        index, X0, X1 = [], [], []
        with open(path, 'r') as file:
            headers, rows = parse_file(file, n=9)
            for *idx, calls in rows:
                index.append(tuple(idx))
                calls = np.asarray([calls]).view('U1')
                X0.append(calls[0::4] == '1')
                X1.append(calls[2::4] == '1')
        return cls(headers=headers, index=index, X0=X0, X1=X1)

    def to_tsv(self, path):
        X = np.full(4 * len(self.X0[0]), '0')
        X[1::4] = '|'
        X[3::4] = '\t'
        X[-1] = '\n'
        vtype = f'U{len(X)}'
        with open(path, 'w') as file:
            self.write_headers(file)
            for idx, X0, X1 in zip(self.index, self.X0, self.X1):
                X[0::4] = np.where(X0, '1', '0')
                X[2::4] = np.where(X1, '1', '0')
                file.write('\t'.join(idx))
                file.write('\t')
                file.write(X.view(vtype)[0])


class AncVCF(BaseVCF):
    """Used for `*.hapcount.tsv` and `*.dosage.tsv` files"""
    fields = ('headers', 'index', 'X')
    index_fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT')
    index_dtypes = (str, int, str, str, str)

    def __init__(self, headers, index, X):
        self.headers = np.asarray(headers)
        self.index = np.asarray(index)
        self.X = np.asarray(X)

    @classmethod
    def from_tsv(cls, path, no_header_char=True):
        index, X = [], []
        with open(path, 'r') as file:
            if no_header_char:
                headers = [file.readline().strip()]
                rows = (r.strip().split('\t') for r in file)
            else:
                headers, rows = parse_file(file)
            for r in rows:
                index.append(r[:5])
                X.append(np.array(r[5:], dtype=np.int8))
        return cls(headers=headers, index=index, X=X)

    def to_tsv(self, path, compat=True):
        _headers = self.headers
        with open(path, 'w') as file:
            # Write headers (if compat, without the header char!)
            if compat:
                self.headers = [h.lstrip('#') for h in self.headers]
            self.write_headers(file)

            buff = np.full(2*self.X.shape[1], '\t')
            buff[-1] = '\n'
            vtype = f'U{len(buff)}'
            for idx, row in zip(self.index, self.X):
                file.write('\t'.join(idx))
                file.write('\t')
                buff[0::2] = row
                file.write(buff.view(vtype)[0])
        self.headers = _headers


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
