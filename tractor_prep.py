import argparse
import builtins
import collections
import contextlib
import copy
import io
import multiprocessing
import os
import pathlib

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
    'unkink_msp',
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
        self.headers = headers
        self.index = index
        self.X0 = X0
        self.X1 = X1

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
        index = np.stack(index)
        X0 = np.stack(X0)
        X1 = np.stack(X1)
        return cls(headers=headers, index=index, X0=X0, X1=X1)

    def to_tsv(self, path):
        m, n = self.X0.shape
        X = np.empty((m, 2*n), dtype=str)
        X[:, 0::2] = self.X0
        X[:, 1::2] = self.X1
        with open(path, 'w') as file:
            for line in self.headers:
                file.write(line)
                file.write('\n')
            for idx, row in zip(self.index, X):
                file.write('\t'.join(idx))
                file.write('\t')
                file.write('\t'.join(row))
                file.write('\n')


class MSPSwitches(BaseVCF):
    fields = ('index', 'switches')
    index_fields = tuple(MSP.index_fields)
    index_dtypes = tuple(MSP.index_dtypes)

    def __init__(self, index, switches):
        self.index = index
        self.switches = switches

    @classmethod
    def from_tsv(cls, path):
        index, switches = [], []
        with open(path, 'r') as file:
            headers, rows = parse_file(file)
            for r in rows:
                index.append(r[:6])
                switches.append(np.array(r[6:]) == 'True')
        return cls(index=np.stack(index), switches=np.stack(switches))

    def to_tsv(self, path):
        with open(path, 'w') as file:
            for idx, row in zip(self.index, self.switches.astype(str)):
                file.write('\t'.join(idx))
                file.write('\t')
                file.write('\t'.join(row))
                file.write('\n')


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
        index = np.stack(index)
        X0 = np.stack(X0)
        X1 = np.stack(X1)
        return cls(headers=headers, index=index, X0=X0, X1=X1)

    def to_tsv(self, path):
        X = np.full(4 * len(self.X0[0]), '0')
        X[1::4] = '|'
        X[3::4] = '\t'
        X[-1] = '\n'
        vtype = f'U{len(X)}'
        with open(path, 'w') as file:
            for line in self.headers:
                file.write(line)
                file.write('\n')
            for idx, X0, X1 in zip(self.index, self.X0, self.X1):
                file.write('\t'.join(idx))
                file.write('\t')
                X[0::4] = np.where(X0, '1', '0')
                X[2::4] = np.where(X1, '1', '0')
                file.write(X.view(vtype)[0])


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


def extract_tracts():
    pass
