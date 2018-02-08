#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
import os
import os.path as op
import sys
import json

from .utils import get_reader_cls_from_file, get_writer_cls_from_file, is_fasta, is_fastq

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("Trim LQ sequences on both ends, where average MapQV of LQ sequences are less than min_qv.")
    parser.add_argument("input_fn", help="Input FASTA or FASTQ filename")
    parser.add_argument("output_fn", help="Output FASTA or FASTQ filename")
    parser.add_argument("--windowsize", help="Compute average MapQV in windows of size", default=100, type=int)
    parser.add_argument("--min_qv", help="Minimum average MapQV in windows to separate HQ and LQ sequences", default=20, type=int)
    return parser

def get_hq_start_end(seq, is_hq_f):
    """
    ..doctest:
        get_hq_start_end('')
        0, 0
        get_hq_start_end('GC')
        0, 2
        get_hq_start_end('aGCc')
        1, 3
    """
    start, end = 0, 0
    i = 0
    while i < len(seq):
        #if seq[i].isupper():
        if is_hq_f(seq[i]):
            start = i
            break
        i += 1

    i = len(seq) - 1
    while i >= 0:
        #if seq[i].isupper():
        if is_hq_f(seq[i]):
            end = i + 1
            break
        i -= 1
    return (start, end)

def get_hq_start_end_fasta(seq):
    def f(c):
        return c.isupper()
    return get_hq_start_end(seq, f)

def get_hq_start_end_fastq(qual, min_qv):
    def f(c):
        return c >= min_qv
    return get_hq_start_end(qual, f)

def trim_fasta(i_fn, o_fn):
    with get_reader_cls_from_file(i_fn)(i_fn) as reader, get_writer_cls_from_file(o_fn)(o_fn) as writer:
        for r in reader:
            hq_start, hq_end = get_hq_start_end_fasta(r.sequence)
            writer.writeRecord(r.name+'/%s_%s' % (hq_start, hq_end), r.sequence[hq_start:hq_end])

def get_avg_value_in_window(values, windowsize):
    """Given a list of integers, return average value in windows, where ret[i] is the average value of
    values[i, min(i+windowsize, len(values))]
    """
    ret = [0] * len(values)
    for i in range(0, len(values)):
        end = min(i+windowsize, len(values))
        ret[i] = sum(values[i:end]) * 1.0 / (end - i)
    return ret

def get_hq_start_end_above_qv(seq, qual, min_qv, windowsize):
    """
    """
    start, end = 0, 0
    avg_qual = get_avg_value_in_window(qual, windowsize)
    return (start, end)

def trim_fastq(i_fn, o_fn, min_qv, windowsize):
    with get_reader_cls_from_file(i_fn)(i_fn) as reader, get_writer_cls_from_file(o_fn)(o_fn) as writer:
        for r in reader:
            avg_qual = get_avg_value_in_window(r.quality, windowsize)
            hq_start, hq_end = get_hq_start_end_fastq(avg_qual, min_qv)
            writer.writeRecord(r.id+'/%s_%s' % (hq_start, hq_end), r.sequence[hq_start:hq_end], r.quality[hq_start:hq_end])


def trim_lq(i_fn, o_fn, min_qv, windowsize):
    if all(is_fasta(fn) for fn in [i_fn, o_fn]):
        trim_fasta(i_fn, o_fn)
    elif all(is_fastq(fn) for fn in [i_fn, o_fn]):
        trim_fastq(i_fn, o_fn, min_qv, windowsize)
    else:
        raise ValueError("Input and output must be both BED or VCF")


def run(args):
    trim_lq(args.input_fn, args.output_fn, args.min_qv, args.windowsize)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
