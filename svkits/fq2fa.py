#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
import os
import os.path as op
import sys
import json

from .utils import is_fastq, is_fasta
from pbcore.io import FastaReader, FastaWriter, FastqReader, FastqWriter

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("Trim LQ sequences on both ends, where average MapQV of LQ sequences are less than min_qv.")
    parser.add_argument("input_fq_fn", help="Input FASTA or FASTQ filename")
    parser.add_argument("output_fa_fn", help="Output FASTA or FASTQ filename")
    parser.add_argument("--windowsize", help="Compute average MapQV in windows of size", default=100, type=int)
    parser.add_argument("--min_qv", help="Minimum average MapQV in windows to separate HQ and LQ sequences", default=20, type=int)
    return parser

def get_reader_cls_from_file(fn):
    if is_fasta(fn):
        return FastaReader
    elif is_fastq(fn):
        return FastqReader
    else:
        raise ValueError("Could not get reader for %s" % fn)

def get_writer_cls_from_file(fn):
    if is_fasta(fn):
        return FastaWriter
    elif is_fastq(fn):
        return FastqWriter
    else:
        raise ValueError("Could not get reader for %s" % fn)

def fq2fa(i_fq_fn, o_fa_fn):
    if not is_fastq(i_fq_fn) or not is_fasta(o_fa_fn):
        raise ValueError("Input %s must be FASTA and output %s must be FASTQ" % (i_fq_fn, o_fa_fn))
    with get_reader_cls_from_file(i_fq_fn)(i_fq_fn) as reader, get_writer_cls_from_file(o_fa_fn)(o_fa_fn) as writer:
        for r in reader:
            writer.writeRecord(r.id, r.sequence)


def run(args):
    fq2fa(args.input_fq_fn, args.output_fa_fn)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
