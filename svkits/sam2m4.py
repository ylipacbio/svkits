#!/usr/bin/env python

from argparse import ArgumentParser
from pbsv.io.linefile import X2PysamReader
import sys
import os
import os.path as op

def sam2m4(fn):
    reader = X2PysamReader(fn)._alignment_file
    for r in reader:
        strand = 1 if r.is_reverse else 0
        fields = [r.query_name, r.reference_name, 'score', 'similarity', 0, r.query_alignment_start, r.query_alignment_end, r.infer_query_length(), strand,  r.reference_start, r.reference_end, r.reference_length]
        print '\t'.join([str(x) for x in fields])

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("Converta SAM/BAM file to m4 and print")
    parser.add_argument("input_fn", type=str, help="Input SAM or BAM filename")
    return parser


def run(args):
    sam2m4(fn=args.input_fn)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
