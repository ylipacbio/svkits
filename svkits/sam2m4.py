#!/usr/bin/env python

from argparse import ArgumentParser
from pbsv.io.linefile import X2PysamReader
import sys
import os
import os.path as op

def sam2m4(fn):
    reader = X2PysamReader(fn)._alignment_file
    for r in reader:
        fields = [r.query_name, r.reference_name, r.query_alignment_start, 'score', 'similarity', r.query_alignment_start, r.query_alignment_end, r.infer_query_length(),  r.reference_start, r.reference_end, r.reference_length]
        print '\t'.join(fields)

#def f2():
#    from pbtranscript.io import GMAPSAMReader
#   p0 = [r for r in GMAPSAMReader(fn)][0]
#    print p0

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser("Trim LQ sequences on both ends, where average MapQV of LQ sequences are less than min_qv.")
    parser.add_argument("input_fn", help="Input SAM or BAM filename")
    return parser


def run(args):
    sam2m4(args.input_fn)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
