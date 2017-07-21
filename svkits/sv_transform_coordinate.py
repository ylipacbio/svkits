#!/usr/bin/env python

from argparse import ArgumentParser
from collections import defaultdict
import os
import os.path as op
import sys
import json

from pbsv.independent.utils import is_bed, is_vcf
from pbsv.io.VcfIO import VcfReader, VcfWriter, VcfRecord, BedReader, BedWriter, BedRecord

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser()
    parser.add_argument("input_sv_fn", help="Input BED or VCF filename")
    parser.add_argument("output_sv_fn", help="Output BED or VCF filename, format must be the same as input")
    return parser

def get_reader_cls_from_file(fn):
    if is_bed(fn):
        return BedReader
    elif is_vcf(fn):
        return VcfReader
    else:
        raise ValueError("Could not get reader for %s" % fn)

def get_writer_cls_from_file(fn):
    if is_bed(fn):
        return BedWriter
    elif is_vcf(fn):
        return VcfWriter
    else:
        raise ValueError("Could not get reader for %s" % fn)

def sv_transform_chrom_coordinate(i_fn, o_fn):
    if not (all(is_bed(fn) for fn in [i_fn, o_fn]) or all(is_vcf(fn) for fn in [i_fn, o_fn])):
        raise ValueError("Input and output must be both BED or VCF")
    with get_reader_cls_from_file(i_fn)(i_fn) as reader, get_writer_cls_from_file(o_fn)(o_fn) as writer:
        for r in reader:
            writer.writeRecord(r.transform_chrom_coordinate())

def run(args):
    sv_transform_chrom_coordinate(args.input_sv_fn, args.output_sv_fn)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
