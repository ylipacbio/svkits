#!/usr/bin/env python
from argparse import ArgumentParser


def get_scan_fasta_parser():
    parser = ArgumentParser("Scan FASTA file and report AT rich region as SV.")
    parser.add_argument("input_fasta", type=str, help="Input FASTA file, must be indexed")
    parser.add_argument("output_bed", type=str, help="Output regions as BED")
    parser.add_argument("--windowsize", default=2000, type=int, help="Scan window size")
    parser.add_argument("--stepsize", default=500, type=int, help="Scan step size")
    parser.add_argument("--slop", default=0, type=int, help="Extend start and end positions in output BED by slop base pairs")
    return parser
