#!/usr/bin/env python
__desc__ = \
"""
Given
    reference genome size G,
    zmw mean subread length L,
    sum of (zmw mean subread length * number of zmws) GX,
    percentage of mappable zmws M%,
simulate mapping zmws to reference genome and compute
    number of reference base pairs that are covered by exactly 0,1,2,3... zmws
"""

import logging
import sys
import argparse
import random
from collections import OrderedDict, Counter
import numpy as np
import os.path as op
from pbcore.util.Process import backticks
from pbcore.io import FastaRecord, FastaWriter, FastaReader
from pbsv.functional.common import SvType
from pbsv.io.VcfIO import BedRecord, BedWriter
from .add_an_indel_to_fasta import fasta_to_ordereddict, get_del_pos, get_ins_pos, del_a_substr_from_read, ins_a_str_to_read

logging.basicConfig(level=logging.DEBUG, format='Log: %(message)s')


def get_parser():
    """return arg parser"""
    desc = __desc__
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("genome_size", metavar='G', help="Reference genome size in base pairs.")
    parser.add_argument("X", help="sum of (mean zmw subread length * number of zmws) over genome size")
    parser.add_argument("mean_srl", metavar='L', help="Zmw mean subread length")
    parser.add_argument("--mappable_zms_percentage", type=int, default=100, help="Percentage of mappable zmws over all zmws")
    parser.add_argument("--out", type=str, default=None, help="Output file")
    return parser

class CovTree(object):

    def __init__(self, G):
        self.G = int(G) # genome size
        self.tree = IntervalTree()
        intv = Interval(0, self.G)
        self.tree.insert(intv)
        #self.cov = dict(self.G: 0)

    def __str__(self):
        return "CovTree(G=%r, nodes=%s)" % (self.G, len(self.tree))

    def add_cov_interval(self, start, end):
        self.tree.find(start, end)



def run(args):
    """Run from args"""
    G, X, L, M = args.genome_size, args.X, args.mean_srl, args.mappable_zmw_percentage
    out = args.out
    logging.info("Locals: %r" % locals())

    number_zmws = X * G / L


def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
