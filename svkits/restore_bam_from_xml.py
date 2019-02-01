#!/usr/bin/env python

"""
Given a dataset xml, generate bam.
Input: a dataset xml, e.g, (movie.subreadset.xml)
Output: output.bam
"""
import sys
import argparse
from pbcore.io import DataSet, SubreadSet
from .independent import system


def get_parser():
    """return arg parser"""
    p = argparse.ArgumentParser(description="""Restore bam from dataset xml.""")
    p.add_argument("in_xml", help="Input pacbio dataset xml file")
    p.add_argument("out_bam", help="Output bam")
    return p


def restore_bam_from_xml(in_xml, out_bam):
    ds = SubreadSet(in_xml)
    ds.consolidate(out_bam, numFiles=1, useTmp=False)
    ds.close()
    system.execute('pbindex {}'.format(out_bam))


def run(args):
    """run main"""
    in_xml = args.in_xml
    out_bam = args.out_bam
    restore_bam_from_xml(in_xml=in_xml, out_bam=out_bam)
    return 0


def main():
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
