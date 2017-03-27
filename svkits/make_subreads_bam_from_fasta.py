#!/usr/bin/env python

"""Make a subreads.bam and subreadset.xml containing
all reads in fasta file.
Input: a fasta file, a fofn of all subreads.bam files
Output: o.subreads.bam and o.subreadset.xml
"""
import sys
import argparse
from .utils import (get_movie2zmws_in_fasta, get_movie2bams_from_fofn, make_subreads_bam, remove_extension)


def get_parser():
    """return arg parser"""
    p = argparse.ArgumentParser(description="""Make a subreads.bam/xml containing all reads of a fasta file""")
    p.add_argument("in_fasta", help="Input FASTA file")
    p.add_argument("in_bam_fofn", help="Input BAM fofn")
    p.add_argument("out_prefix", help="Output prefix")
    p.add_argument("--dry_run", default=False, action="store_true", help="Dry run")
    return p


def make_subreads_bam_of_zmws_in_fasta(in_fasta, in_bam_fofn, out_prefix, dry_run=False):
    """Give a fasta file and a fofn of bam files where reads in fasta file are extracted,
    return subreads.bam file which only contains subreads of zmws in the fasta file.
    e.g., fasta file contains reads ('movie/1/0_100', 'movie2/3/300_400'),
    bam fofn contains ('path_to_movie_subreads_bam', 'path_to_movie2_subreads_bam')
    return output.subreads.bam which contains subreads of zmw 'movie/1' and 'movie2/3',
    for example
        'movie/1/0_100', 'movie/1/150_300', 'movie/1/350_500', ...
        'movie2/3/100_200', 'movie2/3/300_400', 'movie2/3/500_600', ...
    """
    movie2zmws = get_movie2zmws_in_fasta(in_fasta)
    movie2bams = get_movie2bams_from_fofn(in_bam_fofn)
    make_subreads_bam(movie2zmws, movie2bams, out_prefix, dry_run=dry_run)


def run(args):
    """run main"""
    out_prefix = remove_extension(args.out_prefix)
    make_subreads_bam_of_zmws_in_fasta(
        in_fasta=args.in_fasta, in_bam_fofn=args.in_bam_fofn,
        out_prefix=out_prefix, dry_run=args.dry_run)
    return 0


def main():
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
