#!/usr/bin/env python

"""
Add exactcly one indel (either insertion or deletion) to
a reference FASTA file, and make BED of SV calls mapping
subreads of the reference FASTA file to the modified FASTA.
"""

import sys
import argparse
import os.path as op
from pbcore.util.Process import backticks
from pbcore.io import FastaReader, FastaWriter
from pbsv.run import _mktemp
from .utils import mummer_plot, realpath, rmpath


def get_parser():
    """return arg parser"""
    desc = """Mumer dot plot."""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("query_fasta", help="Query FASTA")
    parser.add_argument("target_fasta", help="Target FASTA")
    parser.add_argument("out_prefix", help="Output Prefix")
    parser.add_argument("-l", "--min-match-len", default=20, help="mummer -l")
    parser.add_argument("--query_region", help="Query region to plot, asumming there is only one read in query, e.g., 0-100")
    parser.add_argument("--target_region", help="Target region to plot region, asumming there is only one read in target, e.g., 20-100")
    return parser

def web_path(fn):
    prefix = '/home/UNIXHOME/'
    if realpath(fn).startswith(prefix):
        return 'file://mp-nac01/unixhome/' + realpath(fn)[len(prefix):]
    else:
        return '/Volumes/' + realpath(fn)


def subset_fasta(in_fa, out_fa, region):
    """Assuming there is exactly one read in in_fa, region must be a tuple
    of (start, end), extract read[start:end] and save to out_fa
    """
    from pbsv.libs import Fastafile
    fileobj = Fastafile(in_fa)
    if ':' in region:
        chrom, start_end = region.split(':')
        start, end = [int(x) for x in start_end.split('-')]
    else:
        chrom, start, end = region, 0, dict(zip(fileobj.references, fileobj.lengths))[region]

    try:
        seq = fileobj.fetch(str(chrom), int(start), int(end))
    except Exception:
        raise ValueError("Could not get substring (%s, %s, %s) from %s" % (chrom, start, end, fileobj.filename))

    with FastaWriter(out_fa) as writer:
        writer.writeRecord('%s:%s-%s' % (chrom, start, end), seq)


def run(args):
    query_fasta, target_fasta, out_prefix = args.query_fasta, args.target_fasta, args.out_prefix
    query_region = args.query_region if args.query_region is not None else None
    target_region = args.target_region if args.target_region is not None else None

    new_query_fasta, new_target_fasta = query_fasta, target_fasta
    tmpdir = _mktemp()

    if query_region is not None:
        new_query_fasta = op.join(tmpdir, op.basename(query_fasta)) + '.' + args.query_region + '.fasta'
        subset_fasta(query_fasta, new_query_fasta, query_region)

    if target_region is not None:
        new_target_fasta = op.join(tmpdir, op.basename(target_fasta)) + '.' + args.target_region + '.fasta'
        subset_fasta(target_fasta, new_target_fasta, target_region)

    out_ps, out_png = mummer_plot(new_query_fasta, new_target_fasta, out_prefix, min_match_len=args.min_match_len)
    print "Dot plot of %r %r vs %r %r:\n%r\n%r\n" % (query_fasta, query_region, target_fasta, target_region, out_ps, out_png)
    print "Web Path: %r" % web_path(out_png)

    rmpath(tmpdir)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
