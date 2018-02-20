#!/usr/bin/env python

"""
Scan fasta and report AT rich regions as BED.
"""

import sys
import logging
import itertools
from pbsv.libs import Fastafile
from pbsv.independent import utils
from pbsv.independent.common import SvFmts
from pbsv.io.VcfIO import BedWriter, BedRecord
from .independent import system
from .argsutils import get_scan_fasta_parser

ARTIFICIAL_FMT = '0/1:3:6'
ARTIFICIAL_SAMPLE = 'ScanFasta'
ARTIFICIAL_FMTS = SvFmts.fromDict({ARTIFICIAL_SAMPLE: ARTIFICIAL_FMT})


log = logging.getLogger(__name__)

def scan_fasta(fasta_obj, bedwriter, condition_f, windowsize, stepsize, slop):
    """
    fasta_obj --- pysam.Fastafile
    bedwriter --- pbsv.io.VcfIO.BedWriter
    condition_f --- condition_f(seq_in_window) will return
         * True if a window satisfy a condition
         * False, otherwise
    windowsize --- scan fasta window size
    stepsize --- move window step size
    slop --- Exten start and end positions in BED to both ends by slop base pairs
    """
    for reference in fasta_obj.references:
        seq = fasta_obj.fetch(reference)
        for i in range(0, (len(seq) - windowsize + 1) // stepsize + 2):
            start = i * stepsize
            end = min(len(seq), start + windowsize)
            if start < len(seq) and end <= len(seq):
                if condition_f(seq[start:end]):
                    start = max(0, start - slop)
                    end  = min(end + slop, len(seq))
                    bedobj = BedRecord(chrom=reference.split(' ')[0], start=start, end=end, sv_id='.', sv_type='Inversion',
                                       sv_len=end-start, alt=None, annotations=[], fmts=ARTIFICIAL_FMTS)
                    bedwriter.writeRecord(bedobj)


NUC = 'ACGTN'
DINUC = itertools.combinations_with_replacement(NUC, 2)

from collections import defaultdict

def count_dinuc(seq):
    ret = defaultdict(int)
    for i in range(0, len(seq)-1):
        ret[seq[i:i+2]] += 1
    return ret

def run(args):
    log.info('local args: {}'.format(locals()))
    print('local args: {}'.format(locals()))
    fasta_obj = Fastafile(args.input_fasta)
    bedwriter = BedWriter(args.output_bed, [ARTIFICIAL_SAMPLE])
    def condition_f(seq):
        return count_dinuc(seq)['AT'] >= len(seq) * 2 * (1.0/16)
    scan_fasta(fasta_obj=fasta_obj, bedwriter=bedwriter,
               condition_f=condition_f, windowsize=args.windowsize,
               stepsize=args.stepsize, slop=args.slop)

get_parser = get_scan_fasta_parser

def main(args=sys.argv[1:]):
    """main"""
    run(get_parser().parse_args(args))


if __name__ == "__main__":
    sys.exit(main(args=sys.argv[1:]))
