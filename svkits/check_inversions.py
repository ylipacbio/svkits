#!/usr/bin/env python

"""
check inversions given a BAM file of alignments and a BED of inversions.
Input:
    a BAM file including all alignments, must be sorted and indexed
    a BED file including only inversions
Output:
    a CSV file containing columns below:
      read, read name
      istart, start of an inversion
      iend, end of an inversion
      ilen, length of an inversion
      aligni, does read align to the inversion, 0 no, 1 yes
      alignb, does read align to sequence before the inversion, 0 no, 1 yes
      aligna, does read align to sequence after the inversion, 0 no, 1 yes
      aligni_str, aligned strand of the read to the inversion, excluding before and after alignment
      alignb_str, aligned strand of the read to sequence before the inversion
      aligna_str, aligned strand of the read to sequence after the inversion
      aligni_len, aligned length of the read to the inversion, excluding before and after alignment
      alignb_len, aligned length of the read to sequence before the inversion
      aligna_len, aligned length of the read to sequence after the inversion
      aligni_sim, aligned percentage similarity of the read to the inversion, excluding before and after alignment
      alignb_sim, aligned percentage similarity of the read to sequence before the inversion
      aligna_sim, aligned percentage similarity of the read to sequence after the inversion
"""
import sys
import argparse
from pbsv.io.VcfIO import  BedReader
from pbsv.independent.common import  to_ref_regions
# from pbsv.io.bamstream import BamStream
from pbsv.libs import AlignmentFile
import csv
from collections import namedtuple, defaultdict

def get_parser():
    """return arg parser"""
    p = argparse.ArgumentParser(description="""Report if alignments can report inversions in BED.""")
    p.add_argument("in_bam", help="Input BAM file, must be sorted by chromosome, indexed")
    p.add_argument("in_bed", help="Input BED containing inversions")
    p.add_argument("out_csv", help="Output CSV report")
    return p


ATTRS = ['read', 'istart', 'iend', 'ilen', 'aligni', 'alignb', 'aligna',
         'aligni_str', 'alignb_str', 'aligna_str',
         'aligni_len', 'aligna_len', 'alignb_len',
         'aligni_sim', 'aligna_sim', 'alignb_sim']
DELIMITER = '\t'
Attr = namedtuple('Attr', ATTRS, verbose=False)

def check_align(alnobjs, sv_obj):
    """
    All alnobjs has the same read name.
    """

def yield_alnobjs_groupby_readname(aln_reader, chrom, start, end, slop):
    """Yield (readname, [aln_obj, ..., aln_obj]) tuples, where
     * read must have one or more alignments intersect chrom:(start-slop)-(end+slop) regions.
     * all aln_obj in list are alignments of the read
    """
    ref_lens_d = dict(zip(aln_reader.references, aln_reader.lengths))
    _start = max(0, start)
    _end = min(end, ref_lens_d[chrom])

    d = defaultdict(list) # {readname: [alnobj, ..., alnobj]}
    for aln in aln_reader.fetch(aln_reader, chrom, _start, _end):
        d[aln.query_name].append(aln)

    for readname, alnobjs in d.iteritems():
        yield (readname, alnobjs)


def check_inversion_in_alnobjs_of_a_read(readname, chrom, inv_start, inv_end, alnobjs):
    return Attr(read=readname, istart=inv_start, iend=inv_end, ilen=abs(inv_end - inv_start),
                aligni=0, alignb=0, aligna=0, aligni_str=0, alignb_str=0, aligna_str=0,
                aligni_len=0, alignb_len=0, aligna_len=0, aligni_sim=0, alignb_sim=0, aligna_sim=0)


def _run_checkion(sv_reader, aln_reader, writef):
    for sv in sv_reader:
        assert sv.is_Inversion()
        for readname, alnobjs in yield_alnobjs_groupby_readname(aln_reader, sv.chrom, sv.start, sv.end, slop=1000):
            attr = check_inversion_in_alnobjs_of_a_read(readname, sv.chrom, sv.start, sv.end, alnobjs)
            writef(attr)


def run_checkion(in_bam, in_bed, out_csv):
    """
    """
    sv_reader = BedReader(in_bed)
    aln_reader = AlignmentFile(in_bam)
    csv_writer = csv.writer(open(out_csv, 'w'), delimiter=DELIMITER)

    csv_writer.writerow(DELIMITER.join(ATTRS))
    def writef(attr):
        csv_writer.writerow([attr._asdict()[name] for name in ATTRS])

    _run_checkion(sv_reader, aln_reader, writef)

    sv_reader.close()
    aln_reader.close()
    csv_writer.close()


def run(args):
    """run main"""
    run_checkion(in_bam=args.in_bam, in_bed=args.in_bed, out_csv=args.out_csv)
    return 0


def main():
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
