#!/usr/bin/env python
from __future__ import division

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
from pbsv.libs import AlignmentFile, AlignedSegment
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
Attr1 = namedtuple('Attr1', ['chr', 'istart', 'iend', 'ilen', 'aligni'], verbose=False)


def yield_alnobjs_groupby_readname(aln_reader, chrom, start, end, slop):
    """Yield (readname, [aln_obj, ..., aln_obj]) tuples, where
     * read must have one or more alignments intersect chrom:(start-slop)-(end+slop) regions.
     * all aln_obj in list are alignments of the read
    """
    ref_lens_d = dict(zip(aln_reader.references, aln_reader.lengths))
    _start = max(0, start-slop)
    _end = min(end + slop, ref_lens_d[chrom])

    d = defaultdict(list) # {readname: [alnobj, ..., alnobj]}
    for aln in aln_reader.fetch(chrom, _start, _end):
        d[aln.query_name].append(aln)

    for readname, alnobjs in d.iteritems():
        yield (readname, alnobjs)


#def is_good_aln(alnobj, minlen, chrom, start, end, startwiggle, endwiggle):
#    """This function takes an AlignedSeg as input, returns True if all criteria are satisfied
#    * alnobj's query sequence length is greater than or equals to minlen
#    * alnobj's reference equals to chrom
#    * alnobj's reference start position is within startwiggle bp from start
#    * alnobj's reference end position is within endwiggle bp from end
#    """
#    return ( len(alnobj.query_sequence) >= minlen and
#             alnobj.reference_name == chrom and
#             abs(alnobj.reference_start - start) <= startwiggle and
#             abs(alnobj.reference_end - end) <= endwiggle)


def can_clip_to(alnobj, minlen, chrom, start, end, startwiggle, endwiggle):
    """
    Return True if alnobj can be clipped to chrom:start-end, otherwise False.

    start --- smallest start pos of the clipped alignment, inclusive.
         * If start is None, use alnobj's start position.
         * Otherwise, use max(start, alnobj.reference_start) as start
    end   --- largest end pos of the clipped alignment, exclusive.
         * if end is None, use alnobj's end position.
         * Otherwise, use min(end, alnobj.reference_end) as end
    minlen --- minimum length of clipped alignments, drop short clipped alignments.
    startwiggle --- maximum drift from expected start position, drop far away clipped alignments.
    endwiggle --- maximum drift from expected end position, drop far away clipped alignments.
    """
    _start = max(start, alnobj.reference_start) if start is not None else alnobj.reference_start
    _end   = min(end, alnobj.reference_end) if end is not None else alnobj.reference_end
    return (alnobj.reference_name == chrom and
            (start is not None and abs(_start - start) <= startwiggle) and
            (end is not None and abs(_end - end) <= endwiggle) and
            _end - _start >= minlen)


def clip_alns_to(alnobjs, chrom, start, end, minlen, startwiggle=sys.maxint, endwiggle=sys.maxint):
    """
    Select from AlignedSegment objects (alnobjs), alignments which can clip to chrom:start-end.
    chrom --- reference of the clipped alignment
    start --- smallest start pos of the clipped alignment, inclusive.
         * If start is None, use alnobj's start position.
         * Otherwise, use max(start, alnobj.reference_start) as start
    end   --- largest end pos of the clipped alignment, exclusive.
         * if end is None, use alnobj's end position.
         * Otherwise, use min(end, alnobj.reference_end) as end
    minlen --- minimum length of clipped alignments, drop short clipped alignments.
    startwiggle --- maximum drift from expected start position, drop far away clipped alignments.
    endwiggle --- maximum drift from expected end position, drop far away clipped alignments.
    """
    ret = []
    for alnobj in alnobjs:
        if can_clip_to(alnobj, minlen, chrom, start, end, startwiggle, endwiggle):
            print start, end
            _start = max(start, alnobj.reference_start) if start is not None else alnobj.reference_start
            _end   = min(end, alnobj.reference_end) if end is not None else alnobj.reference_end
            print alnobj.reference_start, alnobj.reference_end, _start, _end
            clipped_alnobj = _clip_aln_to(alnobj, _start, _end)
            ret.append(clipped_alnobj)
    return ret


def __find_aln(alnobjs, get_value_f, max_or_min_f):
    if len(alnobjs) > 0:
        values = [get_value_f(alnobj) for alnobj in alnobjs]
        return alnobjs[values.index(max_or_min_f(values))]
    else:
        return None

def find_longest_aln(alnobjs):
    """return aln of which query_sequence is the longest."""
    return __find_aln(alnobjs, lambda alnobj: alnobj.query_alignment_end - alnobj.query_alignment_start, max)

def find_closet_before_aln(alnobjs, pos):
    """return aln of which reference end is closest to pos"""
    return __find_aln(alnobjs, lambda alnobj: abs(alnobj.reference_end - pos), min)

def find_closet_after_aln(alnobjs, pos):
    """return aln of which reference start is closest to pos"""
    return __find_aln(alnobjs, lambda alnobj: abs(alnobj.reference_start - pos), min)

def get_query_offset_from_cigar(cigartuples, reference_offset):
    """
    ..doctest:
        >>> [get_query_offset_from_cigar(( (4,249), (0, 6), (2, 1), (0, 6), (1,1), (0, 16)), i) for i in range(0, 15)]
        [0, 1, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 11, 12, 14]
    """
    ref_offset = reference_offset
    query_offset = 0
    from pbsv.aln.cigar import DEL, INS, MATCH, REF_SKIP, QUERYOPS, REFOPS
    for op, l in cigartuples:
        diff = ref_offset if ref_offset < l else l
        if op in QUERYOPS:
            query_offset += diff
        if op in REFOPS:
            ref_offset = ref_offset - diff
        elif op == REF_SKIP:
            raise NotImplementedError('Not implemented for REF_SKIP')
        else:
            continue
        if ref_offset == 0:
            break
    if ref_offset != 0:
        raise ValueError("reference_offset {} exceed this alignment's cigar operations!".format(reference_offset))
    return query_offset


#def clip_cigar_to(cigartuples, ref_startoffset, ref_endoffset):
#    ref_offset = ref_startoffset
#    index = 0
#    ret = ()
#    for op, l in cigartuples:
#        diff = ref_offset if ref_offset < l else l
#        if op in QUERYOPS:
#            query_offset += diff
#        if op in REFOPS:
#            ref_offset = ref_offset - diff
#        elif op == REF_SKIP:
#            raise NotImplementedError('Not implemented for REF_SKIP')
#        else:
#            continue
#        if ref_offset == 0:
#            break
#    return ret

class valn(object):
    def __init__(self, query_start, query_end, reference_start, reference_end, is_reverse):
        self.query_start = query_start
        self.query_end = query_end

        self.query_alignment_start = query_start
        self.query_alignment_end = query_end

        self.reference_start = reference_start
        self.reference_end =  reference_end
        self.is_reverse = is_reverse

def _clip_aln_to(alnobj, reference_start, reference_end):
    """
    reference_start -- clip alnobj to reference start
    reference_end -- clip alnobj to reference end
    """
    assert alnobj.reference_start <= reference_start
    assert alnobj.reference_end >= reference_end

    r_startoff, r_endoff = reference_start - alnobj.reference_start, reference_end - alnobj.reference_start
    q_startoff, q_endoff = get_query_offset_from_cigar(alnobj.cigartuples, r_startoff),  get_query_offset_from_cigar(alnobj.cigartuples, r_endoff)
    return valn(q_startoff + alnobj.query_alignment_start, q_endoff + alnobj.query_alignment_start, r_startoff, r_endoff, alnobj.is_reverse)

    ret = AlignedSegment()
    ret.query_name = alnobj.query_name
    ret.query_sequence = alnobj.query_sequence[q_startoff:q_endoff]
    ret.reference_id = alnobj.reference_id
    ret.reference_start = reference_start

    ret.cigar = clip_cigar_to(cigartuples, start, end)
    assert ret.reference_end == reference_end

    ret.next_reference_id = alnobj.next_reference_id
    ret.next_reference_start =  alnobj.next_reference_start
    ret.template_length = alnobj.tempalte_length
    ret.mapping_quality = alnobj.mapping_quality
    if not alnobj.query_qualities:
        ret.query_qualities = None
    else:
        from pysam.libcalignedsegment import array_to_qualitystring
        ret.query_qualities = array_to_qualitystring(alnobj.query_qualities)[q_startoff:q_endoff]

    QV_TAGS = ['ip']
    tags = alnobj.tags[::]
    for index, (tag_name, tag_val) in enumerate(tags):
        if tag_name in QV_TAGS:
            if alnobj.query_alignment_length != len(tag_val):
                raise ValueError("Tag {} length {} != query sequence length {}".format(tag_name, len(tag_val), alnobj.query_alignment_length))
    ret.tags = tags
    return tags


def check_inversion_in_alnobjs_of_a_read(readname, chrom, inv_start, inv_end, alnobjs):
    """
    """
    ilen = abs(inv_end - inv_start)
    if len(alnobjs) >= 3:
        aligni=True
    else:
         aligni=False
    return Attr1(chr=chrom, istart=inv_start, iend=inv_end, ilen=ilen, aligni=aligni)


    # Find alignments to inversion, before inversion and after inversion.
    # =====> before inversion ====> <====inversion====  ====> after inversion====>
    alnobj_i = find_longest_aln(clip_alns_to(alnobjs, chrom, inv_start, inv_end, 0.8*ilen, startwiggle=100, endwiggle=100))
    alnobj_b = find_closet_before_aln(clip_alns_to(alnobjs, chrom, None, inv_start, 100, endwiggle=500), inv_start)
    alnobj_a = find_closet_after_aln(clip_alns_to(alnobjs, chrom, inv_end, None, 100, startwiggle=500), inv_end)

    aligni, alignb, aligna = alnobj_i is not None, alnobj_b is not None, alnobj_a is not None

    def _strand(alnobj):
        """Forward strand 0, reverse strand 1"""
        return -1 if alnobj is None else 1 if alnobj.is_reverse else 0

    def _sim(alnobj):
        return 0
        return alnobj.get_tag(tag='NM') / alnobj.query_alignment_length

    def f(g):
        return (g(alnobj_i), g(alnobj_b), g(alnobj_a))

    aligni_str, alignb_str, aligna_str = f(_strand) # _strand(alnobj_i), _strand(alnobj_b), _strand(alnobj_a)
    aligni_len, alignb_len, aligna_len = f(lambda alnobj: alnobj.query_alignment_end-alnobj.query_alignment_start if alnobj else -1) # len(alnobj_i.query_sequence), len(alnobj_b.query_sequence), len(alnobj_a.query_sequence)
    aligni_sim, alignb_sim, aligna_sim = f(_sim)

    return Attr(read=readname, istart=inv_start, iend=inv_end, ilen=ilen,
                aligni=aligni, alignb=alignb, aligna=aligna,
                aligni_str=aligni_str, alignb_str=alignb_str, aligna_str=aligna_str,
                aligni_len=aligni_len, alignb_len=alignb_len, aligna_len=aligna_len,
                aligni_sim=aligni_sim, alignb_sim=alignb_sim, aligna_sim=aligna_sim)


def _check(sv_reader, aln_reader, writef):
    """
    For each sv obj, fetch srounding alignments from aln_reader,
    group alignments by read name, and for each read name, analyze
    if sv is detected by aligner using attributes in Attr, finally
    call writef to write reports.

    sv_reader --- BedReader, each time yield a sv obj
    aln_reader --- pysam.AlignmentFile
    writef --- a function which writes an Attr obj.
    """
    for sv in sv_reader:
        assert sv.sv_type.is_Inversion
        print repr(sv)
        for readname, alnobjs in yield_alnobjs_groupby_readname(aln_reader, sv.chrom, sv.start, sv.start+sv.sv_len, slop=1000):
            attr = check_inversion_in_alnobjs_of_a_read(readname, sv.chrom, sv.start, sv.start+sv.sv_len, alnobjs)
            writef(attr)


def check(in_bam, in_bed, out_csv):
    """
    For each inversion SV in in_bed, check srounding alignments,
    writer check report to out_csv.

    Input:
    in_bam --- sorted, index bam filename
    in_bed --- inversions in BED
    Output
    out_csv --- inversion check report in csv, columns described in ATTR
    """
    sv_reader = BedReader(in_bed)
    aln_reader = AlignmentFile(in_bam)
    csv_writer = csv.writer(open(out_csv, 'w'), delimiter=DELIMITER)

    #csv_writer.writerow(ATTRS)
    #def writef(attr):
    #    csv_writer.writerow([attr._asdict()[name] for name in ATTRS])

    csv_writer.writerow(['istart', 'iend', 'ilen', 'aligni'])
    def writef1(attr):
        csv_writer.writerow([attr.chr, attr.istart, attr.iend, attr.ilen, attr.aligni])

    _check(sv_reader, aln_reader, writef1)

    sv_reader.close()
    aln_reader.close()


def run(args):
    """run main"""
    check(in_bam=args.in_bam, in_bed=args.in_bed, out_csv=args.out_csv)
    return 0


def main():
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
