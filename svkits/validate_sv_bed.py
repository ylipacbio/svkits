#!/usr/bin/env python
"""Given truth SV calls and your SV calls,
compute number of all SV calls, good SV calls and false positive SV calls.
"""

import sys
import argparse
from .miscio import BedRecord, BedReader


class Constant(object):
    MAX_SV_LEN_DIFF_PERCENTAGE = 25
    MAX_START_DIFF_LEN = 100

MAX_DIFF_SCORE = 999999
def get_diff_score(std_record, out_record,
                   max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                   max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    if not is_similar(std_record, out_record,  max_sv_len_diff_percentage, max_start_diff_len):
        return MAX_DIFF_SCORE
    else:
        return abs(std_record.start - out_record.start) + 3 * abs(std_record.sv_len - out_record.sv_len)

def is_similar(std_record, out_record,
               max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
               max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    """Return True if out_record is a good one according to std_record"""
    if out_record.chrom != std_record.chrom or \
        (out_record.sv_type.val != std_record.sv_type.val):
        return False
    assert max_sv_len_diff_percentage >= 1 and max_sv_len_diff_percentage <= 100
    is_good_sv_len = (abs(std_record.sv_len * (1.0-(max_sv_len_diff_percentage/100.0))) <= abs(out_record.sv_len) and
        abs(out_record.sv_len) <= abs(std_record.sv_len * (1.0 + (max_sv_len_diff_percentage/100.0))))
    is_good_start = (abs(std_record.start - out_record.start) <= max_start_diff_len)
    return is_good_sv_len and is_good_start

def remove_redunant_records(records,
                            max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                            max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    """Return a list of non-redundant BedRecords"""
    ret = []
    for r in records:
        uniq = True
        for p in ret:
            if is_similar(r, p, max_sv_len_diff_percentage, max_start_diff_len):
                uniq = False
                break
        if uniq:
            ret.append(r)
    return ret


class CompareSVCalls(object):
    def __init__(self, std_bed_fn, out_bed_fn,
                 max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                 max_start_diff_len=Constant.MAX_START_DIFF_LEN):
        self.std_bed_fn = std_bed_fn
        self.out_bed_fn = out_bed_fn
        self.max_sv_len_diff_percentage = max_sv_len_diff_percentage
        self.max_start_diff_len = max_start_diff_len

        self.std_records = [r for r in BedReader(self.std_bed_fn)]
        self.out_records = [r for r in BedReader(self.out_bed_fn)]

    @property
    def true_positive_call_ids_to_std_call_ids(self):
        """{index_of_true_positive_call: index_of_corresponding std_call}"""
        ret  = {} # index of true positive call: index of similar ground truth call
        for i, out_record in enumerate(self.out_records):
            diff_scores = [MAX_DIFF_SCORE] * len(self.std_records)
            for j, std_record in enumerate(self.std_records):
                diff_scores[j] = get_diff_score(std_record, out_record, self.max_sv_len_diff_percentage, self.max_start_diff_len)
            if not all([diff_score==MAX_DIFF_SCORE for diff_score in diff_scores]):
                ret[i] = diff_scores.index(min(diff_scores))
        return ret

    @property
    def true_positive_call_ids(self):
        return self.true_positive_call_ids_to_std_call_ids.keys()

    @property
    def true_positive_calls(self):
        """return true positive SV calls."""
        return [self.out_records[i] for i in self.true_positive_call_ids_to_std_call_ids]

    @property
    def n_true_positive_calls(self):
        return len(self.true_positive_call_ids)

    @property
    def n_std_calls(self):
        return len(self.std_records)

    @property
    def n_out_calls(self):
        return len(self.out_records)

    @property
    def false_negative_call_ids(self):
        std_ids = set(range(0, self.n_std_calls))
        tp_ids = set(self.true_positive_call_ids_to_std_call_ids.values())
        return sorted(std_ids.difference(tp_ids))

    @property
    def false_negative_calls(self):
        return [self.std_records[i] for i in self.false_negative_call_ids]

    @property
    def n_false_negative_calls(self):
        return self.n_std_calls - self.n_true_positive_calls

    @property
    def false_positive_call_ids(self):
        """ False positives = (all out calls) - (true positive calls)"""
        return sorted(set(range(0, self.n_out_calls)).difference(set(self.true_positive_call_ids)))

    @property
    def false_positive_calls(self):
        return [self.out_records[i] for i in self.false_positive_call_ids]

    @property
    def n_false_positive_calls(self):
        return len(self.false_positive_calls)

    def pretty_true_positive_calls(self):
        """ e.g., chrI    12578   12578   Insertion   386 --> (0, 0) --> chrI 12764   12764   Insertion   378"""
        ret = []
        for i, j in self.true_positive_call_ids_to_std_call_ids.iteritems():
            out, std = self.out_records[i], self.std_records[j]
            start_diff, len_diff = out.start - std.start, out.sv_len - std.sv_len
            s = '%s --> (i=%s, j=%s, start_diff=%s, len_diff=%s) --> %s' % (record2str(out), i, j, start_diff, len_diff, record2str(std))
            if out.seq != '.' and std.seq != '.':
                from .utils import write_fasta, bed2prefix, circular_align, m42str
                out_seq_fn, std_seq_fn = 'tmp.' + bed2prefix(out) + '.svseq.fasta', 'tmp.' + bed2prefix(std) + '.svseq.fasta'
                write_fasta(out_fa_fn=out_seq_fn, records=[(bed2prefix(out), out.seq)])
                write_fasta(out_fa_fn=std_seq_fn, records=[(bed2prefix(std), std.seq)])
                m4_record = circular_align(out_seq_fn, std_seq_fn)
                if m4_record is None:
                    s += ' --> Not mappable'
                else:
                    s += ' --> %s' % (m42str(m4_record))
            ret.append(s)
        return '\n'.join(ret)

    def len_diffs(self):
        ret = []
        for i, j in self.true_positive_call_ids_to_std_call_ids.iteritems():
            out, std = self.out_records[i], self.std_records[j]
            start_diff, len_diff = out.start - std.start, out.sv_len - std.sv_len
            ret.append(abs(len_diff))
        return ret

    def start_diffs(self):
        ret = []
        for i, j in self.true_positive_call_ids_to_std_call_ids.iteritems():
            out, std = self.out_records[i], self.std_records[j]
            start_diff, len_diff = out.start - std.start, out.sv_len - std.sv_len
            ret.append(abs(start_diff))
        return ret


def get_parser():
    """return arg parser"""
    p = argparse.ArgumentParser(description="""Compare out.bed with std.bed, get recall/sensitivity, false positives""")
    p.add_argument("out_bed", help="Your SV calls BED")
    p.add_argument("std_bed", help="Ground truth SV calls BED")
    p.add_argument("--max_sv_len_diff_percentage", type=float, default=Constant.MAX_SV_LEN_DIFF_PERCENTAGE, help="Maximum percentage difference in lengths allowed to merge two SV calls")
    p.add_argument("--max_start_diff_len", type=int, default=Constant.MAX_START_DIFF_LEN, help="Maximum difference of start positions allowed to merge two SV calls")
    return p

def record2str(r):
    return '\t'.join([str(r.chrom), str(r.start), str(r.end), str(r.sv_type), str(r.sv_len)])

def records2str(records):
    return '\n'.join([record2str(r) for r in records])

def run(args):
    """run main"""
    if not (args.max_sv_len_diff_percentage >= 1 and args.max_sv_len_diff_percentage <= 100):
        raise ValueError("--max_sv_len_diff_percentage must be between 1 and 100")

    cmpsv = CompareSVCalls(std_bed_fn=args.std_bed, out_bed_fn=args.out_bed, max_sv_len_diff_percentage=args.max_sv_len_diff_percentage, max_start_diff_len=args.max_start_diff_len)

    print 'n out calls: %r' % cmpsv.n_out_calls
    print 'n std calls: %r' % cmpsv.n_std_calls
    print 'n true positive: %r' % cmpsv.n_true_positive_calls
    print 'n false negative: %r' % cmpsv.n_false_negative_calls
    print 'n false positive: %r' % cmpsv.n_false_positive_calls

    print '---------------\n'
    print 'true positive calls:\n%s\n' % cmpsv.pretty_true_positive_calls()
    print '---------------\n'
    import numpy as np
    from scipy import stats
    print 'Length difference:\n%s\n%s\n' % (cmpsv.len_diffs(), stats.describe(cmpsv.len_diffs()))
    print 'Start difference:\n%s\n%s\n' % (cmpsv.start_diffs(), stats.describe(cmpsv.start_diffs()))
    print '---------------\n'
    print '---------------\n'
    print 'false nagative calls:\n%s\n' % records2str(cmpsv.false_negative_calls)
    print '---------------\n'
    print 'false positive calls:\n%s\n' % records2str(cmpsv.false_positive_calls)
    print '---------------\n'

    return 0

def main():
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
