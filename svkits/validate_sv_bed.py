#!/usr/bin/env python
"""Given truth SV calls and your SV calls,
compute number of all SV calls, good SV calls and false positive SV calls.
"""

import sys
import argparse
from pbsv.io.VcfIO import BedRecord, BedReader


class Constant(object):
    MAX_SV_LEN_DIFF_PERCENTAGE = 0.25
    MAX_START_DIFF_LEN = 20


def all_calls(out_bed):
    """return number of all calls."""
    return [r for r in BedReader(out_bed)]


def is_good_bed_record(std_record, out_record,
                       max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                       max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    """Return True if out_record is a good one according to std_record"""
    if out_record.chrom != std_record.chrom or \
        (out_record.sv_type.val != std_record.sv_type.val):
        return False
    is_good_sv_len = (std_record.sv_len * (1.0-max_sv_len_diff_percentage) <= out_record.sv_len and
        out_record.sv_len <= std_record.sv_len * (1.0 + max_sv_len_diff_percentage))
    is_good_start = (abs(std_record.start - out_record.start) <= max_start_diff_len)
    return is_good_sv_len and is_good_start


def good_calls(std_bed, out_bed,
               max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
               max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    """return number of good calls"""
    std_records = [r for r in BedReader(std_bed)]
    out_records = [r for r in BedReader(out_bed)]
    ret = []
    for out_record in out_records:
        is_good = False
        for std_record in std_records:
            if is_good_bed_record(std_record, out_record, max_sv_len_diff_percentage, max_start_diff_len):
                is_good = True
                break
        if is_good:
            ret.append(out_record)
    return ret


def n_all_calls(out_bed):
    return len(all_calls(out_bed))


def n_good_calls(std_bed, out_bed,
                 max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                 max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    return len(good_calls(std_bed, out_bed, max_sv_len_diff_percentage, max_start_diff_len))


def n_false_postives(std_bed, out_bed,
                     max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                     max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    """ False positives = (all calls) - (good calls)"""
    return n_all_calls(out_bed) - n_good_calls(std_bed, out_bed, max_sv_len_diff_percentage, max_start_diff_len)


def false_positive_calls(std_bed, out_bed,
                         max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                         max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    std_records = [r for r in BedReader(std_bed)]
    out_records = [r for r in BedReader(out_bed)]
    ret = []
    for out_record in out_records:
        is_good = False
        for std_record in std_records:
            if is_good_bed_record(std_record, out_record, max_sv_len_diff_percentage, max_start_diff_len):
                is_good = True
                break
        if not is_good:
            ret.append(out_record)
    return ret


def remove_redunant_records(records,
                            max_sv_len_diff_percentage=Constant.MAX_SV_LEN_DIFF_PERCENTAGE,
                            max_start_diff_len=Constant.MAX_START_DIFF_LEN):
    """Return a list of non-redundant BedRecords"""
    ret = []
    for r in records:
        uniq = True
        for p in ret:
            if is_good_bed_record(r, p, max_sv_len_diff_percentage, max_start_diff_len):
                uniq = False
                break
        if uniq:
            ret.append(r)
    return ret


def get_parser():
    """return arg parser"""
    p = argparse.ArgumentParser(description="""Compare out.bed with std.bed, get recall/sensitivity, false positives""")
    p.add_argument("out_bed", help="Your SV calls BED")
    p.add_argument("std_bed", help="Ground truth SV calls BED")
    p.add_argument("--max_sv_len_diff_percentage", type=float, default=Constant.MAX_SV_LEN_DIFF_PERCENTAGE, help="Maximum percentage difference in lengths allowed to merge two SV calls")
    p.add_argument("--max_start_diff_len", type=int, default=Constant.MAX_START_DIFF_LEN, help="Maximum difference of start positions allowed to merge two SV calls")
    return p


def run(args):
    """run main"""
    std_bed = args.std_bed
    out_bed = args.out_bed
    print 'all_calls=%r' % n_all_calls(out_bed=out_bed)
    print 'good_calls=%r' % n_good_calls(std_bed=std_bed, out_bed=out_bed)
    print 'false_positives=%r' % n_false_postives(std_bed=std_bed, out_bed=out_bed)
    return 0

def main():
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
