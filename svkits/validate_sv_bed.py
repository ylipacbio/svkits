#!/usr/bin/env python

# All BED files must be coordinate-sorted.
# Consider variants to overlap if:
#    - same type (i.e. both insertions or both deletions)
#    - similar SV length (+/- 25%)
#    - left and right breakpoints agree (+/- 500 bp)

# Introduced variants: chm1.truth.bed
# Call set: chm1.calls.bed
# Recall / sensitivity
import sys
import argparse
from pbcore.util.Process import backticks

"""
# Truth set: hg38.HG00733-truth.bed
# Call set: hg38.HG00733-calls.bed
# Recall / sensitivity
bedtools closest -k 10 -a hg38.HG00733-truth.bed -b hg38.HG00733-calls.bed |
    awk '($4==$12) { print $0; }' |
    awk '{ n=($5>$13)?$13:$5; x=($5>$13)?$5:$13; if ((x-n)/n <= 0.25) { print $0; } }' |
    awk '{ l=($2>$10)?($2-$10):($10-$2); r=($3>$11)?($3-$11):($11-$3); if(l<=500 && r<=500) { print $0; } }' |
    cut -f 1-8 | uniq | wc -l
"""



def _recall_cmd(std_bed, out_bed):
    """Return cmd to compute recall / sensitivity"""
    cmds = [
        """bedtools closest -k 10 -a {std_bed} -b {out_bed}""".format(std_bed=std_bed, out_bed=out_bed),
        """awk '($4==$12) { print $0; }'""",
        """awk '{ n=($5>$13)?$13:$5; x=($5>$13)?$5:$13; if ((x-n)/n <= 0.25) { print $0; } }'""",
        """awk '{ l=($2>$10)?($2-$10):($10-$2); r=($3>$11)?($3-$11):($11-$3); if(l<=500 && r<=500) { print $0; } }'""",
        """cut -f 1-8 | uniq | wc -l"""
    ]
    return ' | '.join(cmds)


def _all_calls_cmd(out_bed):
    """Return cmd to get number of ground truth calls"""
    cmd = """grep -v '^#' {out_bed} | wc -l""".format(out_bed=out_bed)
    return cmd


def _good_calls_cmd(std_bed, out_bed):
    """Return cmd to compute good calls."""
    cmds = [
        """bedtools closest -k 10 -a {out_bed} -b {std_bed}""".format(out_bed=out_bed, std_bed=std_bed),
        """awk '($4==$12) { print $0; }'""",
        """awk '{ n=($5>$13)?$13:$5; x=($5>$13)?$5:$13; if ((x-n)/n <= 0.25) { print $0; } }'""",
        """awk '{ l=($2>$10)?($2-$10):($10-$2); r=($3>$11)?($3-$11):($11-$3); if(l<=500 && r<=500) { print $0; } }'""",
        """cut -f 1-8 | uniq | wc -l # good calls"""
    ]
    return ' | '.join(cmds)


def execute_output(cmd, output_index, out_type):
    """execute a cmd and return output of index output_index"""
    o, dummy_c, dummy_m = backticks(cmd)
    return out_type(o[output_index])


def recall(std_bed):
    """return recall"""
    cmd = _recall_cmd(std_bed=std_bed)
    return execute_output(cmd=cmd, output_index=0, out_type=int)


def all_calls(out_bed):
    """return number of all calls."""
    cmd = _all_calls_cmd(out_bed)
    return execute_output(cmd=cmd, output_index=0, out_type=int)


def good_calls(std_bed, out_bed):
    """return number of good calls"""
    cmd = _good_calls_cmd(std_bed=std_bed, out_bed=out_bed)
    return execute_output(cmd=cmd, output_index=0, out_type=int)


def false_postives(std_bed, out_bed):
    """ False positives = (all calls) - (good calls)"""
    n_all_calls = all_calls(out_bed=out_bed)
    n_good_calls = good_calls(std_bed=std_bed, out_bed=out_bed)
    return n_all_calls - n_good_calls


def get_parser():
    """return arg parser"""
    p = argparse.ArgumentParser(description="""Compare out.bed with std.bed, get recall/sensitivity, false positives""")
    p.add_argument("out_bed", help="Your SV calls BED")
    p.add_argument("std_bed", help="Ground truth SV calls BED")
    return p


def run(args):
    """run main"""
    std_bed = args.std_bed
    out_bed = args.out_bed
    print 'all_calls=%r' % all_calls(out_bed=out_bed)
    print 'good_calls=%r' % good_calls(std_bed=std_bed, out_bed=out_bed)
    print 'false_positives=%r' % false_postives(std_bed=std_bed, out_bed=out_bed)
    return 0

def main():
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
