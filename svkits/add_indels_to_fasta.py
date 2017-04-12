#!/usr/bin/env python

"""
Add indels (insertions and deletions) to a reference FASTA file,
and make BED of SV calls mapping subreads of the reference FASTA file
to the modified FASTA.
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

ARTIFICIAL_FMT = '0/1:3:6'

def to_int_list(s):
    """Convert a comma delimited string to a list of int"""
    return [int(x) for x in s.strip().split(',')]

def get_parser():
    """return arg parser"""
    desc = """Modify a reference fasta by adding indels (insertions
and deletions. Make a truth SV BED file for mapping reads of
the original reference to the modified reference file."""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("in_fasta", help="Input FASTA")
    parser.add_argument("out_fasta", help="Output FASTA with indels")
    parser.add_argument("out_bed", help="Ground truth SV calls in BED if modified fasta were used as reference")

    parser.add_argument("--lens", type=to_int_list, default='50,100',
                        help="Comma delimited string representing indel lengths, e.g., '50,100' means indel lengths are 50, 100")
    parser.add_argument("--numbers", type=to_int_list, default='1,2',
                        help="Comma delimited string representing number of indels for each lengths, e.g., '1,2' means one 50bp indel, and two 100bp indel")
    parser.add_argument("--name_pattern", type=str, default='.+',
                        help="Only work on sequences with name matching name_pattern")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--insert", action='store_true', help="Insert indels")
    group.add_argument("--delete", action='store_true', help="Delete indels")
    return parser


def make_sane(total_l, indel_type, lens, numbers):
    """Sanity check input."""
    def is_match(l, n, t):
        if len(l) != len(n):
            raise ValueError("Indel %r lengths %r does not match numbers %r." % (t, l, n))

    is_match(lens, numbers, indel_type)

    total_bases = sum(l*n for l, n in zip(lens, numbers))

    if total_bases > 0.1 * total_l:
        raise ValueError("Could not insert or delete > 10 percent base pairs of fasta length %s" % (total_l))


def weighted_choose(total_n, p_dict):
    """Choose a total of total_n elements from p_dict.keys(), with
    probability of choosing a key k equals to p_dict[k].
    p_dict is a dict {k: p}, where p is probability of choosing k
    return dict {k: num}
    """
    keys = p_dict.keys()
    r = np.random.choice(keys, p=[p_dict[k] for k in keys], size=total_n)
    c = Counter(r)
    return {k:c[k] for k in keys}

def weighted_choose_multiple(total_n_dict, p_dict):
    """Given total_n_dict = {size: total_n}, p_dict = {k: p}
    return dict = {size-> weighted_choose(total_n, p_dict)}
    """
    ret = dict()
    for k in total_n_dict.keys():
        ret[k] = weighted_choose(total_n=total_n_dict[k], p_dict=p_dict)
    return ret

def swap_dict_dict_k1_k2(dd):
    """Input dd is a dict of dict = {k1: {k2: val}}
    swap keys k1 and k2, return dict = {k2: {k1: val}}
    """
    ret = dict()
    for k1 in dd:
        for k2 in dd[k1]:
            if k2 in ret:
                ret[k2][k1] = dd[k1][k2]
            else:
                ret[k2] = {k1: dd[k1][k2]}
    return ret


def get_del_poses(l_seq, lens, numbers):
    """Get deletion indels' starting positions, must be apart by max(lens) bps"""
    print locals()
    low = int(l_seq * 0.2)
    del_sum = sum([l*n for l,n in zip(lens,numbers)])
    high = l_seq - del_sum
    if low > high:
        raise ValueError("Could not find propoer del posistion, sequence length %r * 0.2 = %r > %r - %r" % (l_seq, low, l_seq, max(lens), high))
    poses = sorted(np.random.randint(low=low, high=high, size=sum(numbers)))
    for index in range(0, len(poses)-1)[::-1]:
        pos, next_pos = poses[index], poses[index+1]
        if next_pos - pos < max(lens): # deletion pos should be sparse
            poses[index] = poses[index] - max(lens)
    if len(poses) > 0 and poses[0] < 0:
       raise ValueError("Could not successfully delete indels (lengths %r, numbers %r) on sequence of length %s, try again!" % (lens, numbers, l_seq))
    return poses


def get_ins_poses(l_seq, lens, numbers):
    """Get insertion indels' starting positions"""
    poses = sorted(np.random.randint(low=0, high=int(l_seq * 0.8), size=sum(numbers)))
    for index in range(1, len(poses)):
        pos, pre_pos = poses[index], poses[index-1]
        if pos - pre_pos < max(lens): # insertion pos should be sparse for convenience of simulation
            poses[index] = poses[index] + max(lens)
        if len(poses) > 0 and poses[len(poses)-1] > l_seq:
            raise ValueError("Could not successfully insert indels (lengths %r, numbers %r) on sequence of length %s, try again!" % (lens, numbers, l_seq))
    return poses


def expand_objects(objs, numbers):
    """e.g., objs=['a', 'b'], numbers = [2,3], return ['a', 'a', 'b', 'b', 'b']"""
    ret = []
    for obj, number in zip(objs, numbers):
        ret.extend([obj for i in range(0, number)])
    return ret

def filter_dict_by_keys(d, key_pattern):
    """Remove items whose keys do not match name_pattern, return filterd dict"""
    ret = OrderedDict()
    import re
    for key in d.keys():
        result = re.match(key_pattern, key)
        if result:
            ret[key] = d[key]
    return ret

def run(args):
    """run main"""
    in_fasta, out_fasta, out_bed = args.in_fasta, args.out_fasta, args.out_bed
    lens, numbers = args.lens, args.numbers
    indel_type = SvType(SvType.Insertion) if args.insert else SvType(SvType.Deletion)
    logging.info("Locals %r" % locals())

    _reads_d = fasta_to_ordereddict(in_fasta)
    reads_d = filter_dict_by_keys(_reads_d, args.name_pattern)
    logging.info("Choosable reads in fasta: %r" % reads_d.keys())

    lens_d = {k:len(reads_d[k]) for k in reads_d.keys()}
    total_l = sum(len(reads_d[k]) for k in reads_d.keys()) # total number of bases in fasta
    # sanity check
    make_sane(total_l=total_l, indel_type=indel_type, lens=lens, numbers=numbers)

    # {sequence/chromosome name: probablity}, e.g., {'chr1': 0.6, 'chr2':0 .2, 'chr3':0.3}
    p_d = {k:(lens_d[k]*1.0/total_l) for k in lens_d.keys()}

    # n_d[chr][size] -> num of indels of size in chromosome
    n_d = swap_dict_dict_k1_k2(weighted_choose_multiple(total_n_dict=dict(zip(lens, numbers)), p_dict=p_d))

    bed_records = []

    fasta_writer = FastaWriter(out_fasta)
    bed_writer =  BedWriter(out_bed)

    for name in n_d.keys():
        # work on sequence/chromosome with name
        new_seq = reads_d[name]
        l_new_seq = len(new_seq)
        new_name = name
        this_lens = n_d[name].keys() # n_d[name]: dict len: number
        this_numbers = [n_d[name][l] for l in this_lens]
        expanded_lens = expand_objects(objs=this_lens, numbers=this_numbers)
        poses = None
        f = None
        logging.info("Selecting sequence %r, sequence length %s, indel lens=%r, indel numbers=%r" % (name, l_new_seq, this_lens, this_numbers))
        if indel_type.is_Deletion: # if indel type is deletion
            poses = get_del_poses(l_seq=l_new_seq, lens=this_lens, numbers=this_numbers)
            logging.info("Delete positions are: %r" % poses)
            f = del_a_substr_from_read
        elif indel_type.is_Insertion: # if indel type is insertion
            poses = get_ins_poses(l_seq=l_new_seq, lens=this_lens, numbers=this_numbers)
            logging.info("Insert positions are: %r" % poses)
            f = ins_a_str_to_read
        else:
            raise ValueError("Unknown indel type %s" % indel_type)

        for l, pos in zip(expanded_lens, poses):
            # add an indel of l bps to sequence/chromosome name
            new_obj, bed_record = f(name=new_name, seq=new_seq, pos=pos, n_bases=l, simple_name=True)
            new_name = new_obj.name
            new_seq = new_obj.sequence
            l_new_seq = len(new_seq)
            bed_records.append(bed_record)

        # write modified sequence
        fasta_writer.writeRecord(FastaRecord(new_name, new_seq))

        # write all bed records
        for bed_record in bed_records:
            bed_writer.writeRecord(bed_record)

    fasta_writer.close()
    bed_writer.close()
    return 0

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
