#!/usr/bin/env python

"""
Add exactcly one indel (either insertion or deletion) to
a reference FASTA file, and make BED of SV calls mapping
subreads of the reference FASTA file to the modified FASTA.
"""

import sys
import argparse
import random
from collections import OrderedDict
import os.path as op
from pbcore.util.Process import backticks
from pbcore.io import FastaRecord, FastaWriter, FastaReader
from pbsv.io.VcfIO import BedRecord, BedWriter


ARTIFICIAL_FMT = '0/1:3:6'

def get_parser():
    """return arg parser"""
    desc = """Modify a reference fasta by either inserting a substring to
a reference sequence or deleting a substring from a reference sequence,
make a truth SV BED file for mapping reads of the original reference to
the modified reference file."""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("in_fasta", help="Input FASTA")
    parser.add_argument("out_fasta", help="Output FASTA with indels")
    parser.add_argument("out_bed", help="Ground truth SV calls in BED if modified fasta were used as reference")
    parser.add_argument("--n_bases", type=int, default=50, help="Number of bases to insert to or delete from a reference sequence")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--insert", action='store_true', help="Insert a string to a sequence in in_fasta, and make a BED file with an Deletion SV call")
    group.add_argument("--delete", action='store_true', help="Delete a substring from a sequence in in_fasta, and make a BED file with an Insertion SV call")
    return parser


def del_a_substr_from_seq(seq, pos, n_bases):
    """Add a deletion of n_bases starting from pos to sequence"""
    if pos < 0 or pos + n_bases > len(seq):
        raise ValueError("Could not delete %r bases starting from pos %r in sequence of length %r" % (n_bases, pos, len(seq)))
    return seq[0:pos] + seq[pos+n_bases:]


def del_a_substr_from_read(name, seq, pos, n_bases, simple_name=False):
    """Delete n_bases bases starting from pos of sequence (seq) of read (name)
    return out_read, bed_record
    """
    out_seq = del_a_substr_from_seq(seq=seq, pos=pos, n_bases=n_bases)
    out_name = name + ' ' + 'ORIGINALLEN=%s;ACTION=DEL:POS=%s;SVLEN=%s;NEWLEN=%s' % (len(seq), pos, n_bases, len(out_seq))
    if simple_name:
        out_name = name.split(' ')[0] + ' ' + ('' if len(name.split(' ')) == 1 else name.split(' ')[1]) + ';d_%s_%s' % (pos, n_bases)
    out_read = FastaRecord(out_name, out_seq)
    bed_record = BedRecord(chrom=name.split(' ')[0], start=pos, end=pos,
                           sv_type='Insertion', sv_len=n_bases, seq=None,
                           fmt=ARTIFICIAL_FMT, annotations=None)
    return (out_read, bed_record)


def make_a_dna_str(n_bases):
    """Make a dna string of length n_bases"""
    seq = ''
    for i in range(0, n_bases):
        seq += random.choice(['A', 'T', 'G', 'C'])
    return seq


def add_a_str_to_seq(seq, pos, n_bases):
    """Simply add an arbitrary string of n_bases to seq at pos"""
    assert isinstance(seq, str)
    s = make_a_dna_str(n_bases)
    out_seq = seq[0:pos] + s + seq[pos:]
    return out_seq

def ins_a_str_to_read(name, seq, pos, n_bases, simple_name=False):
    """Insert a random string of n_bases to pos of sequence (seq) of read (name)
    return out_read, bed_record
    """
    assert isinstance(seq, str)
    out_seq = add_a_str_to_seq(seq=seq, pos=pos, n_bases=n_bases)
    out_name = name + ' ' + 'ORIGINALLEN=%s;ACTION=INS:POS=%s;SVLEN=%s;NEWLEN=%s' % (len(seq), pos, n_bases, len(out_seq))
    if simple_name:
        out_name = name.split(' ')[0] + ' ' + ('' if len(name.split(' ')) == 1 else name.split(' ')[1]) + ';i_%s_%s' % (pos, n_bases)
    out_read = FastaRecord(out_name, out_seq)
    bed_record = BedRecord(chrom=name.split(' ')[0], start=pos, end=pos+n_bases,
                           sv_type='Deletion', sv_len=n_bases, seq=None,
                           fmt=ARTIFICIAL_FMT, annotations=None)
    return (out_read, bed_record)


def fasta_to_ordereddict(in_fasta):
    """Convert a fasta file to OrderedDict([read_name: read_sequence])"""
    d = OrderedDict()
    for r in FastaReader(in_fasta):
        d[r.name] = r.sequence
    return d


def get_del_pos(seq_len, n_del_bases):
    """Get random pos to delete n_del_bases from sequence of length seq_len"""
    return random.randint(0, seq_len - n_del_bases)


def get_ins_pos(seq_len):
    return random.randint(0, seq_len)


def run(args):
    """run main"""
    in_fasta, out_fasta, out_bed = args.in_fasta, args.out_fasta, args.out_bed
    apply_insertion, apply_deletion = args.insert, args.delete
    n_bases = args.n_bases

    reads_d = fasta_to_ordereddict(in_fasta)
    selected_name = random.choice(reads_d.keys())
    selected_seq = reads_d[selected_name]
    out_read, bed_record = None, None
    if apply_deletion:
        if n_bases >= len(selected_seq):
            print "Unable to delete %s bases from a sequence of length %s" % (n_bases, len(selected_seq))
        pos = get_del_pos(len(selected_seq), n_bases)
        out_read, bed_record = del_a_substr_from_read(selected_name, selected_seq, pos, n_bases)
    elif apply_insertion:
        pos = get_ins_pos(len(selected_seq))
        out_read, bed_record = ins_a_str_to_read(selected_name, selected_seq, pos, n_bases)
    else:
        raise ValueError("%s can either insert a string to a reference sequence or delete a substring from a referece sequence." % op.basename(__file__))

    with BedWriter(out_bed) as bed_writer:
        bed_writer.writeRecord(bed_record)

    with FastaWriter(out_fasta) as fasta_writer:
        for name in reads_d:
            if name == selected_name:
                fasta_writer.writeRecord(out_read)
            else:
                fasta_writer.writeRecord(FastaRecord(name, reads_d[name]))
    return 0

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
