from __future__ import absolute_import
from pbcore.io import openDataSet
from pbsv1.independent.utils import execute, is_bam, is_fofn, is_xml

import sys
import argparse
import logging
import os.path as op

FORMATTER = op.basename(__name__) + ':%(levelname)s:'+'%(message)s'
logging.basicConfig(level=logging.INFO, format=FORMATTER)
log = logging.getLogger(__name__)

"""
Given a bam or a PacBio dataset, using `samtools stats` to count number of reads, bases. etc.
"""

class SamStats(object):
    def __init__(self, filename):
        self.filename = filename
        self.contents = [r.strip() for r in open(filename)]
        self.reads_total, self.reads_mapped, self.reads_unmapped = 0, 0, 0
        self.bases_total, self.bases_mapped = 0, 0
        self._parse()

    def _parse(self):
        ATRS = ['reads_total', 'reads_mapped', 'reads_unmapped',
                'bases_total', 'bases_mapped']
        KEYS = ['SN\traw total sequences:', 'SN\treads mapped:', 'SN\treads unmapped:',
                'SN\ttotal length:', 'SN\tbases mapped:']
        for l in self.contents:
            for idx, (atr, key) in enumerate(zip(ATRS, KEYS)):
                if l.startswith(key):
                    value = int(l.split(':')[1].split('\t')[1])
                    setattr(self, atr, value)

    def first_value_of_key(self, key):
        rs = [l.split(':')[1].split('\t')[1] for l in self.contents if l.startswith(key)]
        if len(rs) == 0:
            raise ValueError("Could not find line starts with {} in {}".format(key, self.filename))
        return rs[0]

    def __str__(self):
        return '<SamStats: total reads={}, mapped={}, unmapped={}; total bases={}, mapped={}>'.format(
                self.reads_total, self.reads_mapped, self.reads_unmapped,
                self.bases_total, self.bases_mapped)


def generate_samstats(bam_files, out_prefix):
    samstats_fns = []
    ss_objs = []
    for idx, bam_file in enumerate(bam_files):
        samstats_fn = '{}.{}.samstats'.format(out_prefix, idx)
        samstats_fns.append(samstats_fn)
        cmd = 'samtools stats {} > {}'.format(bam_file, samstats_fn)
        execute(cmd)
        ss = SamStats(samstats_fn)
        ss_objs.append(ss)
    return bam_files, samstats_fns, ss_objs


def write_bam_stat_csv(out_csv, bam_files, samstats_fns, ss_objs):
    log.info("Writing output csv to {}\n".format(out_csv))
    assert len(bam_files) == len(samstats_fns)
    assert len(bam_files) == len(ss_objs)
    headers = ['bam', 'samstats', 'reads_total', 'reads_mapped', 'reads_unmapped', 'bases_total', 'bases_mapped']
    with open(out_csv, 'w') as writer:
        writer.write(','.join(headers) + '\n')
        for idx in range(len(bam_files)):
            ss = ss_objs[idx]
            fields = [bam_files[idx], samstats_fns[idx],
                      ss.reads_total, ss.reads_mapped, ss.reads_unmapped,
                      ss.bases_total, ss.bases_mapped]
            writer.write(','.join([str(x) for x in fields]) + '\n')


def bam_files_from_fofn_or_dataset(i_fn):
    if is_fofn(i_fn):
        bam_files = []
        in_files = [r.strip() for r in i_fn]
        for in_file in open(in_files):
            if is_bam(in_file):
                bam_files.append(in_file)
            elif is_xml(in_file):
                bam_files.extend(openDataSet(in_file).toExternalFiles())
            else:
                raise ValueError("Expected to see bam or fofn or dataset xml, but get {}".format(in_file))
        return bam_files
    elif is_xml(i_fn):
        return openDataSet(i_fn).toExternalFiles()
    elif is_bam(i_fn):
        return [i_fn]
    else:
        raise ValueError("Expected to see bam or dataset xml, but get {}".format(i_fn))


def parse_sum(ss_objs):
    reads_total = 0
    base_total = 0
    for ss in ss_objs:
        reads_total += ss.reads_total
        base_total += ss.bases_total
    log.info("Total number of reads: {}".format(reads_total))
    log.info("Total number of bases: {}".format(base_total))


def load_samstats(samstats_fns):
    ss_objs = []
    for fn in samstats_fns:
        ss = SamStats(fn)
        ss_objs.append(ss)
    return ss_objs


def should_load_samstats(i_fn):
    if i_fn.endswith('.samstats'):
        return True
    if is_fofn(i_fn):
        return all([r.strip().endswith('.samstats') for r in open(i_fn)])
    return False


def samstats_from_fofn(i_fn):
    assert should_load_samstats(i_fn)
    if is_fofn(i_fn):
        return [r.strip() for r in open(i_fn)]
    elif i_fn.endswith('.samstats'):
        return [i_fn]
    raise ValueError("Could not parse samstats file {}".format(i_fn))


def run_main(i_fn, out_prefix):
    out_csv = out_prefix + '.samstats.csv'
    if should_load_samstats(i_fn):
        samstats_fns = samstats_from_fofn(i_fn)
        ss_objs = load_samstats(samstats_fns)
        bam_files = ['fake' for _ in range(len(ss_objs))]
    else:
        bam_files = bam_files_from_fofn_or_dataset(i_fn)
        bam_files, samstats_fns, ss_objs = generate_samstats(bam_files, out_prefix)

    write_bam_stat_csv(out_csv, bam_files, samstats_fns, ss_objs)
    parse_sum(ss_objs)


def run(args):
    log.info("locals: {!r}".format(locals()))
    run_main(args.input_file, args.out_prefix)


def get_parser():
    """Set up and return argument parser."""
    parser = argparse.ArgumentParser("Quick get sam stats of large bam files or datasets using samtools stats")
    parser.add_argument("input_file", help="Input BAM or FOFN or dataset XML")
    parser.add_argument("out_prefix", help="Output prefix string")
    return parser

def main(args=sys.argv[1:]):
    """main"""
    run(get_parser().parse_args(args))


if __name__ == "__main__":
    sys.exit(main(args=sys.argv[1:]))
