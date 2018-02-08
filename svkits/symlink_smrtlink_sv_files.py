#!/usr/bin/env python

"""
Make symlink of sv files to a local directory, including
alignments.bam, alignments.bam.bai, reference genome.fa, genome.fa.fai, structural_variants.bed
"""

from argparse import ArgumentParser
import os
import os.path as op
import sys
import json
from pbcore.io import DataSet
from pbsv.ngmlrmap import ln_cmd
from pbsv.independent.utils import execute, realpath
from pbsv.cli import _mkdir

def ln(src, dst):
    execute(ln_cmd(src, dst))


class SMRTLinkSVFiles(object):
    def __init__(self, smrtlink_job_dir):
        self.root_dir = realpath(smrtlink_job_dir)

    @property
    def entry_points_json(self):
        return op.join(self.root_dir, 'workflow', 'entry-points.json')

    @property
    def reference_xml(self):
        ds = json.load(open(self.entry_points_json, 'r'))
        return str(ds['eid_ref_dataset'])

    @property
    def subreads_xml(self):
        ds = json.load(open(self.entry_points_json, 'r'))
        return str(ds['eid_subread'])

    @property
    def reference_fa(self):
        return DataSet(self.reference_xml).toExternalFiles()[0]

    @property
    def ngmlr_index_1(self):
        return self.reference_fa + '-enc.2.ngm'

    @property
    def ngmlr_index_2(self):
        return self.reference_fa + '-ht-13-2.2.ngm'

    @property
    def alignment_bam(self):
        return op.join(self.root_dir, 'tasks/pbsvtools.tasks.gather_align-1/alignments.bam')

    @property
    def alignment_bam_bai(self):
        return op.join(self.root_dir, 'tasks/pbsvtools.tasks.gather_align-1/alignments.bam.bai')

    @property
    def structural_variants_bed(self):
        return op.join(self.root_dir, 'tasks/pbsvtools.tasks.sort_sv-0/structural_variants.bed')

    @property
    def structural_variants_vcf(self):
        return op.join(self.root_dir, 'tasks/pbsvtools.tasks.sort_sv-0/structural_variants.vcf')


def make_symlinks(sl_dir, out_dir):
    _mkdir(out_dir)
    slfiles = SMRTLinkSVFiles(sl_dir)
    ln(slfiles.subreads_xml, op.join(out_dir, 'subreads.xml'))
    ln(slfiles.alignment_bam, op.join(out_dir, 'alignments.bam'))
    ln(slfiles.alignment_bam_bai, op.join(out_dir, 'alignments.bam.bai'))
    ln(slfiles.reference_fa, op.join(out_dir, 'genome.fa'))
    ln(slfiles.ngmlr_index_1, op.join(out_dir, 'genome.fa-enc.2.ngm'))
    ln(slfiles.ngmlr_index_2, op.join(out_dir, 'genome.fa-ht-13-2.2.ngm'))

    execute('samtools faidx %s' % (op.join(out_dir, 'genome.fa')))
    ln(slfiles.structural_variants_bed, op.join(out_dir, 'structural_variants.bed'))

    with open(op.join(out_dir, '_README'), 'w') as writer:
        writer.write('SMRTLink SV job dir: %s\n' % sl_dir)
        writer.write('Subreads: %s\n' % slfiles.subreads_xml)

def get_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser()
    parser.add_argument("smrtlink_job_dir", help="SMRTLink job directory")
    parser.add_argument("out_dir", help="Output directory")
    return parser

def run(args):
    print args
    make_symlinks(sl_dir=args.smrtlink_job_dir, out_dir=args.out_dir)

def main():
    """main"""
    sys.exit(run(get_parser().parse_args(sys.argv[1:])))

if __name__ == "__main__":
    main()
