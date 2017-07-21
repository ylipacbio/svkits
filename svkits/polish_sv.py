#!/usr/bin/env python
import os.path as op
import sys
import os
from collections import defaultdict

from pbcore.io import DataSet, FastaWriter

from pbsv.libs import Fastafile
from pbsv.io.linefile import X2PysamReader, iter_within_ref_regions
from pbsv.io.VcfIO import BedReader
from pbsv.independent.common import RefRegion
from pbsv.independent.utils import execute, realpath
from pbsv.cli import _mkdir
from pbsv.independent.utils import execute_as_bash


def get_aln_reader(aln_fn, bed_fn):
    # reader = get_aln_reader(aln_fn=aln_fn, bed_fn=bed_fn)
    ref_regions = get_ref_regions_from_bed_file(bed_fn)
    reader = X2PysamReader(aln_fn, ref_regions)
    return reader

def yield_alns_from_bed_file(alnfile_obj, bedreader_obj):
    for bed_record in bedreader_obj:
        ref_region = to_ref_region(bed_record)
        yield (bed_record, get_alns_within_ref_region(alnfile_obj, ref_region))

def to_ref_region(bed_record):
    return RefRegion(bed_record.chrom, bed_record.start, bed_record.end+1)

def get_ref_regions_from_bed_file(bed_fn):
    return [to_ref_region(bed_r) for bed_r in BedReader(bed_fn)]

def yield_ref_region_from_bed_file(bed_fn):
    for bed_r in BedReader(bed_fn):
        yield to_ref_region(bed_r)

def get_alns_within_ref_region(alnfile_obj, ref_region):
    """Return a list of alignments within a ref region."""
    return [aln for aln in iter_within_ref_regions(alnfile_obj, [ref_region])]

def zmw_from_subread(subread):
    """Given a subread 'movie/zmw/start_end', return 'movie/zmw'"""
    return '/'.join(subread.split('/')[0:2])

def get_query_zmws_from_alns(alns):
    """Given a list of alignments, return a list of non-redudant query zmws"""
    return list(set([zmw_from_subread(sr) for sr in get_query_subreads_from_alns(alns)]))

def get_query_subreads_from_alns(alns):
    """Given a list of alignments, return a list of non-redudant query subreads"""
    return list(set([aln.query_name for aln in alns]))

from svkits.utils import get_movie2zmws_from_zmws, make_subreads_bam

def make_subreads_bam_of_zmws(movie2bams, zmws, out_prefix, dry_run=False):
    movie2zmws = get_movie2zmws_from_zmws(zmws)
    return make_subreads_bam(movie2zmws, movie2bams, out_prefix, dry_run=dry_run)

def get_subreads_bam_files_from_xml(in_subreads_xml):
    return [fn for fn in DataSet(in_subreads_xml).toExternalFiles() if fn.endswith('subreads.bam')]

def get_movies2bams_from_subreads_xml(in_subreads_xml):
    subreads_fns = get_subreads_bam_files_from_xml(in_subreads_xml)
    return get_movies2bams_from_subreads_bam_files(subreads_fns)

def get_movies2bams_from_subreads_bam_files(subreads_bam_fns):
    """Input: a list of subreads.bam files.
       Return {movie: bam_fn}
    """
    movie2bams = defaultdict(lambda:set())
    for fn in subreads_bam_fns:
        movie = fn.split('/')[-1].split('.')[0]
        print "movie=%s" % movie
        if movie in movie2bams and not fn == movie2bams[movie]:
            raise ValueError("Movie %s mapping to multiple bam files %r and %r" % (movie, movie2bams[movie], fn))
        movie2bams[movie].add(fn)
    return movie2bams


def make_consensus_script_of_subreads_bam(subreads_bam, o_script_fn, o_consensus_id, o_consensus_fn):
    """Make a script file which when excuted, make a consensus sequence of subreads_bam"""
    cmds = []
    output_prefix = op.join(op.dirname(o_consensus_fn), 'sv_pbdagcon')
    c0 = 'sv_pbdagcon %s %s %s' % (subreads_bam, output_prefix, o_consensus_id)
    output_dagcon_fasta = output_prefix + '_ref.fasta'
    align_bam = op.join(op.dirname(o_consensus_fn), 'sr2_sv_pbdagcon.bam')
    nproc = 16
    c1 = 'blasr %s %s --nproc 16 --hitPolicy randomBest --minMatch 8 --maxMatch 15 --bam --out %s' % (subreads_bam,output_dagcon_fasta, align_bam)
    c2 = 'samtools faidx %s' % (output_dagcon_fasta)
    c3 = 'samtools sort {f1} -o {f2} && mv {f2} {f1} && samtools index {f1} && pbindex {f1}'.format(f1=align_bam, f2=align_bam+'.sorted.bam')
    output_prefix = o_consensus_fn[0:o_consensus_fn.rfind('.')]
    hqlq_out_fa, hqlq_out_fq = output_prefix + '.hqlq.fasta', output_prefix + '.hqlq.fastq'
    out_fa, out_fq = output_prefix + '.fasta', output_prefix + '.fastq'
    c4 = "variantCaller --algorithm=best {aln_bam} --verbose -j{nproc} --minMapQV {minqv} --referenceFilename={ref_fa} -o {out_fa} -o {out_fq}".format(aln_bam=align_bam, ref_fa=output_dagcon_fasta, out_fa=hqlq_out_fa, out_fq=hqlq_out_fq, nproc=nproc, minqv=10)
    c5 = 'trim_lq %s %s --min_qv 30' % (hqlq_out_fq, out_fq) # simply remove lower case sequences on both ends
    c6 = 'fq2fa %s %s' % (out_fq, out_fa)
    cmds = [c0, c1, c2, c3, c4, c5, c6]

    print 'Running %s' % o_script_fn
    execute_as_bash(cmds, o_script_fn)


def make_script_of_pbsv_run(reads_fn, ref_fasta_fn, cfg_fn, o_sv_fn, o_script_fn):
    cmds = []
    tmp_sv_fn = op.join(op.dirname(o_sv_fn), 'use_substr_as_chrom.%s' % op.basename(o_sv_fn))
    c0 = 'pbsv run %s %s %s --cfg_fn %s' % (ref_fasta_fn, reads_fn, tmp_sv_fn, cfg_fn)
    c1 = 'python sv_transform_chrom_coordinate %s %s' % (tmp_sv_fn, o_sv_fn)
    cmds = [c0, c1]
    print 'Running %s' % o_script_fn
    execute_as_bash(cmds, o_script_fn)


def write_fasta(o_fasta_fn, records):
    """Write a list of fasta records [(name, seq), ...,  (name, seq)] to o_fasta_fn"""
    with FastaWriter(o_fasta_fn) as w:
        for r in records:
            w.writeRecord(r[0], r[1])

def substr_fasta(fileobj, chrom, start, end, o_fasta_fn):
    """fetch a substring of reference fasta sequence and save to output fasta file o_fasta_fn"""
    try:
        seq = fileobj.fetch(chrom, start, end)
    except Exception as e:
        raise ValueError("Could not get substring (%s, %s, %s) from %s" % (chrom, start, end, fileobj.filename))
    name = '%s/%s_%s' % (chrom, start, end)
    write_fasta(o_fasta_fn, [(name, seq)])


if __name__ == "__main__":
    #in_dir = "in_sl_228_hg00733_10fold"
    in_dir = 'in_sl_1813_yeast_10fold'
    aln_fn = op.join(in_dir, "alignments.bam")
    subreads_xml_fn = op.join(in_dir, "subreads.xml")
    genome_fa = op.join(in_dir, "genome.fa")
    bed_fn = op.join(in_dir, 'structural_variants.bed')
    bed_fn = op.join(in_dir, 'chrV_116286_116286_Insertion_5899.bed')
    out_dir = 'out'
    reference_fasta_obj = Fastafile(genome_fa)
    REFERENCE_EXTENSION = 200000

    ofile_obj = open('coverage.txt', 'w')
    alnfile_obj = X2PysamReader(aln_fn)._alignment_file
    bedreader_obj = BedReader(bed_fn)
    movie2bams = get_movies2bams_from_subreads_xml(subreads_xml_fn)
    print movie2bams

    POLISH_CFG_FN = '/pbi/dept/secondary/siv/yli/sv/consensus/pbsv.polish.cfg'
    make_bam = False
    i = 0
    for bed_record, alns in yield_alns_from_bed_file(alnfile_obj, bedreader_obj=bedreader_obj):
        i += 1
        if i % 1000 == 0:
            print "i=%s, got %s covering alignemnts for sv %s" % (i, len(alns), bed_record)
        srs = get_query_subreads_from_alns(alns)
        zmws = get_query_zmws_from_alns(alns)
        ofile_obj.write('%s\t%s\t%s\n' % (len(srs), len(zmws), bed_record))

        if True:
            # make a subdirectory
            sv_str = '_'.join([str(x) for x in [bed_record.chrom, bed_record.start, bed_record.end, bed_record.sv_type, bed_record.sv_len]])
            data_dir = realpath(op.join(out_dir, sv_str))
            _mkdir(data_dir)
            out_prefix = op.join(data_dir, 'in')
            subreads_bam = make_subreads_bam_of_zmws(movie2bams=movie2bams, zmws=zmws, out_prefix=out_prefix, dry_run=(not make_bam))
            script_fn = op.join(data_dir, 'make_consensus.sh')
            consensus_fn = op.join(data_dir, 'polished.fasta')
            make_consensus_script_of_subreads_bam(subreads_bam=subreads_bam, o_script_fn=script_fn, o_consensus_id=sv_str, o_consensus_fn=consensus_fn)

            ref_start, ref_end = max(0, bed_record.start - REFERENCE_EXTENSION), bed_record.end + REFERENCE_EXTENSION
            ref_fasta_fn = op.join(data_dir, 'sv_reference_w_extension.fasta')
            substr_fasta(fileobj=reference_fasta_obj, chrom=bed_record.chrom, start=ref_start, end=ref_end, o_fasta_fn=ref_fasta_fn)

            polished_sv_fn = op.join(data_dir, 'polished.sv.bed')
            polish_sv_script_fn = op.join(data_dir, 'pbsv_run_polish.sh')
            make_script_of_pbsv_run(reads_fn=consensus_fn, ref_fasta_fn=ref_fasta_fn, cfg_fn=POLISH_CFG_FN, o_sv_fn=polished_sv_fn, o_script_fn=polish_sv_script_fn)

            if i == 1000:
                break

    ofile_obj.close()
