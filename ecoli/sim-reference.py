#!/usr/bin/env python
import sys
import os.path as op
import logging
from pbsv.functional.utils import mkdir, realpath, execute
from pbsv.io.VcfIO import BedReader

from svkits.utils import execute_cmds, rmpath_cmd, mkdir_cmd, cmds_to_bash_fn

from svkits.validate_sv_bed import n_all_calls, n_good_calls, n_false_postives, false_positive_calls, remove_redunant_records

log = logging.getLogger(op.basename(__file__))

makefile_fn = op.join(op.dirname(__file__), 'Makefile')
global_xml = op.join(op.dirname(__file__), 'global_options.xml')
sv_options_xml = op.join(op.dirname(__file__), 'sv_options.xml')


def make_references_and_bed(in_ref_fa, actions, sv_lens, times, output_dir='sim_reference', dry_run=False):
    """Make multiple modified referece and truth BED files,
    return (out_fa_fns, out_bed_fns)
    """
    mkdir(output_dir)
    out_fa_fns, out_bed_fns = [], []
    idx = 0
    for action in actions:
        for length in sv_lens:
            for time in range(0, times):
                out_ref_dir = op.join(realpath(output_dir), 'action_%s_%s_id_%s' % (action.replace('-', ''), length, idx))
                mkdir(out_ref_dir)
                out_fa = op.join(out_ref_dir, 'reference.fasta')
                out_bed = op.join(out_ref_dir + 'truth.bed')
                out_fa_fns.append(out_fa)
                out_bed_fns.append(out_bed)
                cmd = 'add-an-indel-to-fasta %s %s %s %s --n_bases %s' % (in_ref_fa, out_fa, out_bed, action, length)
                execute_cmds([cmd], dry_run)
                idx += 1
    return out_fa_fns, out_bed_fns

def refset_xml_of_ref_fa(reference_fa):
    """Given a reference fasta, return file path to its referenceset.xml"""
    output_dir = op.join(op.dirname(realpath(reference_fa)))
    name = 'reference'
    refset_xml = op.join(output_dir, name, 'referenceset.xml')
    return refset_xml


def make_fasta_to_reference_cmd(reference_fa):
    """Make cmds to run a reference set from a fasta,
    return referenceset.xml, cmds"""
    log.info("locals=%r" % locals)
    exe = '/pbi/dept/secondary/builds/mainline/current_smrttools_fromsrc_installdir/smrtcmds/bin/fasta-to-reference'
    output_dir = op.join(op.dirname(realpath(reference_fa)))
    name = 'reference'

    c0 = mkdir_cmd(output_dir)
    c1 = rmpath_cmd(op.join(output_dir, name))
    c2 = ' '.join([exe, realpath(reference_fa), output_dir, name])
    cmds = [c0, c1, c2]

    refset_xml = refset_xml_of_ref_fa(reference_fa)
    return refset_xml, cmds

def qsub_make_referenceset_cmd(reference_fa):
    fa_dir = op.dirname(reference_fa)
    reference_ds, cmds = make_fasta_to_reference_cmd(reference_fa)
    sh_fn = op.join(op.dirname(reference_fa), 'make-referenceset.sh')
    cmds_to_bash_fn(cmds, sh_fn)
    execute('chmod +x %s' % sh_fn)
    execute('qu 1 %s' % sh_fn)


def make_pbsmrtpipe_dir(reference_fa, subreadset_xml, dry_run=False):
    """Run pbsmrtpipe.pipelines.sa3_ds_sv using out_fa as reference set"""
    # first make a reference set from a fasta
    fa_dir = op.dirname(reference_fa)
    reference_ds = refset_xml_of_ref_fa(reference_fa)

    # make pbsmrtpipe running dir
    pbsmrtpipe_dir = realpath(op.join(fa_dir, 'pbsmrtpipe'))
    c0 = 'rm -rf %s' % pbsmrtpipe_dir
    c1 = mkdir_cmd(pbsmrtpipe_dir)
    # copy *.xml, Makefile to pbsmrtpipe_dir
    c2 = 'cp %s %s %s %s' % (makefile_fn, global_xml, sv_options_xml, pbsmrtpipe_dir)

    execute_cmds([c0, c1, c2], dry_run)

    # make job.sh
    job_fn = op.join(pbsmrtpipe_dir, 'job.sh')
    job_cmd = \
            'cwd=`pwd` && cd %s\n' % pbsmrtpipe_dir + \
            'pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_sv --debug \\\n' + \
            '-e eid_subread:%s \\\n' % subreadset_xml + \
            '-e eid_ref_dataset:%s \\\n' % reference_ds + \
            '--preset-xml=%s \\\n' % op.join(pbsmrtpipe_dir, op.basename(sv_options_xml)) + \
            '--preset-xml=%s && echo $?\n' % op.join(pbsmrtpipe_dir, op.basename(global_xml)) + \
            'cd $cwd'
    cmds_to_bash_fn([job_cmd], job_fn)
    log.info("job.sh = %s" % job_fn)
    return job_fn

def pbsv_out_bed_of_ref_fa(reference_fa):
    return op.join(op.dirname(reference_fa), 'pbsmrtpipe', 'tasks', 'pbsvtools.tasks.call-0', 'structural_variants.bed')

def make_validate_sv_bed(pbsv_out_bed, truth_bed):
    """validate pbsv_out_bed against truth_bed"""
    n_all = n_all_calls(out_bed=pbsv_out_bed)
    n_good = n_good_calls(std_bed=truth_bed, out_bed=pbsv_out_bed)
    n_fp = n_false_postives(std_bed=truth_bed, out_bed=pbsv_out_bed)
    id = op.basename(truth_bed).split('.')[0] # e.x.,action_delete_50_id_0truth
    sv_type, sv_len = id.split('_')[1:3]
    ret = [id, n_all, n_good, n_fp, sv_type, sv_len, pbsv_out_bed, truth_bed]
    print '\t'.join([str(x) for x in ret])



def run():
    actions = ['--delete', '--insert']
    sv_lens = [50, 100, 200, 500, 1000, 2000, 5000]
    times = 5
    dry_run = False
    ecoli_fa = '/pbi/dept/secondary/siv/yli/sv/ecoli/Ecoli/genome/genome.fa'
    subreadset_xml = '/pbi/dept/secondary/siv/yli/sv/ecoli/subreads/ecoli_10X.subreadset.xml'
    run_pbsmrtpipe_jobs_sh = 'run_pbsmrtpipe_jobs.sh'
    run_validate_sv_bed_sh = 'run_validate_sv_bed.sh'

    # modified references and truth bed
    out_fa_fns, truth_bed_fns = make_references_and_bed(in_ref_fa=ecoli_fa, actions=actions, sv_lens=sv_lens, times=times, output_dir='sim_reference', dry_run=True)

    if False:
        for out_fa, truth_bed in zip(out_fa_fns, truth_bed_fns):
            qsub_make_referenceset_cmd(out_fa)
        return 0

    # make pbsmrtpipe job.sh for each modified reference, and write job.sh path to run_pbsmrtpipe_jobs.sh
    if False:
        job_sh_fns = []
        for out_fa, truth_bed in zip(out_fa_fns, truth_bed_fns):
            job_sh = make_pbsmrtpipe_dir(reference_fa=out_fa, subreadset_xml=subreadset_xml, dry_run=dry_run)
            job_sh_fns.append(job_sh)
        cmds_to_bash_fn(job_sh_fns, run_pbsmrtpipe_jobs_sh)
        return 0

    # evaluate calls done by pbsmrtpipe.sa3_ds_sv
    headers = ['id', 'n_all_calls', 'n_good_calls', 'false_positives', 'sv_type', 'sv_len', 'pbsv_bed', 'truth_bed']
    print '#' + '\t'.join(headers)
    fp_set = []
    for out_fa, truth_bed in zip(out_fa_fns, truth_bed_fns):
        pbsv_out_bed = pbsv_out_bed_of_ref_fa(reference_fa=out_fa)
        make_validate_sv_bed(pbsv_out_bed, truth_bed)
        fp_set.extend(false_positive_calls(out_bed=pbsv_out_bed, std_bed=truth_bed))
    uniq_set = remove_redunant_records(fp_set)
    import pdb
    pdb.set_trace()
    return 0

if __name__ == '__main__':
    sys.exit(run())
