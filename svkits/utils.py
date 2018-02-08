#!/usr/bin/env python

"""
Define utils used in
* make-subreads-bam-from-fasta.py
* make-subreads-bam-from-zmws.py
"""
from collections import defaultdict
import os.path as op
from pbcore.io import FastaReader, FastaWriter, FastqReader, FastqWriter
from pbsv.independent.utils import realpath, rmpath, execute

def mkdir(path):
    return execute('mkdir -p {}'.format(path))

def is_fasta(fn):
    """Return true if a file extension is fa or fasta"""
    return _is_fmt(fn, ["fa", "fasta"])

def is_fastq(fn):
    """Return true if a file extension is fa or fasta"""
    return _is_fmt(fn, ["fq", "fastq"])

def get_reader_cls_from_file(fn):
    """Return FastaReader, FastqReader or raise ValueError"""
    if is_fasta(fn):
        return FastaReader
    elif is_fastq(fn):
        return FastqReader
    else:
        raise ValueError("Could not get reader for %s" % fn)

def get_writer_cls_from_file(fn):
    """Return FastaWriter, FastqWriter or raise ValueError"""
    if is_fasta(fn):
        return FastaWriter
    elif is_fastq(fn):
        return FastqWriter
    else:
        raise ValueError("Could not get reader for %s" % fn)

def autofmt(filename, validfmts, defaultfmt=None):
    """Infer the format of a file from its filename.  As a convention all the
    format to be forced with prefix followed by a colon (e.g. "fmt:filename").

    `validfmts` is a list of acceptable file formats
    `defaultfmt` is the format to use if the extension is not on the valid list

    returns `filename`,`fmt`
    """
    colonix = filename.find(":")
    if colonix != -1:
        extension = filename[:colonix]
        filename = filename[(colonix+1):]
    else:
        extension = None
        for validfmt in validfmts:
            if filename.endswith(validfmt):
                extension = filename[-len(validfmt):]
    return filename, (extension.lower() if extension in validfmts else defaultfmt)

def _is_fmt(fn, validfmts):
    return autofmt(fn, validfmts)[1] in validfmts

def get_movie_and_zmw_from_name(name):
    """Given a string of pacbio zmw name or read name, return movie and zmw"""
    try:
        fs = name.strip().split(' ')[0].split('/')
        movie, zmw = fs[0], fs[1]
        return movie, int(zmw)
    except ValueError:
        raise ValueError("Read %r is not a PacBio read." % name)


def get_movie2zmws_from_zmws(zmws):
    """
    ..doctest:
        d = get_movie2zmws_from_zmws(['movie/0', 'movie1/1'])
        d['movie']
        0
        d['movie1']
        1
    """
    movie2zmws = defaultdict(lambda:set()) #movie --> set of zmws
    for zmw in zmws:
        movie, zmw = get_movie_and_zmw_from_name(zmw)
        movie2zmws[movie].add(zmw)
    return movie2zmws


def get_movie2zmws_in_txt(in_txt):
    """Return movie2zmws {movie: set(zmws)} of zmws or reads in txt"""
    return get_movie2zmws_from_zmws(zmws=[r for r in open(in_txt, 'r')])


def get_movie2zmws_in_fasta(in_fasta):
    """Return movie2zmws {movie: set(zmws)} of reads in fasta"""
    return get_movie2zmws_from_zmws(zmws=[r.name for r in FastaReader(in_fasta)])


def fofn2fns(i_fofn):
    """Get filenames from a fofn or fn"""
    if not i_fofn.endswith('fofn'):
        return [i_fofn]
    fns = []
    for fn in open(i_fofn, 'r'):
        fn = fn.strip()
        if len(fn) == 0 or fn[0] == '#':
            continue
        fns.append(fn)
    return fns


def get_movie2bams_from_fofn(in_bam_fofn):
    """
    Return {movie: bam}
    e.g.,
    if inpupt bam fofn has two bam files of two movies:
    movie1.subreads.bam
    movie2.subreads.bam
    return {'movie1': 'movie1.subreads.bam', 'movie2': 'movie2.subreads.bam'}
    e.g., if input bam fofn has two bam files of the same movie
    movie.1.subreads.bam
    movie.2.subreads.bam
    return {'movie': ['movie.1.subreads.bam', 'movie.2.subreads.bam']}
    """
    movie2bams = defaultdict(lambda:set())
    for fn in fofn2fns(in_bam_fofn):
        movie = fn.split('/')[-1].split('.')[0]
        print "movie=%s" % movie
        movie2bams[movie].add(fn)
    return movie2bams


def rmpath_cmd(path):
    """rm a file path, not directory, command"""
    return 'rm -f %s' % path

def mkdir_cmd(path):
    """mkdir -p path"""
    return 'mkdir -p %s' % path

def cmds_to_bash_fn(cmds, bash_sh_fn):
    """Write a list of cmds to a bash sh"""
    with open(bash_sh_fn, 'w') as writer:
        writer.write('\n'.join(["#!/usr/bin/env bash"] + cmds))
    execute('chmod +x %s' % realpath(bash_sh_fn))
    return realpath(bash_sh_fn)

def merge_bam2xml_cmd(in_bams, out_xml):
    """Merge bam files to xml command """
    cmd = 'dataset merge {out_xml} {in_bams}'.format(out_xml=out_xml, in_bams=' '.join(in_bams))
    return cmd

def consolidate_xml2bam_cmd(in_xml, out_bam, out_xml):
    """Consolidate in_xml to out_bam"""
    cmd = 'dataset consolidate {in_xml} {out_bam} {out_xml}'.format(in_xml=in_xml, out_bam=out_bam, out_xml=out_xml)
    return cmd


def bamsieve_zmw_subreads_cmd(in_subreads_bam, zmws, out_prefix):
    """in_subreads_bam -- input original subreads.bam, from which subreads will be extracted
    zmws -- a list of zmws
    out_prefix -- output prefix
    Return output (cmd, bam)
    """
    wl = ','.join([str(r) for r in list(zmws)])
    out_bam = out_prefix + '.' + op.basename(in_subreads_bam) + ".subset.out.bam"
    cmd = "bamSieve {in_subreads_bam} {out_bam} --whitelist {wl}".format(in_subreads_bam=in_subreads_bam, out_bam=out_bam, wl=wl)
    return cmd, out_bam


def bamsieve_zmw_subreads_cmds(bam2zmws, out_prefix):
    """
    bam2zmws --> {subreads.bam: zmws},
    mapping a subreads.bam file to its white list zmws which will be extracted
    return cmds, out_bams
    """
    cmds, out_bams = [], []
    for in_bam, zmws in bam2zmws.iteritems():
        cmd, out_bam = bamsieve_zmw_subreads_cmd(in_bam, zmws, out_prefix)
        cmds.append(cmd)
        out_bams.append(out_bam)
    return (cmds, out_bams)


def get_bam2zmws(movie2zmws, movie2bams):
    """
    movie2zmws --> {movie: zmws}
    movie2bams --> {movie: bam_files}  # allow *.1.subreads.bam, *.2.subreads.bam
    return bam2zmws{bam: zmws}
    e.g.,
    movie2zmws {'movie': [20,21,31]},
    movie2bams {'movie': ['movie.1.subreads.bam', 'movie.2.subreads.bam']}
    return bam2zmws{'movie.1.subreads.bam': [20,21,31], 'movie.2.subreads.bam': [20,21,31]}
    """
    for movie in movie2zmws.keys():
        if movie not in movie2bams:
            raise ValueError("Could not find subreads bam of movie %s in %r" % (movie, movie2bams.values()))

    bam2zmws = {}
    for movie, zmws in movie2zmws.iteritems():
        in_subreads_bam_files = list(movie2bams[movie])
        for in_subreads_bam in in_subreads_bam_files:
            bam2zmws[in_subreads_bam] = zmws
    return bam2zmws


def execute_cmds(cmds, dry_run):
    """Execute cmds if not dry_run, otherwise,print cmds"""
    for cmd in cmds:
        print cmd
        if not dry_run:
            execute(cmd)


def make_subreads_bam(movie2zmws, movie2bams, out_prefix, dry_run=False):
    """
    movie2zmws --> {movie: [zmw1,zmw2,..]}
    movie2bams --> {movie: [bam1,bam2,..]}
    out_prefix string, e.g., myout
    """
    out_bam = out_prefix + '.subreads.bam'
    out_xml = out_prefix + '.subreadset.xml'
    merged_xml = out_prefix + '.merged.xml'

    # map every bam file to its associated white list zmws
    bam2zmws = get_bam2zmws(movie2zmws=movie2zmws, movie2bams=movie2bams)

    c0 = rmpath_cmd(out_xml)
    c1 = rmpath_cmd(out_bam) # clean up before start

    _cmds, out_bams = bamsieve_zmw_subreads_cmds(bam2zmws=bam2zmws, out_prefix=out_prefix) # bamsieve cmds

    c2 = merge_bam2xml_cmd(in_bams=out_bams, out_xml=merged_xml) # merge cmd
    c3 = consolidate_xml2bam_cmd(in_xml=merged_xml, out_bam=out_bam, out_xml=out_xml) # consolidate cmd
    c4 = rmpath_cmd(merged_xml)

    cmds = [c0, c1] + _cmds + [c2, c3, c4]

    for fn in out_bams:
        cmds.extend([rmpath_cmd(fn), rmpath_cmd(fn+'.pbi')])

    execute_cmds(cmds=cmds, dry_run=dry_run)

    if not dry_run:
        if not op.exists(out_xml):
            raise ValueError("%s does not exist" % out_xml)
        if not op.exists(out_bam):
            raise ValueError("%s does not exist" % out_bam)
    return out_bam


EXTENS = ['.subreadset.xml', '.subreads.bam', '.xml', '.bam', '.fasta', '.fa']


def remove_extension(name):
    """remove known extensions from name"""
    for ext in EXTENS:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def mummer_plot(query_fa, target_fa, out_prefix, min_match_len=50):
    """
    Call mummer and mummerplot to align query_fa to target_fa
    """
    out_ps = "%s.ps" % out_prefix
    out_png = "%s.png" % out_prefix
    print "mummerplot %s vs %s, save to %s|png" % (query_fa, target_fa, out_ps)
    if not op.exists(query_fa) or not op.exists(target_fa):
        raise IOError("mummer inputs %s and %s must exist." % (query_fa, target_fa))

    mums_fn = "%s.mums" % out_prefix
    cmd = "mummer -mum -b -l %s -c %s %s > %s" % (min_match_len, target_fa, query_fa, mums_fn)
    execute(cmd)

    cmd = "mummerplot -f -l -postscript -p %s %s" % (out_prefix, mums_fn)
    execute(cmd)

    cmd = "mummerplot -f -l -png -p %s %s" % (out_prefix, mums_fn)
    execute(cmd)
    return out_ps, out_png

#TODO: use pbsv_polish.utils.*
def bed2prefix(bed_record):
    fields = [bed_record.chrom, bed_record.start, bed_record.end, bed_record.sv_type, bed_record.sv_len]
    return '_'.join([str(x) for x in fields])


#TODO: use pbsv_polish.utils.*
def write_fasta(out_fa_fn, records):
    """Write a list of fasta records [(name, seq), ...,  (name, seq)] to out_fa_fn"""
    from pbcore.io import FastaWriter
    with FastaWriter(out_fa_fn) as w:
        for r in records:
            w.writeRecord(r[0], r[1])

#TODO: use pbsv_polish.utils.*
def _fn2fmtarg(fn):
    fnext2fmtarg = {"m0": "-m 0", "m4": "-m 4", "bam": "--bam"}
    return fnext2fmtarg[autofmt(fn, fnext2fmtarg.keys())[1]]

#TODO: use pbsv_polish.utils.*
def blasr_cmd(query_fn, target_fn, out_fn, nproc=8):
    return "blasr {q} {t} {fmt} --out {out_fn} --nproc {nproc} --maxMatch 15 --bestn 10 --hitPolicy randombest".\
            format(q=query_fn, t=target_fn, out_fn=out_fn, fmt=_fn2fmtarg(out_fn), nproc=nproc)

#TODO: use pbsv_polish.utils.*
def basename_prefix_of_fn(fn):
    """Return base name prefix of a file path, e.g.,
    """
    basename = op.basename(fn)
    return basename[0:basename.rfind('.')] if '.' in basename else basename

def __align(a_fa_fn, b_fa_fn):
    from pbtranscript.io import BLASRM4Reader
    out_m4_fn = op.join(op.dirname(a_fa_fn), basename_prefix_of_fn(a_fa_fn) +'.'+basename_prefix_of_fn(b_fa_fn)+'.m4')
    execute(blasr_cmd(a_fa_fn, b_fa_fn, out_m4_fn))
    alns = [r for r in BLASRM4Reader(out_m4_fn)]
    if len(alns) != 1:
        return None
        #raise ValueError("M4 file %s must contain exactly one alignment"  % out_m4_fn)
    return alns[0]

def rotate_read(read_name, read_seq, break_point):
    new_name = read_name + '__rotate__%s' % (break_point)
    new_read = read_seq[break_point:] +read_seq[0:break_point]
    return new_name, new_read

def get_only_read_in_fa(fa_fn):
    reads = [r for r in FastaReader(fa_fn)]
    if len(reads) != 1:
        raise ValueError("FASTA file %s must contain exactly one read"  % fa_fn)
    return reads[0]

def rotate_only_read_in_fa(fa_fn, break_point):
    read = get_only_read_in_fa(fa_fn)
    new_name, new_read = rotate_read(read_name=read.name, read_seq=read.sequence, break_point=break_point)
    rotate_fn = fa_fn + '.rotate.%s' % (break_point) + '.fasta'
    write_fasta(out_fa_fn=rotate_fn, records=[(new_name, new_read)])
    return rotate_fn


def circular_align(a_fa_fn, b_fa_fn):
    aln = __align(a_fa_fn, b_fa_fn)
    if aln is None:
        return None

    # rotate a_fa_fn sequence
    a_rotate_fn = rotate_only_read_in_fa(a_fa_fn, break_point=aln.qStart)
    # rotate b_fa_fn sequence
    b_rotate_fn = rotate_only_read_in_fa(b_fa_fn, break_point=aln.sStart)

    # align rotated a_fa_fn to rotated b_fa_fn
    raln = __align(a_rotate_fn, b_rotate_fn)
    if raln is None:
        return aln
    return aln if aln.score < raln.score else raln # score the less the better

def m42str(r):
    return "%s %s %s_%s/%s %s_%s/%s"  % (r.identity, r.strand, r.qStart, r.qEnd, r.qLength, r.sStart, r.sEnd, r.sLength)
