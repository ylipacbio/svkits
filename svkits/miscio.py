#!/usr/bin/env python
"""pbsv VCF|BED IO.
Note classes are NOT for generic Vcf|BED IO.
"""
import datetime
import traceback
import re
from collections import defaultdict
from pbsv.independent.iobase import WriterBase, ReaderBase
from pbsv.independent.common import SvType, SvFmt, SvAnnot, SvAnnotations


__all__ = ["VcfRecord",
           "VcfReader",
           "VcfWriter",
           "BedRecord",
           "BedReader",
           "BedWriter"]

PBSV_VCF_META = \
"""##fileformat=VCFv4.3
##fileDate={today}
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVANN,Number=.,Type=String,Description="Repeat annotation of structural variant">
##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Per-sample read depth of this structural variant">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at this position for this sample">
""".format(today=datetime.date.today().strftime('%Y%m%d')) # fileData=YYYYMMDD

PBSV_VCF_HEADER = \
"""#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"""

PBSV_BED_HEADER = \
"""#chrom\tchromStart\tchromEnd\tsvType\tsvLen\tseq\tINFO\tAnnotation\n"""

__WHITESPACE_PATTERN__ = re.compile(r'\s+')
def _str2chrom(s):
    """VCF4.2 CHROM - chromosome: An identifier from the reference genome or an angle-bracketed ID String (<ID>)
    pointing to a contig in the assembly file (cf. the ##assembly line in the header). All entries for a specific
    CHROM should form a contiguous block within the VCF file. The colon symbol (:) must be absent from all
    chromosome names to avoid parsing errors when dealing with breakends. (String, no white-space permitted,
    Required)."""
    s = re.sub(__WHITESPACE_PATTERN__, '_', s)
    s = s.replace(':', '_')
    return s


def _2str(s):
    """Convert s to a string. Return '.' if s is None or an empty string."""
    return '.' if s is None or str(s) == '' else str(s)


def _get_annotations(annotations):
    """Return a list of annotations. If input is None or [], return [SvAnnot.OTHER]"""
    if isinstance(annotations, SvAnnotations):
        return annotations
    if annotations is None or len(annotations) == 0:
        return SvAnnotations([SvAnnot.OTHER])
    else:
        return SvAnnotations([SvAnnot(a) for a in annotations])


def _get_sv_type(sv_type, sv_annot):
    """Return SvType. SAT-1283, only allow SvType = DEL, INS or INV, no DUP"""
    return SvType(sv_type)


def _get_alt(sv_type, sv_annot):
    """Return vcf column 4 alt string e.x., <DEL>, <INS>, <DEL:ME:L1>, <INS:ME:ALU>
    SAT-1283, NO SVTYPE=DUP """
    sv_type = SvType(sv_type)
    sv_annot = SvAnnot(sv_annot)
    if sv_annot.is_ME:
        return '<%s:ME:%s>' % (sv_type.to_vcf_str(), sv_annot)
    else:
        return '<%s>' % (sv_type.to_vcf_str())


def _parse_info(info):
    """parse info, e.x., IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-105;SVANN=OTHER
    Return (sv_type, sv_len, end, sv_annot) """
    try:
        d = defaultdict(str)
        for x in info.split(';'):
            if len(x) > 0 and '=' in x:
                k, v = x.split('=')
                d[k] = v
        return (SvType(d['SVTYPE']), int(d['SVLEN']), int(d['END']), SvAnnotations.fromString(d['SVANN']))
    except (ValueError, KeyError) as e:
        raise ValueError('info %r must have SVTYPE, END, SVLEN and SVANN: %r' %(info, e))

class VcfRecord(object):
    """Vcf Object
    2 321682    . T <DEL>   6 PASS    IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-105;SVANN=OTHER  GT:AD:DP 0/1:3:12
    """
    def __init__(self, chrom, pos, ref, alt, fmt, annotations, end, sv_type, sv_len):
        self.annotations = _get_annotations(annotations)
        self.end = int(end)
        self.sv_type = _get_sv_type(sv_type=sv_type, sv_annot=self.annotations[0])
        self.sv_len = int(sv_len)

        self.chrom = _str2chrom(chrom) # column 0
        self.pos = int(pos) # column 1, start pos of sv in reference
        self.id = '.' # column 2, no id at this point
        self.ref = ref # column 3
        self.alt = alt # column 4, single alt, can either be a sequence, or a symbolic allele e.g., <DEL>
        if alt is None: # infer alt if exact sequence is not provided
            self.alt = _get_alt(sv_type=self.sv_type, sv_annot=self.annotations[0])
        self.qual = '.' # column 5, no qual at this point
        self.filt = 'PASS' # column 6, no filter at this point
        # column 7 info
        # column 8 SvFmt.FMT, GT:AP:DP
        self.fmt = SvFmt.fromString(str(fmt)) # column 9

    @property
    def info(self):
        """vcf column 7 info. e.x., 'IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-105;SVANN=L1' or 'IMPRECISE;SVTYPE=DEL;END=321887;SVLEN=-105;'"""
        infos = ['IMPRECISE', 'SVTYPE=%s' % self.sv_type.to_vcf_str(),
                 'END=%s' % self.end, 'SVLEN=%s' % self.sv_len, "%s" % self.annotations.to_vcf_str()]
        return ';'.join([x for x in infos if len(x) != 0])

    def __str__(self):
        fs = [self.chrom, self.pos, self.id, _2str(self.ref), self.alt, self.qual, self.filt, self.info, SvFmt.FMT, self.fmt]
        return '\t'.join([str(x) for x in fs])

    def __repr__(self):
        return self.__str__()

    @classmethod
    def fromString(cls, s):
        """Create a VcfRecord from a string"""
        fs = s.split('\t')
        if len(fs) < 10:
            raise ValueError("Vcf record {!r} must have 10 columns".format(s))
        if fs[8] != SvFmt.FMT:
            raise ValueError("Format %s is not %s" % (fs[8], SvFmt.FMT))
        sv_type, sv_len, end, annotations = _parse_info(fs[7])
        # support old format in which sv_len of DEL SV is positive
        if sv_type.is_Deletion:
            sv_len = -abs(int(sv_len))
        return VcfRecord(chrom=fs[0], pos=fs[1], ref=fs[3], alt=fs[4], fmt=fs[9],
                         annotations=annotations, end=end,
                         sv_type=sv_type, sv_len=sv_len)


class VcfWriter(WriterBase):
    """ pbsv vcf writer."""
    def __init__(self, f):
        super(VcfWriter, self).__init__(f)
        self._write_meta_lines()
        self._write_comment_lines()

    def _write_meta_lines(self):
        """Write meta info lines to out file."""
        self.file.write(PBSV_VCF_META)

    def _write_comment_lines(self):
        """Write comments/headers to out file."""
        self.file.write(PBSV_VCF_HEADER)

    def writeRecord(self, record):
        """Write Vcf record."""
        if not isinstance(record, VcfRecord):
            raise ValueError("record type %s is not VcfRecord." % type(record))
        else:
            self.file.write("%s\n" % str(record))


class VcfReader(ReaderBase):
    """
    pbsv Vcf Reader
        >>> from pbsv.io import VcfReader
        >>> for record in VcfReader(filename):
        ...     print record
    """
    def __iter__(self):
        for line in self.file:
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                yield VcfRecord.fromString(line)


class BedRecord(object):
    """Bed Object
    chr1\t0\t100\tDeletion\t-100\t.\tGT:AD:DP\t0/1:3:6\t.
    chr1\t1001\t1001\tInsertion\t5\tAAAAA\tGT:AD:DP\t1/1:6:6\tSVANN=ALU
    """
    def __init__(self, chrom, start, end, sv_type, sv_len, seq, fmt, annotations=None):
        self.chrom = chrom # column 0, no need to convert chr1 to 1
        self.start = int(start) # column 1
        self.end = int(end) # column 2
        self.sv_type = SvType(sv_type) # column 3
        self.sv_len = int(sv_len) # column 4
        self.seq = seq # column 5, None, '.' or ({ATGCatgc}+)
        self.info = SvFmt.FMT # column 6, GT:AD:DP
        self.fmt = SvFmt.fromString(str(fmt)) #column 7, 0/1:3:6
        self.annotations = _get_annotations(annotations) # column 8, optional, . or ALU/L1/SVA/TANDEM

    def __str__(self):
        fs = [self.chrom, self.start, self.end, self.sv_type, self.sv_len,
              _2str(self.seq), self.info, self.fmt, self.annotations.to_bed_str()]
        return '\t'.join([str(x) for x in fs])

    @classmethod
    def fromString(cls, s):
        """Create a BedRecord from a string"""
        fs = s.split('\t')
        sv_len = abs(int(fs[4]))
        # support old format in which sv_len of DEL SV is positive
        if 'insertion' in fs[3].lower():
            sv_type = "Insertion"
        elif 'deletion' in fs[3].lower():
            sv_type = "Deletion"
        else:
            raise "SvType Error %s" % sv_type
        return BedRecord(chrom=fs[0], start=fs[1], end=fs[2], sv_type=sv_type,
                         sv_len=sv_len, seq='.', fmt='0/1:3:5', annotations=None)


class BedWriter(WriterBase):
    """ pbsv bed writer."""
    def __init__(self, f):
        super(BedWriter, self).__init__(f)
        self._write_comment_lines()

    def _write_comment_lines(self):
        """Write comments/headers to out file."""
        self.file.write(PBSV_BED_HEADER)

    def writeRecord(self, record):
        """Write Bed record."""
        if not isinstance(record, BedRecord):
            raise ValueError("record type %s is not BedRecord." % type(record))
        else:
            self.file.write("%s\n" % str(record))


class BedReader(ReaderBase):
    """
    pbsv Bed Reader
        >>> from pbsv.io import BedReader
        >>> for record in BedReader(filename):
        ...     print record
    """
    def __iter__(self):
        lineno = 0
        for line in self.file:
            lineno += 1
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                try:
                    fs = line.split()
                    def f(s):
                        return 'Insertion' if 'ins' in s.lower() else 'Deletion' if 'del' in s.lower() else 'None'
                    def g(fs):
                        return '.' if len(fs)<=5 or '#' in fs[5] else fs[5]
                    yield BedRecord(chrom=fs[0], start=fs[1], end=fs[2], sv_type=f(fs[3]), sv_len=fs[4], seq=g(fs), fmt='1/1:1:1')
                except Exception as err:
                    msg = '{}\nError reading line {} of BED input:\n{}\n{!r}'.format(
                            traceback.format_exc(),
                            lineno, line, err)
                    raise Exception(msg)
