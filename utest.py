import os.path as op
from nose.tools import assert_equal, assert_true, assert_false
from collections import defaultdict
from svkits.add_an_indel_to_fasta import *
from svkits.add_indels_to_fasta import *
from svkits.utils import *

OUT_DIR = op.join(op.dirname(__file__), "tests/out/")


def test_weighted_choose():
    r =  weighted_choose(100, {'a':0.5, 'b':0.5})
    assert_equal(sorted(r.keys()), ['a', 'b'])
    assert_equal(r['a'] + r['b'], 100)

def test_weighted_choose_multiple():
    r = weighted_choose_multiple({10: 100, 20: 200}, {'a':0.4, 'b':0.6})
    assert_equal(sorted(r.keys()), [10, 20])
    assert_equal(r[10]['a'] + r[10]['b'], 100)
    assert_equal(r[20]['a'] + r[20]['b'], 200)

def test_swap_dict_dict_k1_k2():
    dd = {10: {'chr1': 5, 'chr2': 6}, 20: {'chr1': 7, 'chr2': 8}}
    r = swap_dict_dict_k1_k2(dd)
    assert_equal(r['chr1'][10], 5)
    assert_equal(r['chr2'][10], 6)
    assert_equal(r['chr1'][20], 7)
    assert_equal(r['chr2'][20], 8)

def test_get_ins_poses():
    poses = get_ins_poses(5000, [100, 200], [2, 5])
    for index in range(1, len(poses)):
        #if not (poses[index] > poses[index-1] + 200):
        #    raise ValueError("poses[index]=%r not > poses[index-1] + 200 = %r + 200" % (poses[index], poses[index-1]))
        assert_true(poses[index] < 5000)

def test_expand_objects():
    assert_equal(expand_objects(['a', 'b'], [2,3]), ['a','a', 'b','b','b'])

def test_filter_dict_by_keys():
    d = {'chr1':'G', 'chr2':'A', 'rand':'C'}
    d1 = filter_dict_by_keys(d, '.+')
    assert_equal(d, d1)

    d2 = filter_dict_by_keys(d, '^chr')
    assert_equal(d2['chr1'], 'G')
    assert_equal(d2['chr2'], 'A')
    assert_false('rand' in d2)

    d3 = filter_dict_by_keys(d, '^((?!rand).)*$')
    assert_true('chr1' in d3)
    assert_true('chr2' in d3)
    assert_false('rand' in d3)

    d['chr_bla'] = 'A'
    d4 = filter_dict_by_keys(d, '^chr(.|..)$')
    assert_true('chr1' in d4)
    assert_true('chr2' in d4)
    assert_false('chr_bla' in d4)

def test_fofn2fns():
    """Get filenames from fofn"""
    o_fn = op.join(OUT_DIR, 'test_fofn2fns.fofn')
    with open(o_fn, 'w') as f:
        f.write("124\n")
    assert_equal(fofn2fns(o_fn) , ['124'])

    o_fn = op.join(OUT_DIR, 'test_fofn2fns')
    assert_equal(fofn2fns(o_fn) , [o_fn])

def test_get_movie2bams_from_fofn():
    fofn = op.join(OUT_DIR, 'test_get_movie2bams_from_fofn.fofn')
    with open(fofn, 'w') as f:
        f.write('/path/movie1.1.subreads.bam\n/path/movie1.2.subreads.bam\nmovie2.subreads.bam')
    o = get_movie2bams_from_fofn(fofn)
    e = {'movie1': set(['/path/movie1.1.subreads.bam', '/path/movie1.2.subreads.bam']), 'movie2':set(['movie2.subreads.bam'])}
    for k in e.keys():
        assert_equal(o[k], e[k])

def test_inv_a_substr_of_seq():
    o = inv_a_substr_of_seq('ABCAAACCC', 3, 5)
    e = 'ABCCCAAAC'
    assert_equal(o, e)

def test_del_a_substr_from_seq():
    o = del_a_substr_from_seq('ABCAAACCC', 2, 1)
    e = 'ABAAACCC'
    assert_equal(o, e)

def test_add_a_str_to_seq():
    o = add_a_str_to_seq('ABCTTT', 3, 10)
    assert_equal(o[0:3], 'ABC')
    assert_equal(o[13:], 'TTT')
