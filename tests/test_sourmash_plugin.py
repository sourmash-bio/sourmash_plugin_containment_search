"""
Tests for sourmash_plugin_contaiment_search.
"""
import os
import pytest
import csv
from pprint import pprint

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed


def test_run_sourmash(runtmp):
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('', fail_ok=True)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    assert runtmp.last_result.status != 0                    # no args provided, ok ;)


def test_run_sourmash_mgsearch(runtmp):
    # mgsearch command exists!
    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('scripts', 'mgsearch', fail_ok=True)

    print(runtmp.last_result.out)
    print(runtmp.last_result.err)

    assert "error: the following arguments are required: query_genome, metagenomes" in runtmp.last_result.err
    assert runtmp.last_result.status != 0                    # no args provided, ok ;)


def test_0_x_podar(runtmp):
    query = utils.get_test_data('0.sig.zip')
    against = utils.get_test_data('SRR606249.k31.sig.zip')
    runtmp.sourmash('scripts', 'mgsearch', query, against)
    
    out = runtmp.last_result.out
    assert "100.0%    54.2       3.1%     SRR606249" in out


def test_0_x_podar_out(runtmp):
    query = utils.get_test_data('0.sig.zip')
    against = utils.get_test_data('SRR606249.k31.sig.zip')

    runtmp.sourmash('scripts', 'mgsearch', query, against, '-o', 'out.csv')
    csvfp = open(runtmp.output('out.csv'), newline='')
    rows = list(csv.DictReader(csvfp))
    assert len(rows) == 1
    row = rows[0]
    pprint(row)

    assert row['scaled'] == '100000'

    columns = ['intersect_bp',
               'match_filename',
               'match_name',
               'match_md5',
               'query_filename',
               'query_name',
               'query_md5',
               'ksize',
               'moltype',
               'scaled',
               'f_query',
               'f_match',
               'f_match_weighted',
               'sum_weighted_found',
               'average_abund',
               'median_abund',
               'std_abund',
               'query_n_hashes',
               'match_n_hashes',
               'match_n_weighted_hashes',
               'jaccard',
               'genome_containment_ani',
               'match_containment_ani',
               'average_containment_ani',
               'max_containment_ani',
               'potential_false_negative'
               ]

    row_keys = set(row.keys())
    col_keys = set(columns)
    assert row_keys == col_keys

    assert round(float(row['f_match']), 3) == 0.010
    assert round(float(row['f_match_weighted']), 3) == 0.031
    assert round(float(row['f_query']), 3) == 1.0
    assert round(float(row['jaccard']), 3) == 0.010
    assert round(float(row['median_abund']), 3) == 53.0
    assert round(float(row['average_abund']), 3) == 54.190
    assert round(float(row['std_abund']), 3) == 14.844
    assert int(row['sum_weighted_found']) == 2276
    assert int(row['match_n_weighted_hashes']) == 73489
    assert int(row['intersect_bp']) == 4200000
    assert int(row['match_n_hashes']) == 4200
    assert int(row['query_n_hashes']) == 42


def test_1_x_podar(runtmp):
    query = utils.get_test_data('1.sig.zip')
    against = utils.get_test_data('SRR606249.k31.sig.zip')
    runtmp.sourmash('scripts', 'mgsearch', query, against)
    
    out = runtmp.last_result.out
    print(out)
    assert "100.0%    45.5       0.4%     SRR606249" in out


def test_1_x_0_require_abundance(runtmp):
    # require abundance!
    query = utils.get_test_data('1.sig.zip')
    against = utils.get_test_data('0.sig.zip')

    with pytest.raises(SourmashCommandFailed):
        runtmp.sourmash('scripts', 'mgsearch', query, against,
                        '--require-abundance', fail_ok=True)
    
    err = runtmp.last_result.err
    print(err)
    assert "must have abundance information" in err


def test_1_x_0_no_abund(runtmp):
    # do not require abundance
    query = utils.get_test_data('1.sig.zip')
    against = utils.get_test_data('0.sig.zip')

    runtmp.sourmash('scripts', 'mgsearch', query, against, '-o', 'out.csv')
    csvfp = open(runtmp.output('out.csv'), newline='')
    rows = list(csv.DictReader(csvfp))
    row = rows[0]
    pprint(row)

    assert round(float(row['f_match']), 3) == 0.0
    assert int(row['intersect_bp']) == 0
    assert int(row['match_n_hashes']) == 42
    assert int(row['query_n_hashes']) == 6

    assert not row['f_match_weighted']
    assert not row['median_abund']
    assert not row['average_abund']
    assert not row['std_abund']
    assert not row['sum_weighted_found']
    assert not row['match_n_weighted_hashes']


def test_manysearch_0_x_podar(runtmp):
    query1 = utils.get_test_data('0.sig.zip')
    query2 = utils.get_test_data('1.sig.zip')
    against = utils.get_test_data('SRR606249.k31.sig.zip')
    runtmp.sourmash('scripts', 'mgmanysearch', '--queries', query1, query2,
                    '--against', against)
    
    out = runtmp.last_result.out
    print(out)
    assert "CP001472.1 Aci...  100.0%    54.2       3.1%     SRR606249" in out
    assert "CP001941.1 Aci...  100.0%    45.5       0.4%     SRR606249" in out
