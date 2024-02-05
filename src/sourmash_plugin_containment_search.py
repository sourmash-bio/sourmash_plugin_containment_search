"""\
Search for genomes in metagenomes, with improved output.

Based on script 'calc-weighted-overlap' from sourmash-slainte.

CTB, 2/4/2024
"""

usage="""
   sourmash scripts mgsearch <genome> <metagenome> [ <metagenomes> ... ]
"""

epilog="""
See https://github.com/xyz for more examples.

Need help? Have questions? Ask at http://github.com/sourmash/issues!
"""

import argparse
import sourmash
import numpy as np
import shutil
import csv

from sourmash import sourmash_args
from sourmash.search import PrefetchResult
from sourmash.cli.utils import (add_ksize_arg, add_moltype_args,
                                add_scaled_arg)
from sourmash.plugins import CommandLinePlugin
from sourmash.logging import notify, error

###

def _get_screen_width():
    # default fallback is 80x24
    (col, rows) = shutil.get_terminal_size()

    return col

#
# CLI plugin - supports 'sourmash scripts mgsearch'
#

class Command_ContainmentSearch(CommandLinePlugin):
    command = 'mgsearch'             # 'scripts <command>'
    description = __doc__       # output with -h
    usage = usage               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, subparser):
        super().__init__(subparser)
        subparser.add_argument('query_genome',
                               help='sketch to look for')
        subparser.add_argument('metagenomes', nargs='+',
                               help='metagenomes to search')
        subparser.add_argument('-o', '--output', default=None,
                               help='output CSV')

        add_ksize_arg(subparser, default=31)
        add_moltype_args(subparser)
        add_scaled_arg(subparser, default=1000)

    def main(self, args):
        super().main(args)

        moltype = sourmash_args.calculate_moltype(args)
        return mgsearch(args.query_genome, args.metagenomes,
                        ksize=args.ksize,
                        moltype=moltype,
                        scaled=args.scaled,
                        output=args.output)


def mgsearch(query_filename, against_list, *,
             ksize=31, moltype='DNA', scaled=1000, output=None):
    screen_width = _get_screen_width()

    query_ss = sourmash.load_file_as_index(query_filename)
    query_ss = query_ss.select(ksize=ksize, moltype=moltype)
    query_ss = list(query_ss.signatures())
    if len(query_ss) > 1:
        error(f"ERROR: can only have one query; {len(query_ss)} found.")
        return -1

    query_ss = query_ss[0]
    print(f"Loaded query signature: {query_ss._display_name(screen_width - 25)}")

    query_mh = query_ss.minhash
    if scaled:
        query_mh = query_mh.downsample(scaled=scaled)
        query_ss.minhash = query_mh

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

    if output:
        out_fp = open(output, 'w', newline='')
        out_w = csv.DictWriter(out_fp, fieldnames=columns)
        out_w.writeheader()
    else:
        out_fp = None
        out_w = None

    # go through metagenomes one by one

    # display stuff
    first = True
    
    for metag_filename in against_list:
        metag = sourmash.load_file_as_signatures(metag_filename,
                                                 ksize=ksize)
        metag = list(metag)
        assert len(metag) == 1
        metag = metag[0]

        # make sure metag has abundance!
        assert metag.minhash.track_abundance, f"'{metag_filename}' must have abundance information"

        result = PrefetchResult(query_ss, metag, threshold_bp=0,
                                 estimate_ani_ci=False)

        # calculate total weighted hashes for use in denominator:
        total_sum_abunds = metag.minhash.sum_abundances

        # get a flattened copy for use in intersections...
        flat_metag = metag.minhash.flatten()

        # other info!

        results_template = dict(match_md5=metag.md5sum(),
                                match_name=metag.name,
                                match_filename=metag_filename,
                                ksize=ksize,
                                moltype=metag.minhash.moltype,
                                scaled=metag.minhash.scaled)

        # now, get weighted containment for query genome
        intersect_mh = query_mh.intersection(flat_metag)
        w_intersect_mh = intersect_mh.inflate(metag.minhash)

        abunds = list(w_intersect_mh.hashes.values())

        if abunds:
            mean = np.mean(abunds)
            median = np.median(abunds)
            std = np.std(abunds)
            overlap_sum_abunds = w_intersect_mh.sum_abundances
            f_sum_abunds = overlap_sum_abunds / total_sum_abunds
        else:
            mean = median = std = 0
            overlap_sum_abunds = 0
            f_sum_abunds = 0

        results_d = dict(intersect_bp=result.intersect_bp,
                         query_filename=query_ss.filename,
                         query_name=query_ss.name,
                         query_md5=query_ss.md5sum(),
                         f_query=result.f_match_query,
                         f_match=result.f_query_match,
                         f_match_weighted=f_sum_abunds,
                         sum_weighted_found=overlap_sum_abunds,
                         query_n_hashes=len(query_mh),
                         match_n_hashes=len(flat_metag),
                         match_n_weighted_hashes=total_sum_abunds,
                         average_abund=mean,
                         median_abund=median,
                         std_abund=std,
                         jaccard=result.jaccard,
                         genome_containment_ani=result.query_containment_ani,
                         match_containment_ani=result.match_containment_ani,
                         average_containment_ani=result.average_containment_ani,
                         max_containment_ani=result.max_containment_ani,
                         potential_false_negative=result.potential_false_negative,
                         )

        results_d.update(results_template)

        if out_w:
            out_w.writerow(results_d)

        if first:
            print("")
            print("p_genome avg_abund   p_metag   metagenome name")
            print("-------- ---------   -------   ---------------")
            first = False

        name = metag._display_name(screen_width - 41)

        f_genome_found = results_d['f_query']
        pct_genome = f"{f_genome_found*100:.1f}"

        f_metag_weighted = results_d['f_match_weighted']
        pct_metag = f"{f_metag_weighted*100:.1f}"

        avg_abund = results_d['average_abund']
        avg_abund = f"{avg_abund:.1f}"

        print(f'{pct_genome:>6}%  {avg_abund:>6}     {pct_metag:>7}%   {name}')

    if out_fp:
        out_fp.close()
