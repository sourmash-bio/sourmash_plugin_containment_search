# sourmash_plugin_containment_search: improved containment search for genomes in metagenomes

This plugin adds a command `sourmash scripts mgsearch` that provides
new & nicer output for searching genomes against metagenomes.

## Installation

```
pip install sourmash_plugin_containment_search
```

## Usage

This command:
```
sourmash scripts mgsearch query.sig metagenome.sig [ metagenome2.sig ...] \
    [ -o output.csv ]
```
will search for the query genome `query.sig` in one or more
`metagenome.sig` files, producing decent human-readable output and
(optionally) useful CSV outputs.

For example,
```
sourmash scripts mgsearch ../sourmash/podar-ref/0.fa.sig ../sourmash/SRR606249.trim.k31.sig.gz
```

produces:
```
Loaded query signature: CP001472.1 Acidobacterium capsulatum ATCC 51196, com...

p_genome avg_abund   p_metag   metagenome name
-------- ---------   -------   ---------------
 100.0%    55.4         3.1%   SRR606249
```

## Backstory: Why this command?

`sourmash search` supports sample search x sample search, broadly -
perhaps too
broadly. [And the output formats aren't that helpful.](https://github.com/sourmash-bio/sourmash/issues/2002)

`sourmash prefetch` supports metagenome overlap search against many genomes,
which is the reverse of this use case. Moreover, [prefetch doesn't provided weighted results](https://github.com/sourmash-bio/sourmash/issues/1828) and its output isn't frendly.

`sourmash gather` has friendly and useful output, but calculates something
different from overlap.

There is also some interest in [reverse containment search](https://github.com/sourmash-bio/sourmash/issues/1198).

The `manysearch` command of
[the sourmash branchwater plugin](https://github.com/sourmash-bio/sourmash_plugin_branchwater/blob/main/doc/README.md#running-manysearch)
also does a nice containment search like this plugin, but it doesn't
provide nice human-readable output and it also doesn't provide
weighted results.

## Advanced info: other implementation details

This command is streaming, in the sense that it will load each metagenome,
calculate the match, and then discard the metagenome.

## CSV output

Each row contains:

* `intersect_bp` - overlap between genome and metagenome.
* `match_filename` - metagenome filename from sketch.
* `match_name` - metagenome name.
* `match_md5` - metagenome md5.
* `query_filename` - genome filename from sketch.
* `query_name` - genome name.
* `query_md5` - genome md5.
* `ksize` - ksize of comparison.
* `moltype` - moltype of comparison.
* `scaled` - scaled of comparison.
* `f_query` - fraction of query (genome) found. "Detection"; roughly matches the number of bases that will be covered by mapped metagenome reads.
* `f_match` - fraction of metagenome found, unweighted.
* `f_match_weighted` - fraction of metagenome found, weighted. Roughly matches the fraction of metagenome reads that will map to this genome.
* `sum_weighted_found` - sum of weights from intersectiong hashes.
* `average_abund` - average abundance of weights intersecting hashes.
* `median_abund` - median abundance of weights from intersecting hashes.
* `std_abund` - std dev of weights from intersecting hashes.
* `total_weighted_hashes` - total number of weighted hashes in metagenome.

## TODO

* write tests
* evaluate whether we should add more columns by looking at prefetch and gather output

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`containment_search` is developed at https://github.com/ctb/sourmash_plugin_containment_search.

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
python -m build
```

followed by `twine upload dist/...`.
