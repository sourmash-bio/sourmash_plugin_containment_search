# sourmash_plugin_containment_search: improved containment search for genomes in metagenomes

This plugin provides two commands `sourmash scripts mgsearch` and
`sourmash scripts mgmanysearch`, that provide new & nicer outputs for
searching for genomes in metagenomes.  It is a plugin for [the sourmash
software](https://sourmash.bio/).

## Background

Reporting the presence and estimating the abundance of queries in data
sets is a core requirement of many bioinformatics analyses -
metagenomics in particular.

This plugin provides two commands that use k-mers to estimate the
presence and abundance of queries in data sets. Use cases include:

* detecting specific genomes in metagenomes;
* estimating the abundance of genomes in metagenomes;
* contig-level abundance estimation in metagenomes for binning;
* strain-level analysis of genomes in metagenomes;

The plugin uses
[FracMinHash-based estimation](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2)
to calculate k-mer
[detection](https://anvio.org/vocabulary/#detection) and estimate
[coverage](https://anvio.org/vocabulary/#coverage) based on k-mer
multiplicity. These numbers correspond closely to mapping-based
detection and coverage.

The plugin outputs Average Nucleotide Identity estimates using the
approach
[described in Rahman Hera et al., 2023](https://genome.cshlp.org/content/33/7/1061)
and implemented in sourmash.

## Installation

To install this plugin, run:
```
pip install sourmash_plugin_containment_search
```

(This will install sourmash if you do not already have it installed.)

## Usage

This plugin enables two commands, `mgsearch` and `mgmanysearch`.

### `mgsearch` - search for a single query in many data sets

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

This plugin will work with all the standard sourmash database types, too.

Note that the metagenomes must have been sketched with `-p abund` to
enable the `avg_abund` and `p_metag` columns.

### `mgmanysearch` - search for many queries in many data sets

This command:
```
sourmash scripts mgmanysearch --queries query1.sig [ query2.sig ... ]\
    --against metagenome.sig [ metagenome2.sig ...] \
    [ -o output.csv ]
```
will search for the queries  `query*.sig` in one or more
`metagenome*.sig` files, producing decent human-readable output and
(optionally) useful CSV outputs.

## Backstory: Why this command?

`sourmash search` supports sample search x sample search, broadly -
perhaps too
broadly. [And the output formats aren't that helpful.](https://github.com/sourmash-bio/sourmash/issues/2002)

`sourmash prefetch` supports metagenome overlap search against many
genomes, which is the reverse of this use case. Moreover,
[prefetch doesn't provided weighted results](https://github.com/sourmash-bio/sourmash/issues/1828)
and its output isn't friendly.

`sourmash gather` has friendly and useful output, but can't be used to
calculate the overlap between a single query genome and many subject
metagenomes.

There is also some interest in
[reverse containment search](https://github.com/sourmash-bio/sourmash/issues/1198).

The `manysearch` command of
[the sourmash branchwater plugin](https://github.com/sourmash-bio/sourmash_plugin_branchwater/blob/main/doc/README.md#running-manysearch)
also does a nice containment search like this plugin, but it doesn't
provide nice human-readable output and it also doesn't provide
weighted results. (`manysearch` is, however, much lower memory &
probably a fair bit faster because it's mostly in Rust.)

## Advanced info: implementation details

This command is streaming, in the sense that it will load each
metagenome, calculate the match, and then discard the metagenome.
Hence its memory usage peaks with the largest metagenome, and its max
should be driven by the size of the query + the size of the largest
metagenome.

## CSV output

Each row contains the following information.

### Comparison details

* `intersect_bp` - overlap between genome and metagenome, estimated by multiplying the number of hashes by the scaled factor used.
* `f_query` - fraction of query (genome) found, aka "detection"; roughly matches the number of bases that will be covered by mapped metagenome reads.
* `f_match` - fraction of metagenome found, unweighted.
* `f_match_weighted` - fraction of metagenome found, weighted. Roughly matches the fraction of metagenome reads that will map to this genome.
* `sum_weighted_found` - sum of weights from intersecting hashes.
* `average_abund` - average abundance of weights intersecting hashes.
* `median_abund` - median abundance of weights from intersecting hashes.
* `std_abund` - std dev of weights from intersecting hashes.
* `jaccard` - (unweighted) Jaccard similarity between sketches.
* `genome_containment_ani` - ANI estimated from the genome containment in the metagenome. Use this for genome ANI estimates.
* `match_containment_ani` - ANI estimated from the metagenome containment in the genome.
* `average_containment_ani` - ANI estimated from the average of the genome and metagenome containments.
* `max_containment_ani` - ANI estimated from the max containment between genome/metagenome.
* `potential_false_negative` - True if the sketch size(s) were too small to give a reliable ANI estimate. False if ANI estimate is reliable.

### Sketch information

* `ksize` - ksize of comparison.
* `moltype` - moltype of comparison.
* `scaled` - scaled of comparison.

### Query (genome) information:

* `query_filename` - genome filename from sketch.
* `query_name` - genome name.
* `query_md5` - genome md5.
* `query_n_hashes` - total number of hashes in the genome.

### Match (metagenome) information:

* `match_filename` - metagenome filename from sketch.
* `match_name` - metagenome name.
* `match_md5` - metagenome md5.
* `match_n_hashes` - total number of hashes in the metagenome.
* `match_n_weighted_hashes` - total number of weighted hashes in metagenome.

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
