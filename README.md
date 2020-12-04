![ONT_logo](/ONT_logo.png)
-----------------------------

spliced_bam2gff - a tool to convert spliced BAM alignments into GFF2 format
===========================================================================
[![bioconda-badge](https://anaconda.org/bioconda/spliced_bam2gff/badges/installer/conda.svg)](https://anaconda.org/bioconda/spliced_bam2gff)

The **spliced_bam2gff** tool converts [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) alignments, produced by spliced aligners (such as [minimap2](https://github.com/lh3/minimap2), [gmap](http://research-pub.gene.com/gmap)), into a [GFF2](https://www.ensembl.org/info/website/upload/gff.html) format.

By default, introns are created based on the `N` cigar feature. Alternatively, if `-d` (i.e. for deletion) is specified, any deletions larger than the limit will be classified as an intron.
The orientation of the GFF2 features is determined by the `XS` strand tag and SAM flags depending on the aligner.

The tool supports splitting the output into *loci* and *bundles of loci* with a minimum number of features, which enables easy parallelisation of downstream analyses.
The generated GFF2 files can be compared to a reference annotation using the [gffcompare](https://github.com/gpertea/gffcompare) tool.

## Installation

The best way to install the tool is from bioconda:

- Make sure you have [miniconda3](https://docs.conda.io/en/latest/miniconda.html) installed.
- Set up the bioconda channel according to the [instructions](://bioconda.github.io/user/install.html).
- Install the tool by issuing `conda install spliced_bam2gff`

## Usage

```
Usage of ./spliced_bam2gff:
  -L string
        Write output partitioned into "loci" to this directory. Turns of output to stdout.
  -M    Input is from minimap2.
  -S    Do NOT discard secondary and supplementary alignments.
  -V    Print out version.
  -b int
        Bundle together loci in batches with at least this size. (default -1)
  -d int
        Classify all deletions larger than this as introns (-1 means off). (default -1)
  -g    Use strand tag as feature orientation then read strand if not available.
  -h    Print out help message.
  -s    Use read strand (from BAM flag) as feature orientation.
  -t int
        Number of cores to use. (default 8)
```

Examples
========

### Convert alignments generated by `minmap2 -x splice`

```
spliced_bam2gff -M ./test_data/sirv_simulated.bam > ./test_data/sirv_simulated_mm2.gff
```

Convert while classifying all deletions larger than 20 as introns:
```
spliced_bam2gff -d 20 -M ./test_data/sirv_simulated.bam > ./test_data/sirv_simulated_mm2.gff
```

Convert to GFF2 and split the output into loci separated by regions with no coverage:
```
spliced_bam2gff -M -L ./test_data/sirv_simulated_mm2_loci -s ./test_data/sirv_simulated.bam
```

Convert to GFF2 and split the output into bundles of loci with at least two thousand features:
```
spliced_bam2gff -b 2000 -M -L ./test_data/sirv_simulated_mm2_loci -s ./test_data/sirv_simulated.bam
```

### Convert alignments generated by `gmap`

```
spliced_bam2gff  -g ./test_data/sirv_errors_gmap.bam > ./test_data/test_out_err.gff
```

Running tests
============

For running tests the following dependencies have to be installed:

- [gffcompare](https://github.com/gpertea/gffcompare)

Which is easy to install using [bioconda](https://bioconda.github.io). 
Look into the `Makefiles` for targets testing the tools on simulated and real data.

Help
====

## Licence and Copyright

(c) 2020 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

## FAQs and tips

- The [GFF2](https://www.ensembl.org/info/website/upload/gff.html) files can be visualised using [IGV](http://software.broadinstitute.org/software/igv).
- The GFF2 files can be converted to GFF3 or GTF using the [gffread](https://bioconda.github.io/recipes/gffread/README.html) utility.

## References and Supporting Information

See the post announcing the tool at the Oxford Nanopore Technologies community [here](https://community.nanoporetech.com/posts/new-transcriptomics-analys).
