![ONT_logo](/ONT_logo.png)

-----------------------------

spliced_bam2gff - a tool to convert spliced BAM alignments into GFF2 format
===========================================================================


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
