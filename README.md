
[![Build Status](https://app.travis-ci.com/benjjneb/dada2.svg?branch=master)](https://app.travis-ci.com/benjjneb/dada2)

# dada2

Exact sample inference from high-throughput amplicon data. Resolves real variants differing by as little as one nucleotide. Visit [the DADA2 website](https://benjjneb.github.io/dada2/index.html) for the most detailed and up-to-date documentation.

### Installation

The dada2 package binaries are available through Bioconductor:

```S
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
```

In order to install dada2 from source (and get the latest and greatest new features) see our [installation from source instructions](https://benjjneb.github.io/dada2/dada-installation.html).

### Documentation

The [tutorial walkthrough of the DADA2 pipeline on paired end Illumina Miseq data](https://benjjneb.github.io/dada2/tutorial.html). 

The [dada2 R package manual](https://www.bioconductor.org/packages/3.6/bioc/manuals/dada2/man/dada2.pdf).

Further documentation is available on [the DADA2 front page](http://benjjneb.github.io/dada2/). 

### DADA2 Articles

[DADA2: High resolution sample inference from Illumina amplicon data. Nature Methods, 2016.](http://dx.doi.org/10.1038/nmeth.3869) [(Open Access link.)](http://rdcu.be/ipGh)

[Bioconductor workflow for microbiome data analysis: from raw reads to community analyses. F1000 Research, 2016.](https://f1000research.com/articles/5-1492)

[Exact sequence variants should replace operational taxonomic units in marker-gene data analysis. ISMEJ, 2017.](http://dx.doi.org/10.1038/ismej.2017.119)

[High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution. Nucleic Acids Research, 2019.](http://dx.doi.org/10.1093/nar/gkz569)

### Other Resources

Planned feature improvements are publicly catalogued at the main DADA2 development site on github, specifically on the "Issues" page for DADA2:

https://github.com/benjjneb/dada2/issues

If the feature you are hoping for is not listed, you are welcome to add it as a feature request "issue" on this page. This request will be publicly available and listed on the page.

Bugs and difficulties in using DADA2 are also welcome on [the issue tracker](https://github.com/benjjneb/dada2/issues).
