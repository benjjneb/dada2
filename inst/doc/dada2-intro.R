## ----filenames, message=FALSE, warning=FALSE-----------------------------
library(dada2); packageVersion("dada2")
fnF1 <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
fnR1 <- system.file("extdata", "sam1R.fastq.gz", package="dada2")

filtF1 <- tempfile(fileext=".fastq.gz")
filtR1 <- tempfile(fileext=".fastq.gz")

## ----filter--------------------------------------------------------------
fastqPairedFilter(c(fnF1, fnR1), fout=c(filtF1, filtR1), maxN=0, trimLeft=10, truncLen=c(240, 200), maxEE=2, compress=TRUE, verbose=TRUE)

## ----derep---------------------------------------------------------------
derepF1 <- derepFastq(filtF1, verbose=TRUE)
derepR1 <- derepFastq(filtR1, verbose=TRUE)

## ----dada, warning=FALSE-------------------------------------------------
dadaF1 <- dada(derepF1, err=inflateErr(tperr1,3), selfConsist=TRUE)
dadaR1 <- dada(derepR1, err=inflateErr(tperr1,3), selfConsist=TRUE)
print(dadaF1)

## ----bimeras-------------------------------------------------------------
bimerasF1 <- isBimeraDenovo(dadaF1, verbose=TRUE)
bimerasR1 <- isBimeraDenovo(dadaR1, verbose=TRUE)

## ----merge---------------------------------------------------------------
merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose=TRUE)
merger1.nochim <- merger1[!bimerasF1[merger1$forward] & !bimerasR1[merger1$reverse],]

## ----sample2, warning=FALSE----------------------------------------------
# Assign filenames
fnF2 <- system.file("extdata", "sam2F.fastq.gz", package="dada2")
fnR2 <- system.file("extdata", "sam2R.fastq.gz", package="dada2")
filtF2 <- tempfile(fileext=".fastq.gz")
filtR2 <- tempfile(fileext=".fastq.gz")
# Filter and Trim
fastqPairedFilter(c(fnF2, fnR2), fout=c(filtF2, filtR2), maxN=0, trimLeft=10, truncLen=c(240, 200), maxEE=2, compress=TRUE, verbose=TRUE)
# Dereplicate
derepF2 <- derepFastq(filtF2, verbose=TRUE)
derepR2 <- derepFastq(filtR2, verbose=TRUE)
# Infer sample composition
dadaF2 <- dada(derepF2, err=inflateErr(tperr1,3), selfConsist=TRUE)
dadaR2 <- dada(derepR2, err=inflateErr(tperr1,3), selfConsist=TRUE)
# Identify chimeras
bimerasF2 <- isBimeraDenovo(dadaF2, verbose=TRUE)
bimerasR2 <- isBimeraDenovo(dadaR2, verbose=TRUE)
# Merge
merger2 <- mergePairs(dadaF2, derepF2, dadaR2, derepR2, verbose=TRUE)
merger2.nochim <- merger2[!bimerasF2[merger2$forward] & !bimerasR2[merger2$reverse],]

## ----make-table----------------------------------------------------------
seqtab <- makeSequenceTable(list(merger1, merger2))

