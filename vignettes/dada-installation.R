## ----bioc-install, eval=FALSE--------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(suppressUpdates = FALSE)
#  biocLite("ShortRead", suppressUpdates = FALSE)

## ----install-packages, eval=FALSE----------------------------------------
#  install.packages("path/to/dada2",
#                   repos = NULL,
#                   type = "source",
#                   dependencies = c("Depends", "Suggests","Imports"))

## ----install-github-example, eval=FALSE----------------------------------
#  install.packages("~/github/dada2",
#                   repos = NULL,
#                   type = "source",
#                   dependencies = c("Depends", "Suggests","Imports"))

## ----packageVersion------------------------------------------------------
packageVersion("dada2")

## ----bioc-install-missing, eval=FALSE------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("missing_package_1")
#  biocLite("missing_package_2")
#  # ... and so on

## ----install-packages-rev2, eval=FALSE-----------------------------------
#  install.packages("path/to/dada2",
#                   repos = NULL,
#                   type = "source",
#                   dependencies = c("Depends", "Suggests","Imports"))

## ----load-dada2, message=FALSE-------------------------------------------
library("dada2")

## ----documentation-example, eval=FALSE, message=FALSE--------------------
#  ?dada
#  ?derepFastq
#  ?importUniques

## ----test-dada-0, eval=FALSE---------------------------------------------
#  score <- matrix(c(5, -4, -4, -4, -4, 5, -4,
#                    -4, -4, -4, 5, -4, -4, -4, -4, 5),
#                  nrow=4, byrow=T)
#  gap_penalty <- -8
#  err_init <- matrix(c(0.991, 0.003, 0.003, 0.003,
#                       0.003, 0.991, 0.003, 0.003, 0.003,
#                       0.003, 0.991, 0.003, 0.003, 0.003,
#                       0.003, 0.991),
#                     nrow=4, byrow=T)
#  test_fastq_file = system.file("extdata", "test-nonunique.fastq.gz", package="dada2")
#  uniques <- derepFastq(fl = test_fastq_file)
#  uniques.dada <- dada(uniques = uniques$uniques, quals = uniques$quals)
#  uniques.dada.consistent <- dada(uniques, err_init, score, gap_penalty, self_consist = TRUE)
#  length(uniques)
#  length(uniques.dada$genotypes)
#  length(uniques.dada.consistent$genotypes)

## ------------------------------------------------------------------------
# Out of `r length(uniques)` original unique sequences,
# did you get `r length(uniques.dada$genotypes)`
# genotypes in `uniques.dada$genotypes` and
# `r length(uniques.dada.consistent$genotypes)` in `uniques.dada.consistent`?

## ----compute-performance-1, eval=FALSE-----------------------------------
#  group_0_25_file = system.file("extdata", "group_0-25.uniques.gz", package="dada2")
#  bu = importUniques(group_0_25_file)
#  length(bu)
#  library("microbenchmark")
#  microbenchmark(dada(bu, err = err_init, ), times=10)

## ----compute-performance-2, eval=FALSE-----------------------------------
#  group_0_30_file = system.file("extdata", "group_0-30.uniques.gz", package="dada2")
#  bu = NULL
#  bu <- importUniques(group_0_30_file)
#  length(bu)
#  microbenchmark(dada(bu, err_init, score, gap_penalty), times=10)

## ----compute-performance-3, eval=FALSE-----------------------------------
#  library("ggplot2")
#  length(bu)
#  # The max number of sequences from `bu` to include in this example
#  N = 1000L
#  # Set random see for repoducibility of example
#  set.seed(711L)
#  # Randomly subsample from `bu`, because `bu`
#  # on its own is large and takes a long time
#  buSub0 = bu[sample(x = length(bu), size = N, replace = FALSE)]
#  # Run and report time elapsed
#  system.time({
#    df <- calibrate_kmers(seqs = names(buSub0),
#                          score = get_dada_opt("SCORE_MATRIX"),
#                          gap =  -8,
#                          max_aligns = 100000)
#  })
#  # Make a summary plot
#  ggplot(df, aes(x=kmer, y=align)) + geom_bin2d()
#  # Som alternative plots
#  # ggplot(df, aes(x=kmer, y=align)) + geom_hex()
#  # ggplot(df, aes(x=kmer, y=align)) + geom_density2d()

