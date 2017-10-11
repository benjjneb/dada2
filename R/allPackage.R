#' DADA2 package
#'
#' The dada2 package is centered around the DADA2 algorithm for accurate high-resolution
#' of sample composition from amplicon sequencing data. The DADA2 algorithm is both more
#' sensitive and more specific than commonly used OTU methods, and resolves amplicon
#' sequence variants (ASVs) that differ by as little as one nucleotide.
#' 
#' The dada2 package also provides a full set of tools for taking raw amplicon sequencing
#' data all the way through to a feature table representing sample composition. Provided
#' facilities include:
#' 
#' \itemize{
#'  \item Quality filtering (\code{\link{filterAndTrim}}, \code{\link{fastqFilter}}, \code{\link{fastqPairedFilter}})
#'  \item Dereplication (\code{\link{derepFastq}})
#'  \item Learn error rates (\code{\link{learnErrors}})
#'  \item Sample Inference (\code{\link{dada}})
#'  \item Chimera Removal (\code{\link{removeBimeraDenovo}}, \code{\link{isBimeraDenovo}}, \code{\link{isBimeraDenovoTable}})
#'  \item Merging of Paired Reads (\code{\link{mergePairs}})
#'  \item Taxonomic Classification (\code{\link{assignTaxonomy}}, \code{\link{assignSpecies}})
#' }
#' 
#' @name dada2-package
#' 
#' @author Benjamin Callahan \email{benjamin.j.callahan@@gmail.com}
#' @author Paul J McMurdie II \email{mcmurdie@@stanford.edu}
#' @author Michael Rosen \email{eigenrosen@@gmail.com}
#' @author Susan Holmes \email{susan@@stat.stanford.edu}
#' 
#' @docType package
#' @keywords package
NA