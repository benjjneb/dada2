################################################################################
#' Construct a sample-by-sequence observation matrix.
#' 
#' This function contructs a sequence table (analogous to an OTU table) from
#' the provided list of samples. 
#' 
#' @param samples (Required). A \code{list} of the samples to include in the sequence table. 
#' Names are propagated to the rownames of the return matrix. Samples can be provided in any
#' format that can be processed by \code{\link{getUniques}}.
#' 
#' @return Integer \code{matrix}.
#' A row for each sample, and a column for each unique sequence across all the samples.
#' Note that the columns are named by the sequence which can make display a little unwieldy.
#' 
#' @seealso \code{\link{dada}}, \code{\link{getUniques}}
#' @export
#' 
#' @examples
#' derep1 <- derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' derep2 <- derepFastq(system.file("extdata", "sam2F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, tperr1)
#' dada2 <- dada(derep2, tperr1)
#' makeSequenceTable(list(sample1=dada1, sample2=dada2))
#' 
makeSequenceTable <- function(samples) {
  if(class(samples) %in% c("dada", "derep")) { samples <- list(samples) }
  if(!is.list(samples)) { stop("Requires a list of samples.") }
  unqs <- lapply(samples, getUniques)
  unqsqs <- unique(do.call(c, lapply(unqs, names)))
  if(length(unique(nchar(unqsqs)))) { message("The sequences being tabled vary in length.") }
  rval <- matrix(0, nrow=length(unqs), ncol=length(unqsqs))
  # Samples are rows, columns are sequences
  colnames(rval) <- unqsqs
  for(i in seq_along(unqs)) {
    rval[i, match(names(unqs[[i]]), colnames(rval))] <- unqs[[i]]
  }
  if(!is.null(names(unqs))) {
    rownames(rval) <- names(unqs)
  }
  return(rval)
}
