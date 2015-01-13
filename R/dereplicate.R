#' Removes duplicate sequences from DNAStringSet object.
#'
#' Borrowed from hiReadsProcessor package in BioC.
#' Given a DNAStringSet object, the function dereplicates reads and 
#' adds counts=X to the definition line to indicate replication. 
#'
#' @param dnaSet DNAStringSet object to dereplicate. 
#'
#' @return DNAStringSet object with names describing frequency of repeat.
#'
#' @seealso \code{\link{replicateReads}}, \code{\link{removeReadsWithNs}}, 
#' \code{\link{findBarcodes}}, \code{\link{splitByBarcode}}
#'
#' @export
#'
#' @examples 
#' dnaSet <- c("CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC", 
#' "CCTGAATCCTGGCAATGTCATCATC", "ATCCTGGCAATGTCATCATCAATGG", 
#' "ATCAGTTGTCAACGGCTAATACGCG", "ATCAATGGCGATTGCCGCGTCTGCA", 
#' "CCGCGTCTGCAATGTGAGGGCCTAA", "GAAGGATGCCAGTTGAAGTTCACAC") 
#' dereplicateReads(dnaSet)
dereplicateReads <- function(dnaSet) {
  if(!is(dnaSet,"DNAStringSet")) {
    dnaSet <- DNAStringSet(dnaSet)
  }
  if(is.null(names(dnaSet))) {
    message("No names attribute found in dnaSet object...", 
            "using artifically generated names")
    names(dnaSet) <- paste("read", 1:length(dnaSet), sep="-")
  }
  dnaSet <- dnaSet[order(dnaSet)]
  counts <- BiocGenerics::table(dnaSet)
  dnaSet <- unique(dnaSet)
  names(dnaSet) <- paste0(names(dnaSet), 
                          "counts=", 
                          as.integer(counts[names(counts)[names(dnaSet)]]))
  return(dnaSet)
}

# Alternative
# Borrowed from ShortRead package in BioC.
## tables
.stringset_tables <- function(x, n=50, ...) {
  if (length(x) == 0) {
    return(list(top=integer(0),
                distribution=data.frame(
                  nOccurrences=integer(0),
                  nReads=integer(0))))
  }
  ## FIXME: two sorts
  srt <- srsort(x)
  r <- srrank(x)
  t <- tabulate(r)
  o <- order(t, decreasing=TRUE)
  ## n most common sequences
  n <- min(n, sum(t!=0))              # remove duplicates
  top <- head(t[o], n)
  names(top) <- as.character(head(srt[o], n))
  ## overall frequency -- equivalent of table(table(sread))
  tt <- tabulate(t)
  nOccurrences <- seq_along(tt)[tt!=0]
  nReads <- tt[tt!=0]
  ## results
  list(top=top,
       distribution=data.frame(
         nOccurrences=nOccurrences,
         nReads=nReads, row.names=NULL))
}

#setMethod(tables, "XStringSet", .stringset_tables)

