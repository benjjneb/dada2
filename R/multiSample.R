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
#' @param orderBy (Optional). String. Default is "abundance".
#' Specifies how the sequences (columns) of the returned table should be ordered (decreasing).
#' Valid values: "abundance", "nsamples", NULL.
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
makeSequenceTable <- function(samples, orderBy = "abundance") {
  if(class(samples) %in% c("dada", "derep", "data.frame")) { samples <- list(samples) }
  if(!is.list(samples)) { stop("Requires a list of samples.") }
  unqs <- lapply(samples, getUniques)
  unqsqs <- unique(do.call(c, lapply(unqs, names)))
  if(length(unique(nchar(unqsqs)))>1) { message("The sequences being tabled vary in length.") }
  rval <- matrix(0, nrow=length(unqs), ncol=length(unqsqs))
  # Samples are rows, columns are sequences
  colnames(rval) <- unqsqs
  for(i in seq_along(unqs)) {
    rval[i, match(names(unqs[[i]]), colnames(rval))] <- unqs[[i]]
  }
  if(!is.null(names(unqs))) {
    rownames(rval) <- names(unqs)
  }
  # Order columns
  if(orderBy == "abundance") {
    rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
  } else if(orderBy == "nsamples") {
    rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
  }
  
  return(rval)
}

################################################################################
#' Collapse together sequences that are identical up to shifts and/or length.
#' 
#' This function takes as input a sequence table and returns a sequence table in which
#' any sequences that are identical up to shifts or length variation, i.e. that have
#' no mismatches or internal indels when aligned, are collapsed together. The most abundant
#' sequences is chosen as the representative of the collapsed sequences.
#' 
#' This function can be thought of as implementing greedy 100-percent OTU clustering, where end-gapping
#' is ignored. The current implementation relies on full alignemnts and is therefore much slower
#' than necessary, other solutions may be required for large sequence tables.
#' 
#' @param seqtab (Required). A samples by sequence matrix, the return of \code{\link{makeSequenceTable}}.
#' 
#' @param minOverlap (Optional). \code{numeric(1)}. Default is 20.
#' The minimum amount of overlap between sequences required to collapse them together.
#' 
#' @return Integer \code{matrix}.
#' A row for each sample, and a column for each collapsed sequence across all the samples.
#' Note that the columns are named by the sequence which can make display a little unwieldy.
#' Columns are in the same order (modulo the removed columns) as in the input matrix.
#' 
#' @seealso \code{\link{makeSequenceTable}}
#' 
collapseNoMismatch <- function(seqtab, minOverlap=20, verbose=FALSE) {
  unqs.srt <- sort(getUniques(seqtab), decreasing=TRUE)
  seqs <- names(unqs.srt) # The input sequences in order of decreasing total abundance
  seqs.out <- character(0) # The output sequences (after collapsing)
  # collapsed will be the output sequence table  
  collapsed <- matrix(0, nrow=nrow(seqtab), ncol=ncol(seqtab))
  colnames(collapsed) <- colnames(seqtab) # Keep input ordering for output table
  rownames(collapsed) <- rownames(seqtab)
  for(query in seqs) {
    added=FALSE
    prefix <- substr(query, 1, minOverlap)
    suffix <- substr(query, nchar(query)-minOverlap+1,nchar(query))
    for(ref in seqs.out) { # Loop over the reference sequences already added to output
      # Prescreen to see if costly alignment worthwhile, this all should possibly be C-side
      if(grepl(prefix, ref) || grepl(suffix, ref)) { 
        if(nwhamming(query,ref,band=-1) == 0) { # No mismatches/indels, join more abundant sequence
          collapsed[,ref] <- collapsed[,ref] + seqtab[,query] 
          added=TRUE
          break
        }
      }
    } # for(ref in seqs.out)
    if(!added) {
      collapsed[,query] <- seqtab[,query]
      seqs.out <- c(seqs.out, query)
    }
  } # for(query in seqs)
  if(!identical(unname(colSums(collapsed)>0), colnames(collapsed) %in% seqs.out)) {
    stop("Mismatch between output sequences and the collapsed sequence table.")
  }
  collapsed <- collapsed[,colnames(collapsed) %in% seqs.out,drop=FALSE]
  
  if(verbose) message("Output ", ncol(collapsed), " collapsed sequences out of ", ncol(seqtab), " input sequences.")
  collapsed
}

# Combines a list of derep-class objects into one single derep object
combineDereps <- function(dereps) {
  if(class(dereps) == "derep") dereps <- list(dereps)
  if(!all(sapply(dereps, function(x) class(x)=="derep"))) Rcpp::stop("Requires derep-class objects.")
  
  combined <- dereps[[1]]
  combined$quals <- sweep(combined$quals, 1, combined$uniques, "*")
  for(derep in dereps[2:length(dereps)]) {
    seen <- names(derep$uniques) %in% names(combined$uniques)
    sqs.seen <- names(derep$uniques)[seen]
    uniques.seen <- derep$uniques[seen]
    quals.seen <- derep$quals[seen,,drop=FALSE]
    quals.seen <- sweep(quals.seen, 1, uniques.seen, "*")
    uniques.new <- derep$uniques[!seen]
    quals.new <- derep$quals[!seen,,drop=FALSE]
    quals.new <- sweep(quals.new, 1, uniques.new, "*")
    
    combined$uniques[sqs.seen] <- combined$uniques[sqs.seen] + uniques.seen
    combined$uniques <- c(combined$uniques, uniques.new)
    combined$quals[sqs.seen,] <- combined$quals[sqs.seen,,drop=FALSE] + quals.seen
    combined$quals <- rbind(combined$quals, quals.new)
    map <- match(names(derep$uniques), names(combined$uniques))
    combined$map <- c(combined$map, map[derep$map])
  }
  combined$quals <- sweep(combined$quals, 1, combined$uniques, "/")
  combined
}

