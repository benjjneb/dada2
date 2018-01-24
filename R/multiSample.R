################################################################################
#' Construct a sample-by-sequence observation matrix.
#' 
#' This function constructs a sequence table (analogous to an OTU table) from
#' the provided list of samples.
#' 
#' @param samples (Required). A \code{list} of the samples to include in the sequence table. 
#' Samples can be provided in any format that can be processed by \code{\link{getUniques}}.
#' Sample names are propagated to the rownames of the sequence table.
#' 
#' @param orderBy (Optional). \code{character(1)}. Default "abundance".
#' Specifies how the sequences (columns) of the returned table should be ordered (decreasing).
#' Valid values: "abundance", "nsamples", NULL.
#' 
#' @return Named integer matrix.
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
  rval <- matrix(0L, nrow=length(unqs), ncol=length(unqsqs))
  # Samples are rows, columns are sequences
  colnames(rval) <- unqsqs
  for(i in seq_along(unqs)) {
    rval[i, match(names(unqs[[i]]), colnames(rval))] <- unqs[[i]]
  }
  if(!is.null(names(unqs))) {
    rownames(rval) <- names(unqs)
  }
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
  }
  
  return(rval)
}

################################################################################
#' Combine together sequences that are identical up to shifts and/or length.
#' 
#' This function takes as input a sequence table and returns a sequence table in which
#' any sequences that are identical up to shifts or length variation, i.e. that have
#' no mismatches or internal indels when aligned, are collapsed together. The most abundant
#' sequence is chosen as the representative of the collapsed sequences. This function can
#' be thought of as implementing greedy 100\% OTU clustering, with end-gapping is ignored.
#' 
#' @param seqtab (Required). A sample by sequence matrix, the return of \code{\link{makeSequenceTable}}.
#' 
#' @param minOverlap (Optional). \code{numeric(1)}. Default 20.
#' The minimum amount of overlap between sequences required to collapse them together.
#' 
#' @param orderBy (Optional). \code{character(1)}. Default "abundance".
#' Specifies how the sequences (columns) of the returned table should be ordered (decreasing).
#' Valid values: "abundance", "nsamples", NULL.
#' 
#' @param vec (Optional). \code{logical(1)}. Default TRUE.
#' Use the vectorized aligner. Should be turned off if sequences exceed 2kb in length.
#' 
#' @param verbose (Optional). \code{logical(1)}. Default FALSE.
#' If TRUE, a summary of the function results are printed to standard output.
#' 
#' @return Named integer matrix.
#' A row for each sample, and a column for each collapsed sequence across all the samples.
#' Note that the columns are named by the sequence which can make display a little unwieldy.
#' Columns are in the same order (modulo the removed columns) as in the input matrix.
#' 
#' @seealso \code{\link{makeSequenceTable}}
#' @export
#' 
#' @examples
#' derep1 <- derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' derep2 <- derepFastq(system.file("extdata", "sam2F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, tperr1)
#' dada2 <- dada(derep2, tperr1)
#' seqtab <- makeSequenceTable(list(sample1=dada1, sample2=dada2))
#' collapseNoMismatch(seqtab)
#' 
collapseNoMismatch <- function(seqtab, minOverlap=20, orderBy="abundance", vec=TRUE, verbose=FALSE) {
  unqs.srt <- sort(getUniques(seqtab), decreasing=TRUE)
  seqs <- names(unqs.srt) # The input sequences in order of decreasing total abundance
  seqs.out <- character(0) # The output sequences (after collapsing)
  # collapsed will be the output sequence table  
  collapsed <- matrix(0L, nrow=nrow(seqtab), ncol=ncol(seqtab))
  colnames(collapsed) <- colnames(seqtab) # Keep input ordering for output table
  rownames(collapsed) <- rownames(seqtab)
  for(query in seqs) {
    added=FALSE
    prefix <- substr(query, 1, minOverlap)
    for(ref in seqs.out) { # Loop over the reference sequences already added to output
      prefix.ref <- substr(ref, 1, minOverlap)
      # Prescreen to see if costly alignment worthwhile, this could all be moved C-side
      if(grepl(prefix, ref, fixed=TRUE) || grepl(prefix.ref, query, fixed=TRUE)) { 
        if(nwhamming(query,ref,vec=vec,band=16) == 0) {  # band is arbitrary since need exact match
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
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      collapsed <- collapsed[,order(colSums(collapsed), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      collapsed <- collapsed[,order(colSums(collapsed>0), decreasing=TRUE),drop=FALSE]
    }
  }
  
  collapsed <- collapsed[,order(colSums(collapsed), decreasing=TRUE)]
  
  if(verbose) message("Output ", ncol(collapsed), " collapsed sequences out of ", ncol(seqtab), " input sequences.")
  collapsed
}

# Combines a list of derep-class objects into one single derep object
#' @importFrom methods as
#' @keywords internal
combineDereps2 <- function(dereps) {
  if(class(dereps) == "derep") dereps <- list(dereps)
  if(!all(sapply(dereps, function(x) class(x)=="derep"))) stop("Requires derep-class objects.")
  maxlen <- max(sapply(dereps, function(x) ncol(x$quals)))

  # Generate the unique sequences and make the output $uniques vector
  sqs.all <- unique(do.call(c, lapply(dereps, getSequences)))
  derepCounts <- integer(length=length(sqs.all))
  names(derepCounts) <- sqs.all
  
  # Make the output $qual matrix with the appropriate size and rownames
  derepQuals <- matrix(0.0, nrow=length(derepCounts), ncol=maxlen)
  rownames(derepQuals) <- sqs.all
  
  # Initialize the $map with appropriate length
  derepMap <- integer(length=sum(sapply(dereps, function(x) length(x$map))))
  
  start.map <- 1
  for(derep in dereps) {
    if(ncol(derep$quals)<maxlen) { derep$quals <- cbind(derep$quals, matrix(NA, nrow=nrow(derep$quals), ncol=(maxlen-ncol(derep$quals)))) }
    derepCounts[names(derep$uniques)] <- derepCounts[names(derep$uniques)] + derep$uniques
    derepQuals[rownames(derep$quals),] <- derepQuals[rownames(derep$quals),] + sweep(derep$quals, 1, derep$uniques, "*")
    map <- match(names(derep$uniques), names(derepCounts))
    derepMap[start.map:(start.map+length(derep$map)-1)] <- map[derep$map]
    start.map <- start.map + length(derep$map)
  }
  
  derepQuals <- sweep(derepQuals, 1, derepCounts, "/")
  
  rval <- list(uniques=derepCounts, quals=derepQuals, map=derepMap)
  rval <- as(rval, "derep")
  rval
}

is.sequence.table <- function(tab) {
  rval <- is.matrix(tab) && all(tab>=0) && 
    !is.null(colnames(tab)) && !is.null(colnames(tab)) &&
    all(sapply(colnames(tab), nchar)>0) &&
    all(sapply(rownames(tab), nchar)>0)
  rval
}

################################################################################
#' Merge two or more sample-by-sequence observation matrices.
#' 
#' This function combines sequence tables together into one merged sequences table.
#' 
#' @param table1 (Required). Named integer matrix. Rownames correspond to samples
#' and column names correspond to sequences. The output of \code{\link{makeSequenceTable}}.
#' 
#' @param table2 (Required). Named integer matrix. Rownames correspond to samples
#' and column names correspond to sequences. The output of \code{\link{makeSequenceTable}}.
#' 
#' @param ... (Optional). Additional sequence tables.
#' 
#' @param repeats (Optional). Default "error".
#'  Specifies how merging should proceed in the presence of repeated sample names.
#'  Valid values: "error", "sum".
#'  If "sum", then samples with the same name are summed together in the merged table.
#' 
#' @param orderBy (Optional). \code{character(1)}. Default "abundance".
#' Specifies how the sequences (columns) of the returned table should be ordered (decreasing).
#' Valid values: "abundance", "nsamples", NULL.
#' 
#' @return Named integer matrix.
#' A row for each sample, and a column for each unique sequence across all the samples.
#' Note that the columns are named by the sequence which can make display unwieldy.
#' 
#' @seealso \code{\link{makeSequenceTable}}
#' @export
#' 
#' @examples
#' 
#' \dontrun{
#'   mergetab <- mergeSequenceTables(seqtab1, seqtab2, seqtab3)
#' }
#' 
mergeSequenceTables <- function(table1, table2, ..., repeats="error", orderBy = "abundance") {
  # Combine passed tables into a list
  tables <- list(table1, table2)
  tables <- c(tables, list(...))
  # Validate tables
  if(!(all(sapply(tables, is.sequence.table)))) {
    stop("At least two valid sequence tables, and no invalid objects, are expected.")
  }
  sample.names <- c(sapply(tables, rownames), recursive=TRUE)
  if(any(duplicated(sample.names))) {
    if(repeats=="error") { stop("Duplicated sample names detected in the rownames.") }
    else { 
      sample.names <- unique(sample.names)
      message("Duplicated sample names detected in the rownames.") 
    }
  }
  seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
  # Make merged table
  rval <- matrix(0L, nrow=length(sample.names), ncol=length(seqs))
  rownames(rval) <- sample.names
  colnames(rval) <- seqs
  for(tab in tables) {
    rval[rownames(tab), colnames(tab)] <- rval[rownames(tab), colnames(tab)] + tab
  }
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      rval <- rval[,order(colSums(rval), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      rval <- rval[,order(colSums(rval>0), decreasing=TRUE),drop=FALSE]
    }
  }
  rval
}
