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
#' seqtab <- makeSequenceTable(list(sample1=dada1, sample2=dada2))
#' 
makeSequenceTable <- function(samples, orderBy = "abundance") {
  if(class(samples) %in% c("dada", "derep", "data.frame")) { samples <- list(samples) }
  if(!is.list(samples)) { stop("Requires a list of samples.") }
  unqs <- lapply(samples, getUniques)
  unqsqs <- unique(do.call(c, lapply(unqs, names)))
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
#' be thought of as implementing greedy 100\% OTU clustering with end-gapping ignored.
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
#' @param identicalOnly (Optional). \code{logical(1)}. Default FALSE.
#' If TRUE, only identical sequences (i.e. duplicates) are collapsed together.
#' 
#' @param vec (Optional). \code{logical(1)}. Default TRUE.
#' Use the vectorized aligner. Should be turned off if sequences exceed 2kb in length.
#' 
#' @param band 	(Optional). \code{numeric(1)}. Default -1 (no banding). The Needleman-Wunsch 
#' alignment can be banded. This value specifies the radius of that band. Set band = -1 
#' to turn off banding.
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
collapseNoMismatch <- function(seqtab, minOverlap=20, orderBy="abundance", identicalOnly=FALSE, vec=TRUE, band=-1, verbose=FALSE) {
  # Collapse identical sequences (duplicates)
  dupes <- duplicated(colnames(seqtab))
  if(any(dupes)) { # Collapse duplicates first
    st <- seqtab[,!dupes,drop=FALSE] # Deduplicated matrix
    for(i in which(dupes)) {
      sq <- colnames(seqtab)[[i]]
      st[,sq] <- st[,sq] + seqtab[,i]
    }
    seqtab <- st # Use deduplicated sequence table going forward
  }
  if(identicalOnly) { return(seqtab) }
  
  # Collapse sequences with no mismatches
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
        if(nwhamming(query,ref,vec=vec,band=band) == 0) {
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
  collapsed <- collapsed[,colnames(collapsed) %in% seqs.out,drop=FALSE]
  
  # Order columns
  if(!is.null(orderBy)) {
    if(orderBy == "abundance") {
      collapsed <- collapsed[,order(colSums(collapsed), decreasing=TRUE),drop=FALSE]
    } else if(orderBy == "nsamples") {
      collapsed <- collapsed[,order(colSums(collapsed>0), decreasing=TRUE),drop=FALSE]
    }
  }
  
  collapsed <- collapsed[,order(colSums(collapsed), decreasing=TRUE),drop=FALSE]
  
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
  
  # Order by decreasing abundance
  ord <- order(derepCounts, decreasing=TRUE)
  derepCounts <- derepCounts[ord]
  derepQuals <- derepQuals[ord,]
  derepMap <- match(seq(length(ord)), ord)[derepMap]
  
  rval <- list(uniques=derepCounts, quals=derepQuals, map=derepMap)
  rval <- as(rval, "derep")
  rval
}

# Check the validity of a putatitve sequence table object, provide useful diagnostic messages if invalid.
is.sequence.table <- function(tab, verbose=TRUE) {
  if(!(is.matrix(tab))) {
    if(verbose) warning("Not a matrix.")
    return(FALSE)
  }
  if(!(all(tab>=0))) {
    if(verbose) warning("Negative entries.")
    return(FALSE)
  }
  if(is.null(colnames(tab))) {
    if(verbose) warning("Column (sequence) names are missing.")
    return(FALSE)
  }
  if(is.null(rownames(tab))) {
    if(verbose) warning("Row (sample) names are missing.")
    return(FALSE)
  }
  if(!all(sapply(colnames(tab), nchar)>0)) {
    if(verbose) warning("Some column (sequence) names are blank.")
    return(FALSE)
  }
  if(!all(sapply(rownames(tab), nchar)>0)) {
    if(verbose) warning("Some row (sample) names are blank.")
    return(FALSE)
  }
  if(any(duplicated(colnames(tab)))) {
    if(verbose) warning("Duplicated column (sequences) names.")
  }
  if(any(duplicated(rownames(tab)))) {
    if(verbose) warning("Duplicated row (sample) names.")
  }
  return(TRUE)
}

################################################################################
#' Merge two or more sample-by-sequence observation matrices.
#' 
#' This function combines sequence tables together into one merged sequences table.
#' 
#' @param table1 (Optional, default=NULL). Named integer matrix. Rownames correspond to samples
#' and column names correspond to sequences. The output of \code{\link{makeSequenceTable}}.
#' 
#' @param table2 (Optional, default=NULL). Named integer matrix. Rownames correspond to samples
#' and column names correspond to sequences. The output of \code{\link{makeSequenceTable}}.
#' 
#' @param ... (Optional). Additional sequence tables.
#' 
#' @param tables (Optional, default=NULL). Either a list of sequence tables, or a list/vector of RDS filenames
#' corresponding to sequence tables. If provided, \code{table1}, \code{table2}, and any
#' additional arguments will be ignored.
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
#' @param tryRC (Optional). \code{logical(1)}. Default FALSE.
#' If tryRC=TRUE, sequences whose reverse complement matches an earlier sequence will be reverse-
#' complemented and merged together with that earlier sequence. This is most useful when different
#' runs sequenced the same gene region in different or mixed orientations. Note, this does not
#' guarantee consistent orientatation from e.g. 5' to 3' on the gene, it just ensures that identical
#' sequences in different orientations are merged.
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
#'   mergetab <- mergeSequenceTables(seqtab1, seqtab2, seqtab3) # unnamed arguments are assumed to be individual sequence tables
#'   input_tables <- list(seqtab1, seqtab2, seqtab3)
#'   mergetab <- mergeSequenceTables(tables=input_tables) # list of sequence tables
#'   files <- c(file1, file2, file3)
#'   mergetab <- mergeSequenceTables(tables=files) # vector of filenames
#' }
#' 
mergeSequenceTables <- function( table1=NULL, table2=NULL, ..., tables=NULL, repeats="error", orderBy = "abundance", tryRC=FALSE) {
  # Convert to sequence tables if necessary
  if (is.null(tables) && (is.null(table1) || is.null(table2))) {
  	stop("Either 'tables' or 'table1' and 'table2' must be provided.")
  } else if (!is.null(tables)) {
		if (is.list(tables)) {
			if(all(sapply(tables, is.character))) {
				tables <- lapply(tables, readRDS)
			}
		} else if (is.character(tables)) {
			tables <- sapply(tables, readRDS)
		} else {
			stop("'tables' must be a list or vector.")
		}
	} else {
		tables <- list(table1, table2)
		tables <- c(tables, list(...))
	}
 
  # Validate tables
  if(length(tables)<2){
    stop("At least two sequence tables are expected")
  }
  tablesValid <- sapply(tables, is.sequence.table)
  if(!(all(tablesValid))) {
    errorMessage <- paste0(names(tables[which(!tablesValid)]), collapse=", ")
    if(length(errorMessage)){
      errorMessage <- paste0(": ", errorMessage)
    }
    stop("Some sequence tables found invalid", errorMessage)
  }
  sample.names <- c(sapply(tables, rownames), recursive=TRUE)
  namesDuplicated <- duplicated(sample.names)
  if(any(namesDuplicated)) {
    if(repeats=="error") {
      stop("Duplicated sample names detected in the sequence table row names: ", paste0(unique(sample.names[which(namesDuplicated)]), collapse=", "))
    }
    else { 
      sample.names <- unique(sample.names)
      message("Duplicated sample names detected in the sequence table row names.")
    }
  }
  seqs <- unique(c(sapply(tables, colnames), recursive=TRUE))
  if(tryRC && length(seqs) > 1) {
    earlier.rc <- c(FALSE, sapply(seq(2, length(seqs)), function(i) rc(seqs[[i]]) %in% seqs[1:(i-1)]))
    rc.seqs <- seqs[earlier.rc]
    if(length(rc.seqs) > 0) {
      message("Reverse complemented sequences detected and re-oriented.")
      for(i in seq_along(tables)) {
        do.rc <- colnames(tables[[i]]) %in% rc.seqs
        if(any(do.rc)) {
          colnames(tables[[i]])[do.rc] <- rc(colnames(tables[[i]])[do.rc])
          tables[[i]] <- collapseNoMismatch(tables[[i]], identicalOnly=TRUE)
        }
      }
      seqs <- seqs[!seqs %in% rc.seqs]
    }
  }
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
