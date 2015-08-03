################################################################################
#' Read and Dereplicate a Fastq file.
#' 
#' This is a custom interface to \code{\link[ShortRead]{FastqStreamer}} 
#' for dereplicating amplicon sequences from a fastq or compressed fastq file,
#' while also controlling peak memory requirement to support large files.
#'
#' @param fl (Required). Character.
#'  The file path to the fastq or fastq.gz file.
#'  Actually, any file format supported by \code{\link[ShortRead]{FastqStreamer}}.
#' 
#' @param n (Optional). A \code{numeric(1)} indicating
#'  the maximum number of records (reads) to parse and dereplicate
#'  at any one time. This controls the peak memory requirement
#'  so that large fastq files are supported.
#'  Defaults is \code{1e6}, one-million reads.
#'  See \code{\link[ShortRead]{FastqStreamer}} for details on this parameter,
#'  which is passed on.
#' 
#' @param verbose (Optional). A \code{logical(1)} indicating
#'  whether to throw any standard R \code{\link{message}}s 
#'  on the intermittent and final status of the dereplication.
#'  Default is \code{FALSE}, no messages.
#'
#' @return Class "derep". A list with three members. 
#' \itemize{
#'  \item{$uniques: Named integer vector. Named by the unique sequence, valued by abundance.}
#'  \item{$quals: Numeric matrix of average quality scores by position for each unique. Uniques are rows, positions are cols.}
#'  \item{$map: Integer vector of length the number of reads, and value the index (in $uniques) of the unique to which that read was assigned.}
#' }
#'
#' @export
#' @import ShortRead 
#'
#' @examples 
#' # Test that chunk-size, `n`, does not affect the result.
#' testFile = system.file("extdata", "test-nonunique.fastq.gz", package="dadac")
#' test1 = derepFastq(testFile, verbose = TRUE)
#' test2 = derepFastq(testFile, 35, TRUE)
#' test3 = derepFastq(testFile, 100, TRUE)
#' all.equal(test1$uniques, test2$uniques[names(test1$uniques)])
#' all.equal(test1$uniques, test3$uniques[names(test1$uniques)])
derepFastq <- function(fl, n = 1e6, verbose = FALSE, qeff=FALSE){
  if(verbose){
    message("Dereplicating sequence entries in Fastq file: ", fl, appendLF = TRUE)
  }
  
  f <- FastqStreamer(fl, n = n)
  suppressWarnings(fq <- yield(f))
  
  out <- qtables2(fq, qeff)
  
  derepCounts <- out$uniques
  derepQuals <- out$cum_quals
  derepMap <- out$map
  while( length(suppressWarnings(fq <- yield(f))) ){
    # A little loop protection
    newniques = alreadySeen = NULL
    # Dot represents one turn inside the chunking loop.
    if(verbose){
      message(".", appendLF = FALSE)
    }
    out <- qtables2(fq, qeff)
    # identify sequences already present in `derepCounts`
    alreadySeen <- names(out$uniques) %in% names(derepCounts)
    # Sum these values, if any
    if(any(alreadySeen)){
      sqnms = names(out$uniques)[alreadySeen]
      derepCounts[sqnms] <- derepCounts[sqnms] + out$uniques[sqnms]
      derepQuals[sqnms,] <- derepQuals[sqnms,] + out$cum_quals[sqnms,]
    }
    # Concatenate the remainder to `derepCounts`, if any
    if(!all(alreadySeen)){
      derepCounts <- c(derepCounts, out$uniques[!alreadySeen])
      derepQuals <- rbind(derepQuals, out$cum_quals[!alreadySeen,,drop=FALSE])
    }
    new2old <- match(names(out$uniques), names(derepCounts)) # map from out$uniques index to derepCounts index
    if(any(is.na(new2old))) warning("Failed to properly extend uniques.")
    derepMap <- c(derepMap, new2old[out$map])
  }
  derepQuals <- derepQuals/derepCounts # Average
  if(qeff) derepQuals <- -10*log10(derepQuals)  # Convert back to effective Q value
  # Sort by decreasing abundance
  ord <- order(derepCounts, decreasing=TRUE)
  derepCounts <- derepCounts[ord]
  derepQuals <- derepQuals[ord,]
  derepMap <- match(derepMap, ord)
  if(verbose){
    message("Encountered ",
            length(derepCounts),
            " unique sequences from ",
            sum(derepCounts),
            " total sequences read.")
  }
  close(f)
  rval <- list(uniques=derepCounts, quals=derepQuals, map=derepMap)
  class(rval) <- "derep"
  return(rval)
}

################################################################################
#  DEPRECATED DEPRECATED DEPRECATED
#' Read and Dereplicate a Fastq file containing multiple samples.
#' 
#' This is a custom interface to \code{\link[ShortRead]{FastqStreamer}} 
#' for dereplicating amplicon sequences from a multi-sample fastq or fastq.gz file.
#' This function relies heavily on the \code{\link[ShortRead]{tables}} method.
#'
#' @param fl (Required). Character.
#'  The file path to the fastq or fastq.gz file.
#' 
#' @param samsubseq (Required). A \code{numeric(2)} specifying the substring
#'  of the id field that will be used to split the fastq reads into samples.
#' 
#' @param n (Optional). A \code{numeric(1)} indicating
#'  the maximum number of records (reads) to parse and dereplicate
#'  at any one time. This controls the peak memory requirement.
#'  Defaults is \code{1e6}, one-million reads.
#'  See \code{\link[ShortRead]{FastqStreamer}} for details on this parameter,
#'  which is passed on.
#' 
#' @param verbose (Optional). A \code{logical(1)} indicating
#'  whether to throw any standard R \code{\link{message}}s 
#'  on the intermittent and final status of the dereplication.
#'  Default is \code{FALSE}, no messages.
#'
#' @return Named integer vector. Named by sequence, valued by number of occurence.
#'
#' @import ShortRead 
#' 
derepMultiSampleFastq <- function(fl, samsubseq, n = 1e6, verbose = FALSE){
  if(verbose){
    message("Dereplicating multi-sample sequences in Fastq file: ", fl, appendLF = TRUE)
    message("Identifying the sample IDs---", appendLF = TRUE)
  }
  
  ## Initial iteration through file to identify the unique sample IDs
  f <- FastqStreamer(fl, n = n)
  totReads <- 0
  # `yield` is the method for returning the next chunk from the stream.
  # Use `fq` as the "current chunk"
  suppressWarnings(fq <- yield(f))
  totReads <- totReads + length(fq)
  # Extract the sample ids
  samids <- subseq(fq@id, samsubseq[[1]], samsubseq[[2]])
  # Find all unique samples in the first chunk
  samids = unique(samids)
  # This loop will stop if/when the end of the file is already reached.
  # If end is already reached, this loop is skipped altogether.
  while( length(suppressWarnings(fq <- yield(f))) ){
    totReads <- totReads + length(fq)
    # Extract the sample ids
    newsamids <- subseq(fq@id, samsubseq[[1]], samsubseq[[2]])
    # Dot represents one turn inside the chunking loop.
    if(verbose){ message(".", appendLF = FALSE) }

    samids <- unique(c(samids, unique(newsamids)))
  }
  if(verbose){
    message("Encountered ",
            length(samids),
            " unique sample IDs from ",
            totReads,
            " total sequences read.")
  }
  close(f)
  samids <- as(samids, "character")
  if(verbose){
    message("Dereplicating the sequences by sample---", appendLF = TRUE)
  }

  # Initialize output list
  derepSamples <- vector("list", length(samids))
  names(derepSamples) <- samids
  ## Second iteration through file to dereplicate 
  f <- FastqStreamer(fl, n = n)
  # `yield` is the method for returning the next chunk from the stream.
  # Use `fq` as the "current chunk"
  suppressWarnings(fq <- yield(f))
  # Calculate the dereplicated counts for the first chunk
  for(sam in samids) {
    derepSamples[[sam]] = tables(fq[subseq(fq@id,1,10) == sam], n = Inf)$top  # WHAT HAPPENS WHEN YOU TABLES ON A ZERO LENGTH XSTRINGSET?
  }
  # This loop will stop if/when the end of the file is already reached.
  # If end is already reached, this loop is skipped altogether.
  while( length(suppressWarnings(fq <- yield(f))) ){
    # A little loop protection
    idrep = alreadySeen = NULL
    # Dot represents one turn inside the chunking loop.
    if(verbose){
      message(".", appendLF = FALSE)
    }
    for(sam in samids) {
      idrep <- tables(fq[subseq(fq@id,1,10) == sam], n = Inf)$top
      # identify sequences already present in `derepCounts`
      alreadySeen = names(idrep) %in% names(derepSamples[[sam]])
      # Sum these values, if any
      if(any(alreadySeen)){
        sqnms = names(idrep)[alreadySeen]
        derepSamples[[sam]][sqnms] <- derepCounts[sqnms] + idrep[sqnms]
      }
      # Concatenate the remainder to `derepCounts`, if any
      if(!all(alreadySeen)){
        derepSamples[[sam]] <- c(derepSamples[[sam]], idrep[!alreadySeen])
      }
    }
  }
  if(verbose){
    message("Encountered ??",
            " unique sequences in ",
            length(samids),
            " samples from ",
            totReads,
            " total sequences read.")
  }
  close(f)
  return(derepSamples)
}
###
#' Internal function to replicate tables functionality while also returning average quals and a map
#' from reads to uniques
#' 
qtables2 <- function(x, qeff = FALSE) {
  # ranks are lexical rank
  srt <- srsort(x) # map from rank to sequence/quality/id
  rnk <- srrank(x) # map from read_i to rank (integer vec of ranks, ties take top rank)
  tab <- tabulate(rnk) # map from rank to abundance (much faster than table)
  
  srtrnk <- srrank(srt)
  unq <- unique(srtrnk)
  uniques <- tab[unq] # abundance of each rank (value)
  names(uniques) <- as.character(sread(srt)[unq]) # sequence of each unique (name)
  
  rnk2unqi <- rep(seq(length(uniques)), tab[tab>0]) # map from rank to uniques index
  map <- rnk2unqi[rnk] # map from read index to unique index
  
  # do matrices
  qmat <- as(quality(srt), "matrix") # map from read_i to quality
  if(qeff) qmat <- 10^(-qmat/10)  # Convert to nominal error probability first
  qmat <- rowsum(qmat, srtrnk, reorder=FALSE)
  rownames(qmat) <- names(uniques)
  list(uniques=uniques, cum_quals=qmat, map=map)
}

##########
#' Write a uniques vector to a FASTA file
#' 
#' Basically a wrapper for writeFastq in the ShortRead package.
#' Format is suitable for uchime.
#' 
#' @param unqs The uniques object. E.g. the output of \code{\link{derepFastq}}.
#' 
#' @param fout The file path of the output file. The file you want to write.
#' 
#' @param ids A character vector of sequence ids, one for each element in \code{unqs}.
#'  Default value is \code{NULL}, in which case an arbitrary ID is assigned.
#'  
#' @param mode A character string flag passed on to \code{\link[ShortRead]{writeFasta}}
#'  indicating the type of file writing mode. Default is \code{"w"}.
#'  
#' @param width The number of characters per line in the file. Default is 20000.
#'  Passed on to \code{\link[ShortRead]{writeFasta}}.
#'  
#' @param ... Additional parameters passed on to \code{\link[ShortRead]{writeFasta}}.
#' 
#' @import ShortRead
#' @export
#' 
uniquesToFasta <- function(unqs, fout, ids=NULL, mode="w", width=20000, ...) {
  unqs <- as.uniques(unqs)
  if(is.null(ids)) {
    ids <- paste0("sq", seq(1, length(unqs)), ";size=", unname(unqs), ";")
  }
  writeFasta(object = ShortRead(sread = DNAStringSet(names(unqs)),
                                id = BStringSet(ids)),
             file = fout, 
             mode = mode, 
             width = width,
             ...)
}
