################################################################################
#' Read and Dereplicate a Fastq file.
#' 
#' This is a custom interface to \code{\link{FastqStreamer}} 
#' for dereplicating amplicon sequences from a fastq or fastq.gz file,
#' while also controlling peak memory requirement to support large files.
#' This function relies heavily on the \code{\link[ShortRead]{tables}} method.
#'
#' @param fl (Required). Character.
#'  The file path to the fastq or fastq.gz file.
#'  Actually, any file format supported by \code{\link{FastqStreamer}}.
#' 
#' @param n (Optional). A \code{numeric(1)} indicating
#'  the maximum number of records (reads) to parse and dereplicate
#'  at any one time. This controls the peak memory requirement
#'  so that large fastq files are supported.
#'  Defaults is \code{1e6}, one-million reads.
#'  See \code{\link{FastqStreamer}} for details on this parameter,
#'  which is passed on.
#' 
#' @param verbose (Optional). A \code{logical(1)} indicating
#'  whether to throw any standard R \code{\link{message}}s 
#'  on the intermittent and final status of the dereplication.
#'  Default is \code{FALSE}, no messages.
#'
#' @return Named integer vector. Named by sequence, valued by number of occurence.
#'
#' @seealso \code{\link{replicateReads}}, \code{\link{removeReadsWithNs}}, 
#' \code{\link{findBarcodes}}, \code{\link{splitByBarcode}}
#'
#' @export
#' @import ShortRead 
#'
#' @examples 
#' # Test that chunk-size, `n`, does not affect the result.
#' testFile = system.file("extdata", "test-nonunique.fastq.gz", package="dadac")
#' test1 = dereplicateFastqReads(testFile, verbose = TRUE)
#' test2 = dereplicateFastqReads(testFile, 35, TRUE)
#' test3 = dereplicateFastqReads(testFile, 100, TRUE)
#' all.equal(test1, test2[names(test1)])
#' all.equal(test1, test3[names(test1)])
derepFastq <- function(fl, n = 1e6, verbose = FALSE){
  # require("ShortRead")
  if(verbose){
    message("Dereplicating sequence entries in Fastq file: ", fl, appendLF = TRUE)
  }
  ## iterating over an entire file using fastq streaming
  f <- FastqStreamer(fl, n = n)
  # `yield` is the method for returning the next chunk from the stream.
  # Use `fq` as the "current chunk"
  suppressWarnings(fq <- yield(f))
  # Calculate the dereplicated counts for the first chunk
  derepCounts = tables(fq, n = Inf)$top
  # This loop will stop if/when the end of the file is already reached.
  # If end is already reached, this loop is skipped altogether.
  while( length(suppressWarnings(fq <- yield(f))) ){
    # A little loop protection
    idrep = alreadySeen = NULL
    # Dot represents one turn inside the chunking loop.
    if(verbose){
      message(".", appendLF = FALSE)
    }
    idrep <- tables(fq, n = Inf)$top
    # identify sequences already present in `derepCounts`
    alreadySeen = names(idrep) %in% names(derepCounts)
    # Sum these values, if any
    if(any(alreadySeen)){
      sqnms = names(idrep)[alreadySeen]
      derepCounts[sqnms] <- derepCounts[sqnms] + idrep[sqnms]
    }
    # Concatenate the remainder to `derepCounts`, if any
    if(!all(alreadySeen)){
      derepCounts <- c(derepCounts, idrep[!alreadySeen])
    }
  }
  if(verbose){
    message("Encountered ",
            length(derepCounts),
            " unique sequences from ",
            sum(derepCounts),
            " total sequences read.")
  }
  close(f)
  return(derepCounts)
}
################################################################################
#' Read and Dereplicate a Fastq file containing multiple samples.
#' 
#' This is a custom interface to \code{\link{FastqStreamer}} 
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
#'  See \code{\link{FastqStreamer}} for details on this parameter,
#'  which is passed on.
#' 
#' @param verbose (Optional). A \code{logical(1)} indicating
#'  whether to throw any standard R \code{\link{message}}s 
#'  on the intermittent and final status of the dereplication.
#'  Default is \code{FALSE}, no messages.
#'
#' @return Named integer vector. Named by sequence, valued by number of occurence.
#'
#' @seealso \code{\link{replicateReads}}, \code{\link{removeReadsWithNs}}, 
#' \code{\link{findBarcodes}}, \code{\link{splitByBarcode}}
#'
#' @export
#' @import ShortRead 
#' 
derepMultiSampleFastq <- function(fl, samsubseq, n = 1e6, verbose = FALSE){
  # require("ShortRead")
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
#' Internal function to replicate tables functionality while also returning average quals
#' 
qtables <- function(x) {
  srt <- srsort(x) # sorts based on sread -- necessary?
  rnk <- srrank(srt) # integer vec of ranks, ties take top rank
  unq <- unique(rnk)
  uniques <- tabulate(rnk)[unq]  # much faster than table
  names(uniques) <- as.character(sread(srt)[unq])
  qmat <- as(quality(srt), "matrix")
  qmat <- rowsum(qmat, rnk)
  rownames(qmat) <- names(uniques)
  list(uniques=uniques, cum_quals=qmat)
}
################################################################################
#' Read and Dereplicate a Fastq, and return average Q score for each unique.
#' 
#' @param fl (Required). Character.
#'  The file path to the fastq or fastq.gz file.
#' 
#' @param n (Optional). A \code{numeric(1)} indicating
#'  the maximum number of records (reads) to parse and dereplicate
#'  at any one time. This controls the peak memory requirement.
#'  Defaults is \code{1e6}, one-million reads.
#'  See \code{\link{FastqStreamer}} for details on this parameter,
#'  which is passed on.
#' 
#' @param verbose (Optional). A \code{logical(1)} indicating
#'  whether to throw any standard R \code{\link{message}}s 
#'  on the intermittent and final status of the dereplication.
#'  Default is \code{FALSE}, no messages.
#'
#' @return A two-element list with `uniques`, an
#'  integer vector named by sequence and valued by number of occurence;
#'  and `quals` the average quality of the dereplicated sequences.
#'
#' @import ShortRead
#'
#' @export
#'
derepFastqWithQual <- function(fl, n = 1e6, verbose = FALSE){
  if(verbose){
    message("Dereplicating sequence entries in Fastq file: ", fl, appendLF = TRUE)
  }

  f <- FastqStreamer(fl, n = n)
  suppressWarnings(fq <- yield(f))
  out <- qtables(fq)
  derepCounts <- out$uniques
  derepQuals <- out$cum_quals
  while( length(suppressWarnings(fq <- yield(f))) ){
    # A little loop protection
    idrep = alreadySeen = NULL
    # Dot represents one turn inside the chunking loop.
    if(verbose){
      message(".", appendLF = FALSE)
    }
    out <- qtables(fq)
    idrep <- out$uniques
    # identify sequences already present in `derepCounts`
    alreadySeen = names(idrep) %in% names(derepCounts)
    # Sum these values, if any
    if(any(alreadySeen)){
      sqnms = names(idrep)[alreadySeen]
      derepCounts[sqnms] <- derepCounts[sqnms] + idrep[sqnms]
      derepQuals[sqnms,] <- derepQuals[sqnms,] + out$cum_quals[sqnms,]
    }
    # Concatenate the remainder to `derepCounts`, if any
    if(!all(alreadySeen)){
      derepCounts <- c(derepCounts, idrep[!alreadySeen])
      derepQuals <- rbind(derepQuals, out$cum_quals[!alreadySeen,])
    }
  }
  derepQuals <- derepQuals/derepCounts # Change to average quals
  if(verbose){
    message("Encountered ",
            length(derepCounts),
            " unique sequences from ",
            sum(derepCounts),
            " total sequences read.")
  }
  close(f)
  return(list(uniques=derepCounts, quals=derepQuals))
}
###
#' Internal function to replicate tables functionality while also returning average quals and a map
#' from reads to uniques
#' 
qtables2 <- function(x) {
  # ranks are lexical rank
  srt <- srsort(x) # map from rank to sequence/quality/id
  rnk <- srrank(x) # map from read_i to rank (integer vec of ranks, ties take top rank)
#  cat("Rnk:", rnk, "\n")
  tab <- tabulate(rnk) # map from rank to abundance (much faster than table)
#  cat("Tab:", tab, "\n")
  
  srtrnk <- srrank(srt)
  #  unq <- unique(rnk)
  unq <- unique(srtrnk)
  uniques <- tab[unq] # abundance of each rank (value)
  names(uniques) <- as.character(sread(srt)[unq]) # sequence of each unique (name)
#  print(uniques)
  
  ##### BUG HERE?
  rnk2unqi <- rep(seq(length(uniques)), tab[tab>0]) # map from rank to uniques index
#  cat("Rnk2unqi:", rnk2unqi, "\n")
  map <- rnk2unqi[rnk] # map from read index to unique index
#  cat("Map:", map, "\n")
  
  # do matrices
  qmat <- as(quality(srt), "matrix") # map from read_i to quality
  qmat <- rowsum(qmat, srtrnk, reorder=FALSE)
  rownames(qmat) <- names(uniques)
  list(uniques=uniques, cum_quals=qmat, map=map)
}
#' @export
#'
derepFastqTest <- function(fl, n = 1e6, verbose = FALSE){
  if(verbose){
    message("Dereplicating sequence entries in Fastq file: ", fl, appendLF = TRUE)
  }
  
  f <- FastqStreamer(fl, n = n)
  suppressWarnings(fq <- yield(f))
  
  ######################## >>>
  out <- qtables2(fq)
  ######################## <<<
  
  derepCounts <- out$uniques
  derepQuals <- out$cum_quals
  derepMap <- out$map
#  print(derepCounts)
#  cat("derepMap:", derepMap, "\n")
  while( length(suppressWarnings(fq <- yield(f))) ){
    # A little loop protection
    newniques = alreadySeen = NULL
    # Dot represents one turn inside the chunking loop.
    if(verbose){
      message(".", appendLF = FALSE)
    }
    ######################## >>>
    out <- qtables2(fq)
    ######################## <<<
    # identify sequences already present in `derepCounts`
    alreadySeen <- names(out$uniques) %in% names(derepCounts)
#    print(alreadySeen)
    # Sum these values, if any
    if(any(alreadySeen)){
#      cat("any\n")
      sqnms = names(out$uniques)[alreadySeen]
      derepCounts[sqnms] <- derepCounts[sqnms] + out$uniques[sqnms]
      derepQuals[sqnms,] <- derepQuals[sqnms,] + out$cum_quals[sqnms,]
#      print(derepCounts)
    }
    # Concatenate the remainder to `derepCounts`, if any
    if(!all(alreadySeen)){
#      cat("!all\n")
      derepCounts <- c(derepCounts, out$uniques[!alreadySeen])
      derepQuals <- rbind(derepQuals, out$cum_quals[!alreadySeen,,drop=FALSE])
#      print(derepCounts)
    }
    ######################## >>>
    new2old <- match(names(out$uniques), names(derepCounts)) # map from out$uniques index to derepCounts index
#  map <- rnk2unqi[rnk] # map from read index to unique index
#    cat("out$map", out$map, "\n")
#    cat("new2old:", new2old, "\n")
    if(any(is.na(new2old))) warning("Failed to properly extend uniques.")
    derepMap <- c(derepMap, new2old[out$map])
#    cat("derepMap:", derepMap, "\n")
    ######################## <<<
  }
  derepQuals <- derepQuals/derepCounts # Change to average quals
  if(verbose){
    message("Encountered ",
            length(derepCounts),
            " unique sequences from ",
            sum(derepCounts),
            " total sequences read.")
  }
  close(f)
  return(list(uniques=derepCounts, quals=derepQuals, map=derepMap))
}
################################################################################
#' Load .uniques file
#' 
#' This is a custom interface to read.table for loading .uniques format
#'  files. The .uniques format is a delimited text file (tab-delimited by default)
#'  with the abundance in column 1 and the ASCII sequence in column 2.
#'
#' @param fl (Required). \code{character(1)}.
#'  The file path to the .uniques file,
#'  a delimited text table file containing unique sequences and their abundances.
#'  
#' @param colClasses (Optional). \code{character}.
#'  The classes of the columns in the delimited text file.
#'  Defaults to \code{c("integer", "character")}.
#'  See \code{\link[utils]{read.table}} for more information.
#' 
#' @param ... (Optional). Additional arguments passed on to \code{\link[utils]{read.table}}.
#'
#' @return Named integer vector. Named by sequence, valued by abundance.
#'
#' @export
#' 
importUniques <- function(fl, colClasses = c("integer", "character"), ...){
  unqs <- read.table(fl, colClasses = colClasses, ...)
  # Check that things worked as expected.
  if(ncol(unqs) != 2) stop(paste("Unexpected number of columns:", ncol(unqs)))
  if(class(unqs[,2]) != "character") stop(paste("Integer/sequence pairs required(1)."))
  if(class(unqs[,1]) != "integer") {
    if(class(unqs[,1]) != "numeric") {
      stop(paste("Integer/sequence pairs required(2)."))
    }
    foo <- as.integer(unqs[,1])
    if(all.equal(foo, unqs[,1])) {
      unqs[,1] <- as.integer(unqs[,1])
    } else {
      stop(paste("Integer/sequence pairs required(3)."))
    }
  }
  # Checks are done, create the uniques vector.
  rvec <- unqs[,1]
  names(rvec) <- unqs[,2]
  return(rvec)
}
################################################################################
