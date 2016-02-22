################################################################################
#' Read and Dereplicate a Fastq file.
#' 
#' This is a custom interface to \code{\link[ShortRead]{FastqStreamer}} 
#' for dereplicating amplicon sequences from fastq or compressed fastq files,
#' while also controlling peak memory requirement to support large files.
#'
#' @param fls (Required). \code{character}.
#'  The file path(s) to the fastq or fastq.gz file(s).
#'  Actually, any file format supported by \code{\link[ShortRead]{FastqStreamer}}.
#' 
#' @param n (Optional). \code{numeric(1)}.
#'  the maximum number of records (reads) to parse and dereplicate
#'  at any one time. This controls the peak memory requirement
#'  so that large fastq files are supported.
#'  Default is \code{1e6}, one-million reads.
#'  See \code{\link[ShortRead]{FastqStreamer}} for details on this parameter,
#'  which is passed on.
#' 
#' @param verbose (Optional). \code{logical(1)}.
#'  Whether or not to throw standard R \code{\link{message}}s 
#'  on the intermittent and final status of the dereplication.
#'  Default is \code{FALSE}, no messages.
#'
#' @return A \code{\link{derep-class}} object or list of such objects. 
#'
#' @export
#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead yield
#'
#' @examples 
#' # Test that chunk-size, `n`, does not affect the result.
#' testFastq = system.file("extdata", "sam1F.fastq.gz", package="dada2")
#' derep1 = derepFastq(testFastq, verbose = TRUE)
#' derep1.35 = derepFastq(testFastq, 35, TRUE)
#' all.equal(getUniques(derep1), getUniques(derep1.35)[names(getUniques(derep1))])
#' 
derepFastq <- function(fls, n = 1e6, verbose = FALSE){
  if(!is.character(fls)) {
    stop("Filenames must be provided as in character format.")
  }
  rval <- list()
  for(i in seq_along(fls)) {
    fl <- fls[[i]]
    if(verbose){
      message("Dereplicating sequence entries in Fastq file: ", fl, appendLF = TRUE)
    }
    
    f <- FastqStreamer(fl, n = n)
    suppressWarnings(fq <- yield(f))
    
    out <- qtables2(fq, FALSE)
    
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
      out <- qtables2(fq, FALSE)
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
  ###  if(qeff) derepQuals <- -10*log10(derepQuals)  # Convert back to effective Q value
    # Sort by decreasing abundance
    ord <- order(derepCounts, decreasing=TRUE)
    derepCounts <- derepCounts[ord]
    derepQuals <- derepQuals[ord,,drop=FALSE]
    derepMap <- match(derepMap, ord)
    if(verbose){
      message("Encountered ",
              length(derepCounts),
              " unique sequences from ",
              sum(derepCounts),
              " total sequences read.")
    }
    close(f)
    if(sum(tabulate(nchar(names(derepCounts)))>0) > 1) {
      warning("Not all sequences were the same length. dada(...) requires same length input sequences. See truncLen parameter in ?fastqFilter or ?fastqPairedFilter.")
    }
    derepO <- list(uniques=derepCounts, quals=derepQuals, map=derepMap)
    derepO <- as(derepO, "derep")
    rval[[i]] <- derepO
  }
  if(length(rval) == 1) rval <- rval[[1]]
  return(rval)
}
###
#' Internal tables function
#' 
#' Internal function to replicate ShortRead::tables functionality while also returning average quals and a map
#' from reads to uniques
#' 
#' @importFrom ShortRead srsort
#' @importFrom ShortRead srrank
#' @importFrom ShortRead sread
#' 
#' @keywords internal
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
#' A wrapper for writeFastq in the ShortRead package.
#' Default output format is compatible with uchime.
#' 
#' @param unqs The uniques vector. E.g. the $uniques from the output of \code{\link{derepFastq}}.
#' 
#' @param fout The file path of the output file. The file you want to write.
#' 
#' @param ids A character vector of sequence ids, one for each element in \code{unqs}.
#'  Default value is \code{NULL}, in which case a uchime-compatible ID is assigned.
#'  
#' @param mode A character string flag passed on to \code{\link[ShortRead]{writeFasta}}
#'  indicating the type of file writing mode. Default is \code{"w"}.
#'  
#' @param width The number of characters per line in the file. Default is 20000.
#'  Passed on to \code{\link[ShortRead]{writeFasta}}.
#'  
#' @param ... Additional parameters passed on to \code{\link[ShortRead]{writeFasta}}.
#' 
#' @importFrom ShortRead writeFasta
#' @importFrom ShortRead ShortRead
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings BStringSet
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' outfile <- tempfile(fileext=".fasta")
#' uniquesToFasta(getUniques(derep1), outfile)
#' uniquesToFasta(getUniques(derep1), outfile, ids=paste0("Sequence", seq(length(getUniques(derep1)))))
#' 
uniquesToFasta <- function(unqs, fout, ids=NULL, mode="w", width=20000, ...) {
  unqs <- getUniques(unqs)
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
