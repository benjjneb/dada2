################################################################################
#' Read in and dereplicate a fastq file.
#' 
#' A custom interface to \code{\link[ShortRead]{FastqStreamer}} 
#' for dereplicating amplicon sequences from fastq or compressed fastq files,
#' while also controlling peak memory requirement to support large files.
#'
#' @param fls (Required). \code{character}.
#'  The file path(s) to the fastq or fastq.gz file(s).
#'  Actually, any file format supported by \code{\link[ShortRead]{FastqStreamer}}.
#' 
#' @param n (Optional). \code{numeric(1)}.
#'  The maximum number of records (reads) to parse and dereplicate
#'  at any one time. This controls the peak memory requirement
#'  so that large fastq files are supported.
#'  Default is \code{1e6}, one-million reads.
#'  See \code{\link[ShortRead]{FastqStreamer}} for details on this parameter,
#'  which is passed on.
#' 
#' @param qualityType (Optional). \code{character(1)}.
#'  The quality encoding of the fastq file(s). "Auto" (the default) means to 
#'  attempt to auto-detect the encoding. This may fail for PacBio files with
#'  uniformly high quality scores, in which case use "FastqQuality". This
#'  parameter is passed on to \code{\link[ShortRead]{readFastq}}; see
#'  information there for details.
#' 
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, throw standard R \code{\link{message}}s 
#'  on the intermittent and final status of the dereplication.
#'
#' @return A \code{\link{derep-class}} object or list of such objects. 
#'
#' @export
#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead yield
#' @importFrom methods as
#'
#' @examples 
#' # Test that chunk-size, `n`, does not affect the result.
#' testFastq = system.file("extdata", "sam1F.fastq.gz", package="dada2")
#' derep1 = derepFastq(testFastq, verbose = TRUE)
#' derep1.35 = derepFastq(testFastq, 35, TRUE)
#' all.equal(getUniques(derep1), getUniques(derep1.35)[names(getUniques(derep1))])
#' 
derepFastq <- function(fls, n = 1e6, verbose = FALSE, qualityType = "Auto"){
  if(!is.character(fls)) { stop("Filenames must be provided in character format.") }
  rval <- list()
  for(i in seq_along(fls)) {
    fl <- fls[[i]]
    if(verbose){
      message("Dereplicating sequence entries in Fastq file: ", fl, appendLF = TRUE)
    }
    
    f <- FastqStreamer(fl, n = n)
    suppressWarnings(fq <- yield(f, qualityType = qualityType))
    
    out <- qtables2(fq, FALSE) ###ITS
    
    derepCounts <- out$uniques
    derepQuals <- out$cum_quals
    derepMap <- out$map
    while( length(suppressWarnings(fq <- yield(f, qualityType = qualityType))) ){
      # A little loop protection
      newniques = alreadySeen = NULL
      # Dot represents one turn inside the chunking loop.
      if(verbose){
        message(".", appendLF = FALSE)
      }
      out <- qtables2(fq, FALSE)
      # Augment quality matrices with NAs as needed to match ncol
      if(ncol(out$cum_quals) > ncol(derepQuals)) {
        derepQuals <- cbind(derepQuals, matrix(NA, nrow=nrow(derepQuals), ncol=(ncol(out$cum_quals)-ncol(derepQuals))))
      } else if(ncol(out$cum_quals) < ncol(derepQuals)) {
        out$cum_quals <- cbind(out$cum_quals, matrix(NA, nrow=nrow(out$cum_quals), ncol=(ncol(derepQuals)-ncol(out$cum_quals))))
      }
      # identify sequences already present in `derepCounts`
      alreadySeen <- names(out$uniques) %in% names(derepCounts)
      # Sum these values, if any
      if(any(alreadySeen)){
        sqnms = names(out$uniques)[alreadySeen]
        derepCounts[sqnms] <- derepCounts[sqnms] + out$uniques[sqnms] ###ITS
        derepQuals[sqnms,] <- derepQuals[sqnms,] + out$cum_quals[sqnms,] ###ITS
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
    derepO <- list(uniques=derepCounts, quals=derepQuals, map=derepMap)
    derepO <- as(derepO, "derep")
    rval[[i]] <- derepO
  }
  if(length(rval) == 1) {
    rval <- rval[[1]]
  } else {
    if(is.null(names(fls))) {
      names(rval) <- basename(fls)
    } else {
      names(rval) <- names(fls)
    }
  }
  return(rval)
}
###
#' Internal tables function
#' 
#' Internal function to replicate ShortRead::tables functionality while also returning average quals and a map
#' from reads to uniques
#' 
#' @param x ShortReadQ.
#'  The ShortReadQ-class object to table (or dereplicate).
#' 
#' @param qeff \code{logical(1)}.
#'  Calculate average quality by first transforming to expected error rate.
#' 
#' @return List.
#'  Matches format of derep-class object.
#' 
#' @importFrom ShortRead srsort
#' @importFrom ShortRead srrank
#' @importFrom ShortRead sread
#' @importFrom methods as
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
#' @param unqs (Required).
#'  A \code{\link{uniques-vector}} or any object that can be coerced
#'  into one with \code{\link{getUniques}}.
#'   
#' @param fout (Required).
#'  The file path of the output file.
#' 
#' @param ids (Optional). \code{character}. Default NULL.
#'  A vector of sequence ids, one for each element in \code{unqs}.
#'  If NULL, a uchime-compatible ID is assigned.
#'  
#' @param mode (Optional). Default "w".
#'  Passed on to \code{\link[ShortRead]{writeFasta}}
#'  indicating the type of file writing mode. Default is \code{"w"}.
#'  
#' @param width (Optional). Default 20000.
#'  The number of characters per line in the file. Default is effectively one line
#'  per sequence. Passed on to \code{\link[ShortRead]{writeFasta}}.
#'  
#' @param ... Additional parameters passed on to \code{\link[ShortRead]{writeFasta}}.
#' 
#' @return NULL.
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
#' uniquesToFasta(derep1, outfile)
#' uniquesToFasta(derep1, outfile, ids=paste0("Sequence", seq(length(getSequences(derep1)))))
#' 
uniquesToFasta <- function(unqs, fout, ids=NULL, mode="w", width=20000, ...) {
  unqs <- getUniques(unqs, collapse=FALSE)
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

######## DEREPFASTA DEREPFASTA DEREPFASTA ################
#' derepFasta creates a derep-class object from a fasta file, by
#' creating a corresponding fastq file with a uniform quality score
#' and calling derepFastq.
#' 
#' @param fls (Required). \code{character}.
#'  The file path(s) to the fasta or gzipped fasta file(s).
#' 
#' @param ... (Optional).
#'  Additional arguments passed on to \code{\link{derepFastq}}
#' 
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings writeXStringSet
#' 
#' @keywords internal
#' 
derepFasta <- function(fls, ...){
  if(!is.character(fls)) {
    stop("Filenames must be provided in character format.")
  }
  fastqs <- character(length(fls))
  for(i in seq_along(fls)) {
    fastq <- tempfile()
    fl <- fls[[i]]
    foo <- readDNAStringSet(fl)
    writeXStringSet(foo, fastq, format="fastq")
    fastqs[[i]] <- fastq
  }
  
  derepFastq(fastqs, ...)
}

#' Writes a named character vector of DNA sequences to a fasta file.
#' Values are the sequences, and names are used for the id lines.
#'
#' @seealso \code{\link[Biostrings]{writeXStringSet}}
#'
#' @param object (Required). A named \code{character} vector.
#' @param file (Required). The output file.
#' @param mode (Optional). Default "w". Append with "a".
#' @param width (Optional). Default 20000L. Maximum line length before newline.
#' @param ... (Optional). Additonal arguments passed to \code{\link[Biostrings]{writeXStringSet}}.
#' @return NULL.
#' @rdname writeFasta
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings writeXStringSet
#' 
setMethod("writeFasta", "character", function(object, file, mode="w", width=20000L, ...){
  append = mode == "a"
  seqs <- DNAStringSet(object)
  if(is.null(names(seqs))) { names(seqs) <- as.character(seq(length(seqs))) }
  writeXStringSet(seqs, file, ..., append = append, width=width, format = "fasta")
})

