#' Removes primers and orients reads in a consistent direction.
#' 
#' Removes primer(s) and orients the reads in input fastq file(s) (can be compressed).
#' Reads that do not contain the primer(s) are discarded.
#' Intended for use with PacBio CCS data.
#' Faster external solutions such as cutadapt or trimmomatic are recommended for short-read data.
#' 
#' @param fn (Required). \code{character}.
#' The path(s) to the input fastq file(s). Can be compressed.
#'   
#' @param fout (Required). \code{character}.
#'  The path(s) to the output fastq file(s) corresponding to the \code{fwd} input files.
#'  If directory containing the file does not exist, it will be created.
#'  Output files are gzip compressed by default.
#'   
#' @param primer.fwd (Required). \code{character}.
#'  The forward primer sequence expected to be at the beginning of the sequenced amplicon.
#'  Can contain IUPAC ambiguous nucleotide codes.
#' 
#' @param primer.rev (Optional). Default NULL.
#'  The reverse primer sequence expected to be at the end of the sequenced amplicon.
#'  Can contain IUPAC ambiguous nucleotide codes.
#'  NOTE: `primer.rev` should be provided in the orientation that would appear in a DNA sequence
#'  starting at the forward primer and being read towards the reverse primer. Thus, it is
#'  often necessary to reverse-complement the reverse primer sequence before providing it to
#'  this function.
#' 
#' @param max.mismatch (Optional). Default 2.
#'  The number of mismatches to tolerate when matching reads to primer sequences.
#'  See \code{\link[Biostrings]{vmatchPattern}} for details.
#'  
#' @param allow.indels (Optional). Default FALSE.
#'  If TRUE, indels ared allowed when matching the primer sequences to the read. If FALSE,
#'  no indels are allowed. Note that when `allow.indels=TRUE`, primer matching is significantly
#'  slower, currently about 4x slower.
#' 
#' @param trim.fwd (Optional). Default TRUE.
#'  If TRUE, reads are trimmed to the end of the forward primer, i.e. the forward
#'  primer and any preceding sequence are trimmed off.
#'  
#' @param trim.rev (Optional). Default TRUE.
#'  If TRUE, reads are trimmed to the beginning of the reverse primer, i.e. the reverse
#'  primer and any subsequent sequence are trimmed off.
#'  
#' @param orient (Optional). Default TRUE.
#'  If TRUE, reads are re-oriented if the reverse complement of the read is a better match to the
#'  provided primer sequence(s). This is recommended for PacBio CCS reads, which come in a random
#'  mix of forward and reverse-complement orientations.
#' 
#' @param compress (Optional). Default TRUE.
#'  If TRUE, the output fastq file(s) are gzipped.
#' 
#' @param verbose (Optional). Default FALSE.
#'  Whether to output status messages.  
#' 
#' @return Integer matrix. Returned invisibly (i.e. only if assigned to something).
#'  Rows correspond to the input files, columns record the number of reads.in and reads.out after
#'  discarding reads that didn't match the provided primers.
#' 
#' @importFrom Biostrings matchPattern
#' @importFrom Biostrings vmatchPattern
#' @importFrom ShortRead sread
#' @importFrom ShortRead reverseComplement
#' @importFrom ShortRead readFastq
#' @importFrom XVector rev
#' @importFrom methods as
#' @importFrom BiocGenerics end
#' @importFrom BiocGenerics width
#' @importFrom BiocGenerics start
#' 
#' @export
#' 
#' @examples
#' F27 <- "AGRGTTYGATYMTGGCTCAG"
#' R1492 <- "RGYTACCTTGTTACGACTT"
#' fn <- system.file("extdata", "samPBprimers.fastq.gz", package="dada2")
#' fn.noprime <- tempfile(fileext=".fastq.gz")
#' removePrimers(fn, fn.noprime, primer.fwd=F27, primer.rev=rc(R1492), orient=TRUE, verbose=TRUE)
#' 
# Further testing warranted for trimming of partial reverse primers, and sometimes present reverse primers
removePrimers <- function(fn, fout, 
                          primer.fwd, primer.rev=NULL, max.mismatch=2, 
                          allow.indels=FALSE, ### require.fwd=TRUE, require.rev=TRUE, 
                          trim.fwd=TRUE, trim.rev=TRUE, orient=TRUE,
                          compress=TRUE, verbose = FALSE) {
  # Check and enforce filepaths
  if(length(fn) != length(fout)) stop("Every input file must have a corresponding output file.")
  if(allow.indels) message("Primer matching with indels allowed is currently significantly (~4x) slower.")
  require.fwd <- TRUE; require.rev <- TRUE ###
  first.multi.msg <- TRUE
  odirs <- unique(dirname(fout))
  for(odir in odirs) {
    if(!dir.exists(odir)) { 
      message("Creating output directory: ", odir)
      dir.create(odir, recursive=TRUE, mode="0777")
    }
  }
  if(!all(file.exists(fn))) stop("Some input files do not exist.")
  fn <- normalizePath(fn, mustWork=TRUE)
  fout <- suppressWarnings(normalizePath(fout, mustWork=FALSE))
  if(any(duplicated(fout))) stop("All output files must be distinct.")
  if(any(fout %in% fn)) stop("Output files must be distinct from the input files.")
  # Check and enforce primers
  if(!is.character(primer.fwd)) stop("Primer sequences must be provided as base R strings.")
  if(is.null(primer.rev)) { 
    has.rev <- FALSE
  } else {
    has.rev <- TRUE
    if(!is.character(primer.rev)) stop("Primer sequences must be provided as base R strings.")
  }
  fixed.fwd <- C_isACGT(primer.fwd)
  if(has.rev) fixed.rev <- C_isACGT(primer.rev)
  rval <- matrix(0L, nrow=length(fn), ncol=2)
  colnames(rval) <- c("reads.in", "reads.out")
  rownames(rval) <- basename(fn)
  for(i in seq_along(fn)) {
    # Read in file and init filtering stats
    fq <- readFastq(fn[[i]])
    inseqs <- length(fq)
    outseqs <- 0
    rval[i,c("reads.in", "reads.out")] <- c(inseqs, outseqs)
    # Match patterns
    if(allow.indels) { # Use slower matchPattern because it supports indels
      match.fwd <- lapply(sread(fq), function(x) matchPattern(primer.fwd, x, max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.fwd))
    } else { # Use faster vmatchPattern that doesn't support indels
      match.fwd <- as(vmatchPattern(primer.fwd, sread(fq), max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.fwd), "list")
    }
    if(has.rev) {
      if(allow.indels) { # Use slower matchPattern because it supports indels
        match.rev <- lapply(sread(fq), function(x) matchPattern(primer.rev, x, max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.rev))
      } else { # Use faster vmatchPattern that doesn't support indels
        match.rev <- as(vmatchPattern(primer.rev, sread(fq), max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.rev), "list")
      }
    }
    # If orient, match reverse complement as well
    if(orient) {
      fq.rc <- reverseComplement(fq)
      if(allow.indels) { # Use slower matchPattern because it supports indels
        match.fwd.rc <- lapply(sread(fq.rc), function(x) matchPattern(primer.fwd, x, max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.fwd))
      } else { # Use faster vmatchPattern that doesn't support indels
        match.fwd.rc <- as(vmatchPattern(primer.fwd, sread(fq.rc), max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.fwd), "list")
      }
      if(has.rev) {
        if(allow.indels) { # Use slower matchPattern because it supports indels
          match.rev.rc <- lapply(sread(fq.rc), function(x) matchPattern(primer.rev, x, max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.rev))
        } else { # Use faster vmatchPattern that doesn't support indels
          match.rev.rc <- as(vmatchPattern(primer.rev, sread(fq.rc), max.mismatch=max.mismatch, with.indels=allow.indels, fixed=fixed.rev), "list")
        }
      }
    }
    # Tally up hits
    # Check for possible mis-oriented primer sequences?
    hits.fwd <- sapply(match.fwd, length)
    if(has.rev) hits.rev <- sapply(match.rev, length)
    if(!require.fwd) stop("Currently, only require.fwd=TRUE is supported.")
    if(has.rev && !require.rev) stop("Currently, only require.rev=TRUE is supported when a reverse primer sequence is provided.")
    if(require.fwd && sum(hits.fwd) == 0) { filt.print(inseqs, outseqs); next }
    if(has.rev && require.rev && sum(hits.rev) == 0) { filt.print(inseqs, outseqs); next }
    if(any(hits.fwd>1) || (has.rev && any(hits.rev>1))) {
      if(verbose && first.multi.msg) {
        message("Multiple matches to the primer(s) in some sequences. Using the longest possible match.")
        first.multi.msg <- FALSE
      }
      match.fwd[hits.fwd>1] <- sapply(match.fwd[hits.fwd>1], `[`, 1)
      if(has.rev) match.rev[hits.rev>1] <- sapply(match.rev[hits.rev>1], function(x) rev(x)[1])
    }
    if(orient) {
      hits.fwd.rc <- sapply(match.fwd.rc, length)
      if(has.rev) hits.rev.rc <- sapply(match.rev.rc, length)
      if(any(hits.fwd.rc>1) || (has.rev && any(hits.rev.rc>1))) {
        if(verbose && first.multi.msg) {
          message("Multiple matches to the primer(s) in some sequences. Using the longest possible match.")
          first.multi.msg <- FALSE
        }
        match.fwd.rc[hits.fwd.rc>1] <- sapply(match.fwd.rc[hits.fwd.rc>1], `[`, 1)
        if(has.rev) match.rev.rc[hits.rev.rc>1] <- sapply(match.rev.rc[hits.rev.rc>1], function(x) rev(x)[1])
      }
    }
    # If orient, replace non-matches with rc matches where they exist
    if(orient) {
      flip <- !hits.fwd & hits.fwd.rc
      if(any(flip) && verbose) cat(sum(flip), "sequences out of", length(flip), "are being reverse-complemented.\n")
      fq[flip] <- fq.rc[flip]
      match.fwd[flip] <- match.fwd.rc[flip]
      hits.fwd <- sapply(match.fwd, length)
      if(has.rev) {
        match.rev[flip] <- match.rev.rc[flip]
        hits.rev <- sapply(match.rev, length)
      }
    }
    # If require, remove sequences w/o forward and reverse hits
    keep <- rep(TRUE, length(fq))
    if(require.fwd) keep <- keep & (hits.fwd > 0)
    if(has.rev && require.rev) keep <- keep & (hits.rev > 0)
    if(!all(keep)) {
      fq <- fq[keep]
      match.fwd <- match.fwd[keep]
      if(has.rev) match.rev <- match.rev[keep]
    }
    # If trim, narrow to the desired subsequence
    if(trim.fwd) {
      first <- sapply(match.fwd, end) + 1
    } else {
      first <- rep(1L, length(fq))
    }
    if(has.rev && trim.rev) {
      last <- sapply(match.rev, start) - 1
    } else {
      last <- width(fq)
    }
    keep <- last > first
    if(!all(keep)) first <- first[keep]; last <- last[keep]; fq <- fq[keep]
    fq <- narrow(fq, first, last) # Need to handle zero case gracefully, w/ informative error
    # Delete fout if it already exists (since writeFastq doesn't overwrite)
    if(file.exists(fout[[i]])) {
      if(file.remove(fout[[i]])) {
        if(verbose) message("Overwriting file:", fout[[i]])
      } else {
        stop("Failed to overwrite file:", fout[[i]])
      }
    }
    writeFastq(fq, fout[[i]], "w", compress=compress)
    outseqs <- length(fq)
    rval[i,c("reads.in", "reads.out")] <- c(inseqs, outseqs)
    if(verbose) filt.print(inseqs, outseqs)
  }
  if(all(rval[,"reads.out"]==0)) {
    warning("No reads passed the primer detection.")
  } else if(any(rval[,"reads.out"]==0)) {
    message("Some input samples had no reads pass the primer detection.")
  }
  return(invisible(rval))
}  

filt.print <- function(inseqs, outseqs) {
  outperc <- round(outseqs * 100 / inseqs, 1)
  outperc <- paste(" (", outperc, "%)", sep="")
  message("Read in ", inseqs, ", output ", outseqs, outperc, " filtered sequences.", sep="")
}

#' Filter and trim fastq file(s).
#' 
#' Filters and trims an input fastq file(s) (can be compressed)
#' based on several user-definable criteria, and outputs fastq file(s)
#' (compressed by default) containing those trimmed reads which passed the filters. Corresponding
#' forward and reverse fastq file(s) can be provided as input, in which case filtering
#' is performed on the forward and reverse reads independently, and both reads must pass for
#' the read pair to be output.
#' 
#' \code{filterAndTrim} is a multithreaded convenience interface for the \code{\link{fastqFilter}}
#' and \code{\link{fastqPairedFilter}} filtering functions.
#' Note that error messages and tracking are not handled gracefully when using the multithreading
#' functionality. If errors arise, it is recommended to re-run without multithreading to
#' troubleshoot the issue.
#' 
#' @param fwd (Required). \code{character}.
#'  The file path(s) to the fastq file(s), or the directory containing the fastq file(s). 
#'  Compressed file formats such as .fastq.gz and .fastq.bz2 are supported.
#'     
#' @param filt (Required). \code{character}.
#'  The path(s) to the output filtered file(s) corresponding to the \code{fwd} input files, or a directory
#'  that will contain those files.
#'  If containing directory does not exist, it will be created.
#'   
#' @param rev (Optional). Default NULL.
#'  The file path(s) to the reverse fastq file(s) from paired-end sequence data corresponding to those
#'  provided to the \code{fwd} argument, or the directory containing those fastq file(s). 
#'  Compressed file formats such as .fastq.gz and .fastq.bz2 are supported.
#'  If NULL, the \code{fwd} files are processed as single-reads.
#' 
#' @param filt.rev (Optional). Default NULL, but required if \code{rev} is provided.
#'  The path(s) to the output filtered file(s) corresponding to the \code{rev} input files, or a directory
#'  that will contain those files.
#'  If containing directory does not exist, it will be created.
#'   
#' @param compress (Optional). Default TRUE.
#'  If TRUE, the output fastq file(s) are gzipped.
#' 
#' \strong{FILTERING AND TRIMMING PARAMETERS ---------}   
#' 
#' \strong{Note:} When filtering paired reads...
#' If a length 1 vector is provided, the same parameter value is used for the forward and reverse reads.
#' If a length 2 vector is provided, the first value is used for the forward reads, and the second 
#'   for the reverse reads.
#' 
#' @param truncQ (Optional). Default 2.
#'  Truncate reads at the first instance of a quality score less than or equal to \code{truncQ}.
#'  
#' @param truncLen (Optional). Default 0 (no truncation).
#'  Truncate reads after \code{truncLen} bases. Reads shorter than this are discarded.
#'  
#' @param trimLeft (Optional). Default 0.
#'  The number of nucleotides to remove from the start of each read. If both \code{truncLen} and 
#'  \code{trimLeft} are provided, filtered reads will have length \code{truncLen-trimLeft}.
#'  
#' @param trimRight (Optional). Default 0.
#'  The number of nucleotides to remove from the end of each read. If both \code{truncLen} and 
#'  \code{trimRight} are provided, truncation will be performed after \code{trimRight} is enforced.
#'  
#' @param maxLen (Optional). Default Inf (no maximum).
#'  Remove reads with length greater than maxLen. maxLen is enforced \strong{before} trimming and truncation.   
#'  
#' @param minLen (Optional). Default 20.
#'  Remove reads with length less than minLen. minLen is enforced \strong{after} trimming and truncation.   
#'  
#' @param maxN (Optional). Default 0.
#'  After truncation, sequences with more than \code{maxN} Ns will be discarded. 
#'  Note that \code{\link{dada}} does not allow Ns.
#'  
#' @param minQ (Optional). Default 0.
#'  After truncation, reads contain a quality score less than \code{minQ} will be discarded.
#'
#' @param maxEE (Optional). Default \code{Inf} (no EE filtering).
#'  After truncation, reads with higher than \code{maxEE} "expected errors" will be discarded.
#'  Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
#'  
#' @param rm.phix (Optional). Default TRUE.
#'  If TRUE, discard reads that match against the phiX genome, as determined by \code{\link{isPhiX}}.
#'
#' @param rm.lowcomplex (Optional). Default 0.
#'  If greater than 0, reads with an effective number of kmers less than this value will be removed.
#'  The effective number of kmers is determined by \code{\link{seqComplexity}} using a Shannon information
#'  approximation. The default kmer-size is 2, and therefore perfectly random sequences will approach an
#'  effective kmer number of 16 = 4 (nucleotides) ^ 2 (kmer size).
#'
#' @param orient.fwd (Optional). Default NULL.
#'  A character string present at the start of valid reads. Only allows unambiguous nucleotides. 
#'  This string is compared to the start of each read, and the reverse complement of each read.
#'  If it exactly matches the start of the read, the read is kept.
#'  If it exactly matches the start of the reverse-complement read, the read is reverse-complemented and kept.
#'  Otherwise the read if filtered out.
#'  For paired reads, the string is compared to the start of the forward and reverse reads, and if it matches
#'  the start of the reverse read the reaads are swapped and kept.
#'  The primary use of this parameter is to unify the orientation of amplicon sequencing libraries that
#'  are a mixture of forward and reverse orientations, and that include the forward primer on the reads.
#'
#' @param matchIDs (Optional). Default FALSE. Paired-read filtering only.
#'  Whether to enforce matching between the id-line sequence identifiers of the forward and reverse fastq files.
#'    If TRUE, only paired reads that share id fields (see below) are output.
#'    If FALSE, no read ID checking is done.
#'  Note: \code{matchIDs=FALSE} essentially assumes matching order between forward and reverse reads. If that
#'    matched order is not present future processing steps may break (in particular \code{\link{mergePairs}}).
#'
#' @param id.sep (Optional). Default "\\s" (white-space). Paired-read filtering only.
#'  The separator between fields in the id-line of the input fastq files. Passed to the \code{\link{strsplit}}.
#' 
#' @param id.field (Optional). Default NULL (automatic detection). Paired-read filtering only.
#'  The field of the id-line containing the sequence identifier.
#'  If NULL (the default) and matchIDs is TRUE, the function attempts to automatically detect
#'    the sequence identifier field under the assumption of Illumina formatted output.
#'
#' @param multithread (Optional). Default is FALSE.
#'  If TRUE, input files are filtered in parallel via \code{\link[parallel]{mclapply}}.
#'  If an integer is provided, it is passed to the \code{mc.cores} argument of \code{\link[parallel]{mclapply}}.
#'  Note that the parallelization here is by forking, and each process is loading another fastq file into
#'  memory. This option is ignored in Windows, as Windows does not support forking, with \code{mc.cores} set to 1.
#'	If memory is an issue, execute in a clean environment and reduce the chunk size \code{n} and/or
#'  the number of threads.
#'   
#' @param n (Optional). Default \code{1e5}.
#' The number of records (reads) to read in and filter at any one time. 
#' This controls the peak memory requirement so that very large fastq files are supported. 
#' See \code{\link[ShortRead]{FastqStreamer}} for details.
#'
#' @param OMP (Optional). Default TRUE.
#'  Whether or not to use OMP multithreading when calling \code{\link[ShortRead]{FastqStreamer}}. 
#'  Should be set to FALSE if calling this function within a parallelized chunk of code.
#'  If \code{multithread=TRUE}, this argument will be coerced to FALSE.
#'  
#' @param qualityType (Optional). \code{character(1)}.
#'  The quality encoding of the fastq file(s). "Auto" (the default) means to
#'  attempt to auto-detect the encoding. This may fail for PacBio files with
#'  uniformly high quality scores, in which case use "FastqQuality". This
#'  parameter is passed on to \code{\link[ShortRead]{readFastq}}; see
#'  information there for details.
#' 
#' @param verbose (Optional). Default FALSE.
#'  Whether to output status messages.  
#' 
#' @return Integer matrix. Returned invisibly (i.e. only if assigned to something).
#'  Rows correspond to the input files, columns record the reads.in and reads.out after filtering.
#' 
#' @seealso 
#'  \code{\link{fastqFilter}}
#'  \code{\link{fastqPairedFilter}}
#'  \code{\link[ShortRead]{FastqStreamer}}
#' 
#' @importFrom parallel mcmapply
#' @importFrom parallel detectCores
#' @importFrom methods as
#' @importFrom methods is
#' 
#' @export
#' 
#' @examples
#' testFastqs = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"), 
#'                system.file("extdata", "sam2F.fastq.gz", package="dada2"))
#' filtFastqs <- c(tempfile(fileext=".fastq.gz"), tempfile(fileext=".fastq.gz"))
#' filterAndTrim(testFastqs, filtFastqs, maxN=0, maxEE=2, verbose=TRUE)
#' filterAndTrim(testFastqs, filtFastqs, truncQ=2, truncLen=200, rm.phix=TRUE, rm.lowcomplex=8)
#' 
filterAndTrim <- function(fwd, filt, rev=NULL, filt.rev=NULL, compress=TRUE,
                        truncQ=2, truncLen=0, trimLeft=0, trimRight=0, maxLen=Inf, minLen=20, 
                        maxN=0, minQ=0, maxEE=Inf, rm.phix=TRUE, rm.lowcomplex=0, orient.fwd=NULL,
                        matchIDs=FALSE, id.sep="\\s", id.field=NULL,
                        multithread=FALSE, n = 1e5, OMP=TRUE, qualityType = "Auto", verbose = FALSE) {
  PAIRED <- FALSE
  # Validate inputs
  if(!(is.character(fwd) && is.character(filt))) stop("File paths must be provided as character vectors.")
  if(length(fwd)==1 && dir.exists(fwd)) fwd <- parseFastqDirectory(fwd)
  if(!all(file.exists(fwd))) stop("Some input files do not exist.")
  if(length(filt)==1 && length(fwd)>1) filt <- file.path(filt, basename(fwd)) # Interpret filt as a directory
  if(length(fwd) != length(filt)) stop("Every input file must have a corresponding output file.")
  odirs <- unique(dirname(filt))
  for(odir in odirs) {
    if(!dir.exists(odir)) { 
      message("Creating output directory: ", odir)
      dir.create(odir, recursive=TRUE, mode="0777")
    }
  }
  fwd <- normalizePath(fwd, mustWork=TRUE)
  filt <- suppressWarnings(normalizePath(filt, mustWork=FALSE))
  if(any(duplicated(filt))) stop("All output files must be distinct.")
  if(any(filt %in% fwd)) stop("Output files must be distinct from the input files.")
  if(!is.null(rev)) {
    PAIRED <- TRUE
    if(is.null(filt.rev)) stop("Output files for the reverse reads are required.")
    if(!(is.character(rev) && is.character(filt.rev))) stop("File paths (rev/filt.rev) must be provided as character vectors.")
    if(length(rev)==1 && dir.exists(rev)) rev <- parseFastqDirectory(rev)
    if(!all(file.exists(rev))) stop("Some input files (rev) do not exist.")
    if(length(rev) != length(fwd)) stop("Paired forward and reverse input files must correspond.")
    if(length(filt.rev)==1 && length(rev)>1) filt.rev <- file.path(filt.rev, basename(rev)) # Interpret filt.rev as a directory
    if(length(rev) != length(filt.rev)) stop("Every input file (rev) must have a corresponding output file (filt.rev).")
    odirs <- unique(dirname(filt.rev))
    for(odir in odirs) {
      if(!dir.exists(odir)) { 
        message("Creating output directory:", odir)
        dir.create(odir, recursive=TRUE, mode="0777") 
      }
    }
    rev <- suppressWarnings(normalizePath(rev, mustWork=TRUE))
    filt.rev <- suppressWarnings(normalizePath(filt.rev, mustWork=FALSE))
    if(any(duplicated(c(filt, filt.rev)))) stop("All output files must be distinct.")
    if(any(c(filt,filt.rev) %in% c(fwd, rev))) stop("Output files must be distinct from the input files.")
  }
  # Parse multithreading
  if(multithread){
    OMP <- FALSE
    ncores <- detectCores()
    if(is.numeric(multithread)) ncores <- multithread
    if(is.na(ncores)) ncores <- 1
    if(ncores>1) verbose <- FALSE
  }else{
    ncores <- 1
  }
  
  # Filter and Trim
  # Switch between using mcmapply or clusterMap depending on availability of forking
  if(PAIRED) {
    if(!multithread | .Platform$OS.type == "unix"){
      rval <- mcmapply(fastqPairedFilter,
                       mapply(c, fwd, rev, SIMPLIFY=FALSE), mapply(c, filt, filt.rev, SIMPLIFY=FALSE),
                       MoreArgs = list(truncQ=truncQ, truncLen=truncLen, trimLeft=trimLeft, trimRight=trimRight,
                                       maxLen=maxLen, minLen=minLen, maxN=maxN, minQ=minQ, maxEE=maxEE,
                                       rm.phix=rm.phix, rm.lowcomplex=rm.lowcomplex, orient.fwd=orient.fwd,
                                       matchIDs=matchIDs, id.sep=id.sep, id.field=id.field, n=n, OMP=OMP,
                                       qualityType=qualityType, compress=compress, verbose=verbose),
                       mc.cores=ncores, mc.silent=TRUE)
    }else if(.Platform$OS.type == "windows"){
      cl <- makeCluster(ncores, type = "PSOCK")
      rval <- clusterMap(cl = cl,
                         fun = fastqPairedFilter,
                         mapply(c, fwd, rev, SIMPLIFY=FALSE), mapply(c, filt, filt.rev, SIMPLIFY=FALSE),
                         MoreArgs = list(truncQ=truncQ, truncLen=truncLen, trimLeft=trimLeft, trimRight=trimRight, 
                                         maxLen=maxLen, minLen=minLen, maxN=maxN, minQ=minQ, maxEE=maxEE, 
                                         rm.phix=rm.phix, rm.lowcomplex=rm.lowcomplex, orient.fwd=orient.fwd,
                                         matchIDs=matchIDs, id.sep=id.sep, id.field=id.field, n=n, OMP=OMP, 
                                         qualityType=qualityType, compress=compress, verbose=verbose),
                         .scheduling = "dynamic")
      stopCluster(cl)
    }
  } else {
    if(!multithread | .Platform$OS.type == "unix"){
      rval <- mcmapply(fastqFilter, 
                       fwd, filt, 
                       MoreArgs = list(truncQ=truncQ, truncLen=truncLen, trimLeft=trimLeft, trimRight=trimRight,
                                       maxLen=maxLen, minLen=minLen, maxN=maxN, minQ=minQ, maxEE=maxEE, 
                                       rm.phix=rm.phix, rm.lowcomplex=rm.lowcomplex, orient.fwd=orient.fwd,
                                       n=n, OMP=OMP, qualityType=qualityType, compress=compress, verbose=verbose),
                       mc.cores=ncores, mc.silent=TRUE)
    }else if(.Platform$OS.type == "windows"){
      cl <- makeCluster(ncores, type = "PSOCK")
      rval <- clusterMap(cl = cl,
                         fun = fastqFilter,
                         fwd, filt,
                         MoreArgs = list(truncQ=truncQ, truncLen=truncLen, trimLeft=trimLeft, trimRight=trimRight,
                                         maxLen=maxLen, minLen=minLen, maxN=maxN, minQ=minQ, maxEE=maxEE, 
                                         rm.phix=rm.phix, rm.lowcomplex=rm.lowcomplex, orient.fwd=orient.fwd,
                                         n=n, OMP=OMP, qualityType=qualityType, compress=compress, verbose=verbose),
                         mc.cores=ncores, mc.silent=TRUE)
      stopCluster(cl)
    }
  }
  
  # Check if expected matrix was returned, if not there are errors
  if(!is(rval, "matrix")) {
    if(is(rval, "list")) { # Mix of errors and not
      rval <- unlist(rval[sapply(rval, is.character)])
    }
    if(length(rval)>5) rval <- rval[1:5]
    stop("These are the errors (up to 5) encountered in individual cores...\n", rval)
  }
  # Check if all input files generated a return (to catch poorly behaving out-of-memory errors)
  if(ncol(rval) != length(fwd)) {
    stop("Some input files were not processed, perhaps due to memory issues. Consider lowering ncores.")
  }
  colnames(rval) <- basename(fwd)
  if(all(rval["reads.out",]==0)) {
    warning("No reads passed the filter. Please revisit your filtering parameters.")
  } else if(any(rval["reads.out",]==0)) {
    message("Some input samples had no reads pass the filter.")
  }
  return(invisible(t(rval)))
}
#' Filter and trim a fastq file.
#' 
#' fastqFilter takes an input fastq file (can be compressed), filters it based on several
#' user-definable criteria, and outputs those reads which pass the filter
#' to a new fastq file (also can be compressed). Several functions in the \code{ShortRead}
#' package are leveraged to do this filtering.
#' 
#' @param fn (Required). The path to the input fastq file.
#'   
#' @param fout (Required). The path to the output file.
#'  Note that by default (\code{compress=TRUE}) the output fastq file is gzipped.
#' 
#' @param truncQ (Optional). Default 2.
#'  Truncate reads at the first instance of a quality score less than or equal to \code{truncQ}.
#'  
#' @param truncLen (Optional). Default 0 (no truncation).
#'  Truncate reads after \code{truncLen} bases. Reads shorter than this are discarded.
#'  
#' @param maxLen (Optional). Default Inf (no maximum).
#'  Remove reads with length greater than maxLen. maxLen is enforced on the raw reads.   
#'  
#' @param minLen (Optional). Default 20.
#'  Remove reads with length less than minLen. minLen is enforced after all other trimming and truncation.   
#'  
#' @param trimLeft (Optional). Default 0.
#'  The number of nucleotides to remove from the start of each read. If both \code{truncLen} and 
#'  \code{trimLeft} are provided, filtered reads will have length \code{truncLen-trimLeft}.
#'  
#' @param trimRight (Optional). Default 0.
#'  The number of nucleotides to remove from the end of each read. If both \code{truncLen} and 
#'  \code{trimRight} are provided, truncation will be performed after \code{trimRight} is enforced.
#'  
#' @param maxN (Optional). Default 0.
#'  After truncation, sequences with more than \code{maxN} Ns will be discarded. 
#'  Note that \code{\link{dada}} currently does not allow Ns.
#'  
#' @param minQ (Optional). Default 0.
#'  After truncation, reads contain a quality score below minQ will be discarded.
#'
#' @param maxEE (Optional). Default \code{Inf} (no EE filtering).
#'  After truncation, reads with higher than maxEE "expected errors" will be discarded.
#'  Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
#'  
#' @param rm.phix (Optional). Default TRUE.
#'  If TRUE, discard reads that match against the phiX genome, as determined by 
#'  \code{\link{isPhiX}}.
#'
#' @param rm.lowcomplex (Optional). Default 0.
#'  If greater than 0, reads with an effective number of kmers less than this value will be removed.
#'  The effective number of kmers is determined by \code{\link{seqComplexity}} using a Shannon information
#'  approximation. The default kmer-size is 2, and therefore perfectly random sequences will approach an
#'  effective kmer number of 16 = 4 (nucleotides) ^ 2 (kmer size).
#'
#' @param orient.fwd (Optional). Default NULL.
#'  A character string present at the start of valid reads. Only allows unambiguous nucleotides. 
#'  This string is compared to the start of each read, and the reverse complement of each read.
#'  If it exactly matches the start of the read, the read is kept.
#'  If it exactly matches the start of the reverse-complement read, the read is reverse-complemented and kept.
#'  Otherwise the read if filtered out.
#'  The primary use of this parameter is to unify the orientation of amplicon sequencing libraries that
#'  are a mixture of forward and reverse orientations, and that include the forward primer on the reads.
#'  
#' @param n (Optional). The number of records (reads) to read in and filter at any one time. 
#'  This controls the peak memory requirement so that very large fastq files are supported. 
#'  Default is \code{1e6}, one-million reads. See \code{\link[ShortRead]{FastqStreamer}} for details.
#'
#' @param OMP (Optional). Default TRUE.
#'  Whether or not to use OMP multithreading when calling \code{\link[ShortRead]{FastqStreamer}}. 
#'  Set this to FALSE if calling this function within a parallelized chunk of code 
#'  (eg. within \code{\link[parallel]{mclapply}}).
#' 
#' @param qualityType (Optional). \code{character(1)}.
#'  The quality encoding of the fastq file(s). "Auto" (the default) means to
#'  attempt to auto-detect the encoding. This may fail for PacBio files with
#'  uniformly high quality scores, in which case use "FastqQuality". This
#'  parameter is passed on to \code{\link[ShortRead]{readFastq}}; see
#'  information there for details.
#' 
#' 
#' @param compress (Optional). Default TRUE.
#'  Whether the output fastq file should be gzip compressed.
#' 
#' @param verbose (Optional). Default FALSE.
#'  Whether to output status messages.
#' 
#' @param ... (Optional). Arguments passed on to \code{\link{isPhiX}}.
#' 
#' @return \code{integer(2)}.
#'  The number of reads read in, and the number of reads that passed the filter and were output.
#' 
#' @seealso 
#'  \code{\link{fastqPairedFilter}}
#'  \code{\link[ShortRead]{FastqStreamer}}
#'  \code{\link[ShortRead]{trimTails}}
#' 
#' @export
#' 
#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead yield
#' @importFrom ShortRead writeFastq
#' @importFrom ShortRead trimTails
#' @importFrom ShortRead encoding
#' @importFrom ShortRead narrow
#' @importFrom IRanges narrow
#' @importFrom Biostrings quality
#' @importFrom Biostrings width
#' @importFrom Biostrings end
#' @importFrom methods as
#' 
#' @examples
#' testFastq = system.file("extdata", "sam1F.fastq.gz", package="dada2")
#' filtFastq <- tempfile(fileext=".fastq.gz")
#' fastqFilter(testFastq, filtFastq, maxN=0, maxEE=2)
#' fastqFilter(testFastq, filtFastq, trimLeft=10, truncLen=200, maxEE=2, verbose=TRUE)
#' 
fastqFilter <- function(fn, fout, truncQ = 2, truncLen = 0, maxLen = Inf, minLen = 20, trimLeft = 0, trimRight = 0, maxN = 0, minQ = 0, maxEE = Inf, rm.phix = TRUE, rm.lowcomplex = 0, orient.fwd = NULL, n = 1e6, OMP = TRUE, qualityType = "Auto", compress = TRUE, verbose = FALSE, ...){
  if(!OMP) {
    ompthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
    on.exit(.Call(ShortRead:::.set_omp_threads, ompthreads))
  }
  if(any(sapply(list(truncQ, truncLen, maxLen, minLen, trimLeft, trimRight, maxN, minQ, maxEE), length) > 1)) {
    stop("Filtering and trimming arguments should be of length 1 when processing single-end (rather than paired-end) data.")
  }
  
  start <- max(1, trimLeft + 1, na.rm=TRUE)
  end <- truncLen
  if(end < start) { end = NA }
  end <- end - start + 1
  
  if(fn == fout) { stop("The output and input files must be different.") }
  
  ## iterating over an entire file using fastq streaming
  f <- FastqStreamer(fn, n = n)
  on.exit(close(f))
  
  # Delete fout if it already exists (since writeFastq doesn't overwrite)
  if(file.exists(fout)) {
    if(file.remove(fout)) {
      if(verbose) message("Overwriting file:", fout)
    } else {
      stop("Failed to overwrite file:", fout)
    }
  }
  
  first=TRUE
  inseqs = 0
  outseqs = 0
  while( length(suppressWarnings(fq <- yield(f, qualityType = qualityType))) ){
    inseqs <- inseqs + length(fq)
    
    # Enforce and orient on orient.fwd
    if(!is.null(orient.fwd)) {
      if(!C_isACGT(orient.fwd)) stop("Non-ACGT characters detected in orient.fwd")
      barlen <- nchar(orient.fwd)
      fq.rc <- reverseComplement(fq)
      keepF <- narrow(sread(fq),1,barlen) == orient.fwd
      keepR <- narrow(sread(fq.rc),1,barlen) == orient.fwd & !keepF
      fq <- ShortReadQ(sread=c(sread(fq[keepF]), sread(fq.rc[keepR])),
                       quality=c(quality(quality(fq[keepF])), quality(quality(fq.rc[keepR]))),
                       id=c(id(fq[keepF]), id(fq.rc[keepR])))
    }
    # Enforce maxLen
    if(is.finite(maxLen)) { fq <- fq[width(fq) <= maxLen] }
    # Trim left
    fq <- fq[width(fq) >= start]
    fq <- narrow(fq, start = start, end = NA)
    # Trim right
    if(trimRight > 0) {
      fq <- fq[width(fq) > trimRight]
      fq <- narrow(fq, start=NA, end=width(fq)-trimRight)
    }
    # Trim on truncQ 
    # Convert numeric quality score to the corresponding ascii character
    enc <- encoding(quality(fq))
    if(is.numeric(truncQ)) {
      ind <- which(enc==truncQ)
      if(length(ind) != 1) stop("Encoding for this truncQ value not found.")
      truncQ <- names(enc)[[ind]]
    }
    if(length(fq) > 0) fq <- trimTails(fq, 1, truncQ)
    truncQ <- enc[truncQ] # Convert back to integer
    # Filter any with less than required length
    if(!is.na(end)) { fq <- fq[width(fq) >= end] }
    # Truncate to truncLen
    fq <- narrow(fq, start = 1, end = end)
    # Enforce minLen
    fq <- fq[width(fq) >= minLen]
    
    # Filter based on minQ and Ns and maxEE
    fq <- fq[.nFilter(fq, maxN)]
    keep <- rep(TRUE, length(fq))
    qq <- as(quality(fq), "matrix")
    if(minQ > truncQ) keep <- keep & (apply(qq, 1, min, na.rm=TRUE)>minQ) # Prob a faster trimTails trick
    if(maxEE < Inf) {
      keep <- keep & C_matrixEE(qq) <= maxEE
    }
    fq <- fq[keep]
    
    # Remove phiX
    if(rm.phix) {
      is.phi <- isPhiX(as(sread(fq), "character"), ...)
      fq <- fq[!is.phi]
    }
    
    # Remove low complexity
    if(rm.lowcomplex > 0) {
      sqcmplx <- seqComplexity(sread(fq), ...)
      fq <- fq[sqcmplx >= rm.lowcomplex]
    }
    
    outseqs <- outseqs + length(fq)
    
    if(first) {
      writeFastq(fq, fout, "w", compress=compress)
      first=FALSE
    } else {
      writeFastq(fq, fout, "a", compress=compress)
    }
  }
  
  if(verbose) {
    outperc <- round(outseqs * 100 / inseqs, 1)
    outperc <- paste(" (", outperc, "%)", sep="")
    message("Read in ", inseqs, ", output ", outseqs, outperc, " filtered sequences.", sep="")
  }
  
  if(outseqs==0) {
    message(paste("The filter removed all reads:", fout, "not written."))
    file.remove(fout)
  }
  
  return(invisible(c(reads.in=inseqs, reads.out=outseqs)))
}
#' Filters and trims paired forward and reverse fastq files.
#' 
#' fastqPairedFilter filters pairs of input fastq files (can be compressed) based on several
#' user-definable criteria, and outputs those read pairs which pass the filter in \strong{both} directions
#' to two new fastq file (also can be compressed). Several functions
#' in the \code{ShortRead} package are leveraged to do this filtering. The filtered forward/reverse reads
#' remain identically ordered.
#' 
#' @param fn (Required). A \code{character(2)} naming the paths to the (forward,reverse) fastq files.
#'   
#' @param fout (Required). A \code{character(2)} naming the paths to the (forward,reverse) output files.
#'  Note that by default (\code{compress=TRUE}) the output fastq files are gzipped.
#' 
#' \strong{FILTERING AND TRIMMING ARGUMENTS}   
#' 
#' If a length 1 vector is provided, the same parameter value is used for the forward and reverse reads.
#' If a length 2 vector is provided, the first value is used for the forward reads, and the second 
#'   for the reverse reads.
#' 
#' @param truncQ (Optional). Default 2.
#'  Truncate reads at the first instance of a quality score less than or equal to \code{truncQ}.
#'  
#' @param truncLen (Optional). Default 0 (no truncation).
#'  Truncate reads after \code{truncLen} bases. Reads shorter than this are discarded.
#'  
#' @param maxLen (Optional). Default Inf (no maximum).
#'  Remove reads with length greater than maxLen. maxLen is enforced on the raw reads.   
#'  
#' @param minLen (Optional). Default 20.
#'  Remove reads with length less than minLen. minLen is enforced after all other trimming and truncation.   
#'  
#' @param trimLeft (Optional). Default 0.
#'  The number of nucleotides to remove from the start of each read. If both \code{truncLen} and 
#'  \code{trimLeft} are provided, filtered reads will have length \code{truncLen-trimLeft}.
#'  
#' @param trimRight (Optional). Default 0.
#'  The number of nucleotides to remove from the end of each read. If both \code{truncLen} and 
#'  \code{trimRight} are provided, truncation will be performed after \code{trimRight} is enforced.
#'  
#' @param maxN (Optional). Default 0.
#'  After truncation, sequences with more than \code{maxN} Ns will be discarded. 
#'  Note that \code{\link{dada}} currently does not allow Ns.
#'  
#' @param minQ (Optional). Default 0.
#'  After truncation, reads contain a quality score below minQ will be discarded.
#'
#' @param maxEE (Optional). Default \code{Inf} (no EE filtering).
#'  After truncation, reads with higher than maxEE "expected errors" will be discarded.
#'  Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
#'  
#' @param rm.phix (Optional). Default TRUE.
#'  If TRUE, discard reads that match against the phiX genome, as determined by 
#'  \code{\link{isPhiX}}.
#'  
#' @param rm.lowcomplex (Optional). Default 0.
#'  If greater than 0, reads with an effective number of kmers less than this value will be removed.
#'  The effective number of kmers is determined by \code{\link{seqComplexity}} using a Shannon information
#'  approximation. The default kmer-size is 2, and therefore perfectly random sequences will approach an
#'  effective kmer number of 16 = 4 (nucleotides) ^ 2 (kmer size).
#'
#' @param orient.fwd (Optional). Default NULL.
#'  A character string present at the start of valid reads. Only allows unambiguous nucleotides. 
#'  This string is compared to the start of the forward and reverse reads. 
#'  If it exactly matches the start of the forward read, the read is kept.
#'  If it exactly matches the start of the reverse read, the fwd/rev reads are swapped.
#'  Otherwise the read if filtered out.
#'  The primary use of this parameter is to unify the orientation of amplicon sequencing libraries that
#'  are a mixture of forward and reverse orientations, and that include the forward primer on the reads.
#'
#' \strong{ID MATCHING ARGUMENTS}   
#' 
#'  The following optional arguments enforce matching between the sequence identification
#'  strings in the forward and reverse reads, and can automatically detect and match ID fields in 
#'  Illumina format, e.g: EAS139:136:FC706VJ:2:2104:15343:197393. ID matching is not required
#'  when using standard Illumina output fastq files.
#' 
#' @param matchIDs (Optional). Default FALSE.
#'  Whether to enforce matching between the id-line sequence identifiers of the forward and reverse fastq files.
#'    If TRUE, only paired reads that share id fields (see below) are output.
#'    If FALSE, no read ID checking is done.
#'  Note: \code{matchIDs=FALSE} essentially assumes matching order between forward and reverse reads. If that
#'    matched order is not present future processing steps may break (in particular \code{\link{mergePairs}}).
#'
#' @param id.sep (Optional). Default "\\s" (white-space).
#'  The separator between fields in the id-line of the input fastq files. Passed to the \code{\link{strsplit}}.
#' 
#' @param id.field (Optional). Default NULL (automatic detection).
#'  The field of the id-line containing the sequence identifier.
#'  If NULL (the default) and matchIDs is TRUE, the function attempts to automatically detect
#'    the sequence identifier field under the assumption of Illumina formatted output.
#'
#' @param n (Optional). The number of records (reads) to read in and filter at any one time.
#'  This controls the peak memory requirement so that very large fastq files are supported.
#'  Default is \code{1e6}, one-million reads. See \code{\link[ShortRead]{FastqStreamer}} for details.
#'  
#' @param OMP (Optional). Default TRUE.
#'  Whether or not to use OMP multithreading when calling \code{\link[ShortRead]{FastqStreamer}}. 
#'  Set this to FALSE if calling this function within a parallelized chunk of code 
#'  (eg. within \code{\link[parallel]{mclapply}}).
#' 
#' @param qualityType (Optional). \code{character(1)}.
#'  The quality encoding of the fastq file(s). "Auto" (the default) means to
#'  attempt to auto-detect the encoding. This parameter is passed on to
#'  \code{\link[ShortRead]{readFastq}}; see information there for details.
#' 
#' @param compress (Optional). Default TRUE.
#'  Whether the output fastq files should be gzip compressed.
#' 
#' @param verbose (Optional). Default FALSE.
#'  Whether to output status messages.  
#' 
#' @param ... (Optional). Arguments passed on to \code{\link{isPhiX}} or \code{\link{seqComplexity}}.
#' 
#' @return \code{integer(2)}.
#'  The number of reads read in, and the number of reads that passed the filter and were output.
#' 
#' @seealso 
#' \code{\link{fastqFilter}}
#' \code{\link[ShortRead]{FastqStreamer}}
#' \code{\link[ShortRead]{trimTails}}
#' 
#' @export
#' 
#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead yield
#' @importFrom ShortRead writeFastq
#' @importFrom ShortRead trimTails
#' @importFrom ShortRead encoding
#' @importFrom ShortRead append
#' @importFrom ShortRead ShortReadQ
#' @importFrom ShortRead narrow
#' @importFrom IRanges narrow
#' @importFrom Biostrings quality
#' @importFrom Biostrings width
#' @importFrom Biostrings end
#' @importFrom methods as
#' 
#' @examples
#'
#' testFastqF = system.file("extdata", "sam1F.fastq.gz", package="dada2")
#' testFastqR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
#' filtFastqF <- tempfile(fileext=".fastq.gz")
#' filtFastqR <- tempfile(fileext=".fastq.gz")
#' fastqPairedFilter(c(testFastqF, testFastqR), c(filtFastqF, filtFastqR), maxN=0, maxEE=2)
#' fastqPairedFilter(c(testFastqF, testFastqR), c(filtFastqF, filtFastqR), trimLeft=c(10, 20),
#'                     truncLen=c(240, 200), maxEE=2, rm.phix=TRUE, rm.lowcomplex=5, kmerSize=2)
#' 
fastqPairedFilter <- function(fn, fout, maxN = c(0,0), truncQ = c(2,2), truncLen = c(0,0), maxLen=c(Inf, Inf), minLen=c(20, 20), trimLeft = c(0,0), trimRight=c(0,0), minQ = c(0,0), maxEE = c(Inf, Inf), rm.phix = c(TRUE, TRUE), rm.lowcomplex = c(0, 0), matchIDs = FALSE, orient.fwd=NULL, id.sep = "\\s", id.field = NULL, n = 1e6, OMP=TRUE, qualityType="Auto", compress = TRUE, verbose = FALSE, ...){
  if(!OMP) {
    ompthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
    on.exit(.Call(ShortRead:::.set_omp_threads, ompthreads))
  }
  # Warning: This assumes that forward/reverse reads are in the same order unless matchIDs=TRUE
  if(!is.character(fn) || length(fn) != 2) stop("Two paired input file names required.")
  if(!is.character(fout) || length(fout) != 2) stop("Two paired output file names required.")
  
  if(any(duplicated(c(fn, fout)))) { stop("The output and input file names must be different.") }
  
  for(var in c("maxN", "truncQ", "truncLen", "maxLen", "minLen", "trimLeft", "trimRight", "minQ", "maxEE", "rm.phix", "rm.lowcomplex")) {
    if(length(get(var)) == 1) { # Double the 1 value to be the same for F and R
      assign(var, c(get(var), get(var)))
    }
    if(length(get(var)) != 2) {
      stop(paste("Input variable", var, "must be length 1 or 2 (Forward, Reverse)."))
    }
  }
  
  startF <- max(1, trimLeft[[1]] + 1, na.rm=TRUE)
  startR <- max(1, trimLeft[[2]] + 1, na.rm=TRUE)
  
  endF <- truncLen[[1]]
  if(endF < startF) { endF = NA }
  endF <- endF - startF + 1
  endR <- truncLen[[2]]
  if(endR < startR) { endR = NA }
  endR <- endR - startR + 1
  
  ## iterating over forward and reverse files using fastq streaming
  fF <- FastqStreamer(fn[[1]], n = n)
  on.exit(close(fF))
  fR <- FastqStreamer(fn[[2]], n = n)
  on.exit(close(fR), add=TRUE)
  
  if(file.exists(fout[[1]])) {
    if(file.remove(fout[[1]])) {
      if(verbose) message("Overwriting file:", fout[[1]])
    } else {
      stop("Failed to overwrite file:", fout[[1]])
    }
  }
  if(file.exists(fout[[2]])) {
    if(file.remove(fout[[2]])) {
      if(verbose) message("Overwriting file:", fout[[2]])
    } else {
      stop("Failed to overwrite file:", fout[[2]])
    }
  }
  
  first=TRUE
  remainderF <- ShortReadQ(); remainderR <- ShortReadQ()
  casava <- "Undetermined"
  inseqs = 0; outseqs = 0
  while( TRUE ) {
    suppressWarnings(fqF <- yield(fF, qualityType = qualityType))
    suppressWarnings(fqR <- yield(fR, qualityType = qualityType))
    if(length(fqF) == 0 && length(fqR) == 0) { break } # Loop Logic
    
    inseqs <- inseqs + length(fqF)
    
    if(matchIDs) {
      if(first) { 
        if(is.null(id.field)) {
          # Determine the sequence identifier field. Looks for a single 6-colon field (CASAVA 1.8+ id format)
          # or a single 4-colon field (earlier format). Fails if it doesn't find such a field.
          id1 <- as.character(id(fqF)[[1]])
          id.fields <- strsplit(id1, id.sep)[[1]]
          ncolon <- sapply(gregexpr(":", id.fields), length)
          ncoltab <- table(ncolon)
          if(max(ncolon) == 6 && ncoltab["6"] == 1) { # CASAVA 1.8+ format
            casava <- "Current"
            id.field <- which(ncolon == 6)
          } else if (max(ncolon) == 4 && ncoltab["4"] == 1) { # CASAVA <=1.7 format
            casava <- "Old"
            id.field <- which(ncolon == 4)
          } else { # Couldn't unambiguously find the seq id field
            stop("Couldn't automatically detect the sequence identifier field in the fastq id string.")
          }
        }
      } else { # !first
        # Prepend the unmatched sequences from the end of previous chunks
        # Need ShortRead::append or the method is not dispatched properly
        fqF <- append(remainderF, fqF)
        fqR <- append(remainderR, fqR)
      }
    } else { # !matchIDs
      if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files: ", length(fqF), ", ", length(fqR), ".")
    }
    
    # Enforce id matching (ASSUMES SAME ORDERING IN F/R, BUT ALLOWS DIFFERENTIAL MEMBERSHIP)
    # Keep the tail of unmatched sequences (could match next chunk)
    if(matchIDs) {
      idsF <- sapply(strsplit(as.character(id(fqF)), id.sep), `[`, id.field)
      idsR <- sapply(strsplit(as.character(id(fqR)), id.sep), `[`, id.field)
      if(casava == "Old") { # Drop the index number/pair identifier (i.e. 1=F, 2=R)
        idsF <- sapply(strsplit(idsF, "#"), `[`, 1)
      }
      lastF <- max(c(0,which(idsF %in% idsR)))
      lastR <- max(c(0,which(idsR %in% idsF)))
      if(lastF < length(fqF)) {
        remainderF <- fqF[(lastF+1):length(fqF)]
      } else {
        remainderF <- ShortReadQ() 
      }
      if(lastR < length(fqR)) {
        remainderR <- fqR[(lastR+1):length(fqR)]
      } else {
        remainderR <- ShortReadQ() 
      }
      fqF <- fqF[idsF %in% idsR]
      fqR <- fqR[idsR %in% idsF]
    }
    
    # Enforce orient.fwd
    if(!is.null(orient.fwd)) {
      if(!C_isACGT(orient.fwd)) stop("Non-ACGT characters detected in orient.fwd")
      barlen <- nchar(orient.fwd)
      keepF <- narrow(sread(fqF),1,barlen) == orient.fwd
      keepR <- (narrow(sread(fqR),1,barlen) == orient.fwd) & !keepF
      fq <- ShortReadQ(sread=c(sread(fqF[keepF]), sread(fqR[keepR])), 
                       quality=c(quality(quality(fqF[keepF])), quality(quality(fqR[keepR]))), 
                       id=c(id(fqF[keepF]), id(fqR[keepR])))
      fqR <- ShortReadQ(sread=c(sread(fqR[keepF]), sread(fqF[keepR])), 
                       quality=c(quality(quality(fqR[keepF])), quality(quality(fqF[keepR]))), 
                       id=c(id(fqR[keepF]), id(fqF[keepR])))
      fqF <- fq
      rm(fq)
    }
    # Enforce maxLen
    if(is.finite(maxLen[[1]]) || is.finite(maxLen[[2]])) {
      keep <- width(fqF) <= maxLen[[1]] & width(fqR) <= maxLen[[2]]
      fqF <- fqF[keep]
      fqR <- fqR[keep]
    }
    # Trim left
    keep <- (width(fqF) >= startF & width(fqR) >= startR)
    fqF <- fqF[keep]
    fqF <- narrow(fqF, start = startF, end = NA)
    fqR <- fqR[keep]
    fqR <- narrow(fqR, start = startR, end = NA)
    # Trim right
    if(trimRight[[1]] > 0) {
      keep <- width(fqF) > trimRight[[1]]
      fqF <- fqF[keep]; fqR <- fqR[keep]
      fqF <- narrow(fqF, start=NA, end=width(fqF)-trimRight[[1]])
    }
    if(trimRight[[2]] > 0) {
      keep <- width(fqR) > trimRight[[2]]
      fqF <- fqF[keep]; fqR <- fqR[keep]
      fqR <- narrow(fqR, start=NA, end=width(fqR)-trimRight[[2]])
    }
    # Trim on truncQ
    # Convert numeric quality score to the corresponding ascii character
    encF <- encoding(quality(fqF))
    encR <- encoding(quality(fqR))
    if(is.numeric(truncQ)) {
      indF <- which(encF==truncQ[[1]])
      indR <- which(encR==truncQ[[2]])
      if(!(length(indF) == 1 && length(indR) == 1)) stop("Encoding for this truncQ value not found.")
      truncQ <- c(names(encF)[[indF]], names(encR)[[indR]])
    }
    if(length(fqF) > 0) {
      rngF <- trimTails(fqF, 1, truncQ[[1]], ranges=TRUE)
      fqF <- narrow(fqF, 1, end(rngF)) # have to do it this way to avoid dropping the zero lengths
    }
    if(length(fqR) > 0) {
      rngR <- trimTails(fqR, 1, truncQ[[2]], ranges=TRUE)
      fqR <- narrow(fqR, 1, end(rngR)) # have to do it this way to avoid dropping the zero lengths
    }
    truncQ <- c(encF[truncQ[1]], encR[truncQ[2]]) # convert back to integer
    # And now filter any with length zero in F or R
    # May want to roll this into the next length cull step...
    keep <- (width(fqF) > 0 & width(fqR) > 0)
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    
    # Filter any with less than required length
    keep <- rep(TRUE, length(fqF))
    if(!is.na(endF)) { keep <- keep & (width(fqF) >= endF) }
    if(!is.na(endR)) { keep <- keep & (width(fqR) >= endR) }
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    # Truncate to truncLen
    fqF <- narrow(fqF, start = 1, end = endF)
    fqR <- narrow(fqR, start = 1, end = endR)
    # Enforce minLen
    keep <- width(fqF) >= minLen[[1]] & width(fqR) >= minLen[[2]]
    fqF <- fqF[keep]
    fqR <- fqR[keep]

    # Filter based on minQ and Ns and maxEE
    keep <- .nFilter(fqF, maxN[[1]]) & .nFilter(fqR, maxN[[2]])
    fqF <- fqF[keep]; fqR <- fqR[keep]
    keep <- rep(TRUE, length(fqF))
    qmat <- as(quality(fqF), "matrix")
    if(minQ[[1]] > truncQ[[1]]) suppressWarnings(keep <- keep & (apply(qmat, 1, min, na.rm=TRUE)>minQ[[1]]))
    if(maxEE[[1]] < Inf) keep <- keep & C_matrixEE(qmat) <= maxEE[[1]]
    qmat <- as(quality(fqR), "matrix")
    if(minQ[[2]] > truncQ[[2]]) suppressWarnings(keep <- keep & (apply(qmat, 1, min, na.rm=TRUE)>minQ[[2]]))
    if(maxEE[[2]] < Inf) keep <- keep & C_matrixEE(qmat) <= maxEE[[2]]
    fqF <- fqF[keep]; fqR <- fqR[keep]
    rm(qmat)
    
    if(length(fqF) != length(fqR)) stop("Filtering caused mismatch between forward and reverse sequence lists: ", length(fqF), ", ", length(fqR), ".")
    
    # Remove phiX
    if(rm.phix[[1]] && rm.phix[[2]]) {
      is.phi <- isPhiX(as(sread(fqF), "character"), ...)
      is.phi <- is.phi | isPhiX(as(sread(fqR), "character"), ...)
    } else if(rm.phix[[1]] && !rm.phix[[2]]) {
      is.phi <- isPhiX(as(sread(fqF), "character"), ...)
    } else if(!rm.phix[[1]] && rm.phix[[2]]) {
      is.phi <- isPhiX(as(sread(fqR), "character"), ...)
    }
    if(any(rm.phix)) {
      fqF <- fqF[!is.phi]
      fqR <- fqR[!is.phi]
    }
    
    # Remove low complexity
    if(rm.lowcomplex[[1]] && rm.lowcomplex[[2]]) {
      is.lowc <- (seqComplexity(sread(fqF), ...) < rm.lowcomplex[[1]])
      is.lowc <- is.lowc | (seqComplexity(sread(fqR), ...) < rm.lowcomplex[[2]])
    } else if(rm.lowcomplex[[1]] && !rm.lowcomplex[[2]]) {
      is.lowc <- (seqComplexity(sread(fqF), ...) < rm.lowcomplex[[1]])
    } else if(!rm.lowcomplex[[1]] && rm.lowcomplex[[2]]) {
      is.lowc <- (seqComplexity(sread(fqR), ...) < rm.lowcomplex[[2]])
    }
    if(rm.lowcomplex[[1]] || rm.lowcomplex[[2]]) {
      fqF <- fqF[!is.lowc]
      fqR <- fqR[!is.lowc]
    }

    outseqs <- outseqs + length(fqF)    
    
    if(first) {
      writeFastq(fqF, fout[[1]], "w", compress = compress)
      writeFastq(fqR, fout[[2]], "w", compress = compress)
      first=FALSE
    } else {
      writeFastq(fqF, fout[[1]], "a", compress = compress)
      writeFastq(fqR, fout[[2]], "a", compress = compress)
    }
  }
  
  if(outseqs==0) {
  }
  
  if(verbose) {
    outperc <- round(outseqs * 100 / inseqs, 1)
    outperc <- paste(" (", outperc, "%)", sep="")
    message("Read in ", inseqs, " paired-sequences, output ", outseqs, outperc, " filtered paired-sequences.", sep="")
  }
  
  if(outseqs==0) {
    message(paste("The filter removed all reads:", fout[[1]], "and", fout[[2]], "not written."))
    file.remove(fout[[1]])
    file.remove(fout[[2]])
  }
  
  return(invisible(c(reads.in=inseqs, reads.out=outseqs)))
}
################################################################################
#' Determine if input sequence(s) match the phiX genome.
#' 
#' This function compares the word-profile of the input
#' sequences to the phiX genome, and the reverse complement of the phiX genome. If
#' enough exactly matching words are found, the sequence is flagged.
#' 
#' @param seqs (Required). A \code{character} vector of A/C/G/T sequences.
#' 
#' @param wordSize (Optional). Default 16.
#'  The size of the words to use for comparison.
#'   
#' @param minMatches (Optional). Default 2.
#' The minimum number of words in the input sequences that must match the phiX genome
#' (or its reverse complement) for the sequence to be flagged.
#' 
#' @param nonOverlapping (Optional). Default TRUE.
#' If TRUE, only non-overlapping matching words are counted.
#' 
#' @param ... (Optional). Ignored.
#' 
#' @return \code{logical(1)}.
#'  TRUE if sequence matched the phiX genome.
#'
#' @seealso 
#'  \code{\link{fastqFilter}}, \code{\link{fastqPairedFilter}}
#'  
#' @export
#' 
#' @importFrom ShortRead readFasta
#' @importFrom methods as
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' sqs1 <- getSequences(derep1)
#' is.phi <- isPhiX(sqs1)
#' is.phi <- isPhiX(sqs1, wordSize=20,  minMatches=1)
#' 
isPhiX <- function(seqs, wordSize=16, minMatches=2, nonOverlapping=TRUE, ...) {
  seqs <- getSequences(seqs)
  sq.phix <- as(sread(readFasta(system.file("extdata", "phix_genome.fa", package="dada2"))), "character")
  rc.phix <- rc(sq.phix)
  hits <- C_matchRef(seqs, sq.phix, wordSize, nonOverlapping)
  hits.rc <- C_matchRef(seqs, rc.phix, wordSize, nonOverlapping)
  return((hits >= minMatches) | (hits.rc >= minMatches))
}
################################################################################
#' Determine if input sequence(s) are low complexity.
#'
#' This function calculates the kmer
#' complexity of input sequences. Complexity is quantified as the Shannon
#' richness of kmers, which can be thought of as the
#' effective number of kmers if they were all
#' at equal frequencies. If a window size is provided, the minimum Shannon
#' richness observed over sliding window along the sequence is returned.
#'
#' This function can be used to identify potentially artefactual or undesirable
#' low-complexity sequences, or sequences with low-complexity regions, as are
#' sometimes observed in Illumina sequencing runs. When such artefactual
#' sequences are present, the Shannon kmer
#' richness values returned by this function will typically show a clear
#' bimodal signal.
#' 
#' Kmers with non-ACGT characters are ignored. Also note that no correction is
#' performed for sequence lengths. This is important when using longer kmer
#' lengths, where 4^wordSize approaches the length of the sequence, as shorter
#' sequences will then have a lower effective richness simply due to their
#' being too little sequence to sample all the possible kmers.
#'
#' @param seqs (Required). A \code{character} vector of A/C/G/T sequences, or
#' any object coercible by \code{\link{getSequences}}.
#'
#' @param kmerSize (Optional). Default 2.
#'  The size of the kmers (or "oligonucleotides" or "words") to use.
#'
#' @param window (Optional). Default NULL.
#' The width in nucleotides of the moving window. If NULL the whole sequence is used.
#'
#' @param by (Optional). Default 5.
#' The step size in nucleotides between each moving window tested.
#'
#' @param ... (Optional). Ignored.
#' 
#' @return \code{numeric}.
#'  A vector of minimum kmer complexities for each sequence.
#'
#' @seealso
#'  \code{\link{plotComplexity}}
#'  \code{\link[Biostrings]{oligonucleotideFrequency}}
#'
#' @export
#'
#' @importFrom Biostrings oligonucleotideFrequency
#' @importFrom Biostrings DNAStringSet
#' @importFrom IRanges narrow
#' @importFrom Biostrings width
#' @importFrom methods is
#'
#' @examples
#' sq.norm <- "TACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGA"
#' sq.lowc <- "TCCTTCTTCTCCTCTCTTTCTCCTTCTTTCTTTTTTTTCCCTTTCTCTTCTTCTTTTTCTTCCTTCCTTTTTTC"
#' sq.part <- "TTTTTCTTCTCCCCCTTCCCCTTTCCTTTTCTCCTTTTTTCCTTTAGTGCAGTTGAGGCAGGCGGAATTCGTGG"
#' sqs <- c(sq.norm, sq.lowc, sq.part)
#' seqComplexity(sqs)
#' seqComplexity(sqs, kmerSize=3, window=25)
#'
seqComplexity <- function(seqs, kmerSize=2, window=NULL, by=5, ...) {
  if(!is.null(window) && kmerSize >= window) stop("The window over which kmer frequency is calculated must be larger than the kmerSize.")
  if(!is(seqs, "DNAStringSet")) {
    seqs <- getSequences(seqs)
    if(!(all(C_isACGT(seqs)))) warning("Not all sequences were A/C/G/T only, which can interfere with the calculation of the Shannon richness.")
    seqs <- DNAStringSet(seqs)
  }
  si.max <- 4**kmerSize
  if(is.null(window)) {
    si <- apply(oligonucleotideFrequency(seqs, width=kmerSize), 1, sindex)
  } else {
    si <- rep(si.max, length(seqs))
    for(i in seq(1, max(width(seqs))-window, by=by)) {
      keep <- (width(seqs) >= (i+window-1))
      wind <- narrow(seqs[keep],i,i+window-1)
      si.i <- apply(oligonucleotideFrequency(wind, width=kmerSize), 1, sindex)
      si[keep] <- pmin(si[keep], si.i)
    }
  }
  si
}
## Helper function to calculate the effective shannon richness
## Which is basically the exponential of the usual Shannon diversity
sindex <- function(x) {
  y <- x/sum(x)
  y <- y[y>0]
  exp(sum(-y*log(y)))
}

## Helper function to calculate number of Ns in fastq file
## Near drop-in replacement for ShortRead::nFilter (and uses that functions internal code)
##   but avoids the srFilter step within that function that is causing issues with different Matrix package versions
## Note important differences in return format (logical instead of SRFilterResult) and argument format (x, max.n) from ShortRead::nread
## But functionally should act identical, and does in testing. E.g.:
## fq <- readFastq(fn[[1]]); maxN <- 1
## keep.old <- ShortRead::nFilter(maxN)(fq) ### Note this errors with Matrix 1.3.3
## keep.new <- dada2:::.nFilter(fq, maxN)
## identical(fq[keep.old], fq[keep.new]) ## TRUE
#'
#' @importFrom ShortRead sread
#' @importFrom Biostrings alphabetFrequency
#' @importFrom methods is
#'
.nFilter <- function(x, max.n=0) {
  if (is(x, "ShortRead")) 
    alphabetFrequency(sread(x), baseOnly=TRUE)[,"other"] <= max.n
  else alphabetFrequency(x, baseOnly=TRUE)[,"other"] <= max.n
}
