#' fastqFilter filters and trims fastq files.
#' 
#' fastqFilter takes an input fastq file (can be compressed), filters it based on several
#' user-definable criteria, and outputs those reads which pass the filter and their associated
#' qualities to a new fastq file (also can be compressed). Several functions in the ShortRead
#' package are leveraged to do this filtering.
#' 
#' fastqFilter replicates most of the functionality of the fastq_filter command in usearch
#' (http://www.drive5.com/usearch/manual/cmd_fastq_filter.html). 
#' 
#' @param fn (Required). A character string naming the path to the fastq file, or an R connection.
#'   
#' @param fout (Required). A character string naming the path to the output file, or an R connection.
#' 
#' @param truncQ (Optional). Truncate reads at the first instance of a quality score less than
#'    or equal to truncQ. Default is 2, a special quality score indicating the end of good quality
#'    sequence in Illumina 1.8+. Can provide truncQ as an integer or appropriate ascii encoding. 
#'  
#' @param truncLen (Optional). A \code{numeric(1)} Truncate after truncLen bases, reads shorter than
#'    this are discarded.
#'  
#' @param trimLeft (Optional). Remove trimLeft nucleotides from the start of each read. If both
#'    truncLen and trimLeft are used, all filtered reads will have length truncLen-trimLeft.
#'  
#' @param maxN (Optional). After truncation, sequences with more than maxN Ns will be discarded.
#'    Default is 0. Currently dada() does not allow Ns.
#'  
#' @param minQ (Optional). After truncation, reads contain a quality score below minQ will be discarded.
#'
#' @param maxEE (Optional). After truncation, reads with higher than maxEE "expected errors" will be discarded.
#'  Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
#'  
#' @param n (Optional). The number of records (reads) to read in and filter at any one time. 
#'  This controls the peak memory requirement so that very large fastq files are supported. 
#'  Default is \code{1e6}, one-million reads. See \code{\link{FastqStreamer}} for details.
#'
#' @param compress (Optional). A \code{logical(1)} indicating whether the output should be gz compressed.
#' 
#' @param verbose (Optional). A \code{logical(1)}. If TRUE, some status messages are displayed.  
#'     
#' @seealso 
#'  \code{\link{fastqPairedFilter}}
#' 
#'  \code{\link[ShortRead]{FastqStreamer}}
#'  
#'  \code{\link[ShortRead]{srFilter}}
#'  
#'  \code{\link[ShortRead]{trimTails}}
#' 
#' @export
#' 
#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead yield
#' @importFrom ShortRead writeFastq
#' @importFrom ShortRead trimTails
#' @importFrom ShortRead nFilter
#' @importFrom Biostrings quality
#' @importFrom Biostrings narrow
#' @importFrom Biostrings width
#' 
fastqFilter <- function(fn, fout, truncQ = 2, truncLen = 0, trimLeft = 0, maxN = 0, minQ = 0, maxEE = Inf, n = 1e6, compress = TRUE, verbose = FALSE){
  # See also filterFastq in the ShortRead package
  start <- max(1, trimLeft + 1)
  end <- truncLen
  if(end < start) { end = NA }
  
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
  while( length(suppressWarnings(fq <- yield(f))) ){
    inseqs <- inseqs + length(fq)
    
    # Trim on truncQ 
    # Convert numeric quality score to the corresponding ascii character
    if(is.numeric(truncQ)) {
      enc <- encoding(quality(fq))
      ind <- which(enc==truncQ)
      if(length(ind) != 1) stop("Encoding for this truncQ value not found.")
      truncQ <- names(enc)[[ind]]
    }
    fq <- trimTails(fq, 1, truncQ)
    
    # Filter any with less than required length
    if(!is.na(end)) { fq <- fq[width(fq) >= end] }
    # Trim left and truncate to truncLen
    fq <- narrow(fq, start = start, end = end)
    
    # Filter based on minQ and Ns and maxEE
    keep <- rep(TRUE, length(fq))
    keep <- keep & minQFilter(minQ)(fq)
    keep <- keep & nFilter(maxN)(fq)
    if(maxEE < Inf) keep <- keep & maxEEFilter(maxEE)(fq)
    fq <- fq[keep]
    
    outseqs <- outseqs + length(fq)
    
    if(first) {
      writeFastq(fq, fout, "w", compress=compress)
      first=FALSE
    } else {
      writeFastq(fq, fout, "a", compress=compress)
    }
  }
  
  if(verbose) {
    cat("Read in", inseqs, "sequences, outputted", outseqs, "filtered sequences.")
  }
}

#' fastqPairedFilter filters and trims paired forward and reverse fastq files.
#' 
#' fastqPairedFilter takes in two input fastq file (can be compressed), filters them based on several
#' user-definable criteria, and outputs those reads which pass the filter in both directions along
#' with their associated qualities to two new fastq file (also can be compressed). Several functions
#' in the ShortRead package are leveraged to do this filtering. The filtered forward/reverse reads
#' remain identically ordered.
#' 
#' fastqPairedFilter replicates most of the functionality of the fastq_filter command in usearch
#' (http://www.drive5.com/usearch/manual/cmd_fastq_filter.html) but applied to forward and reverse
#' reads simultaneously.
#' 
#' @param fn (Required). A \code{character(2)} naming the paths to the forward/reverse fastq files.
#'   
#' @param fout (Required). A \code{character(2)} naming the path to the output file.
#' 
#' FURTHER ARGUMENTS can be provided as a length 1 or length 2 vector. If provided as a length 1
#' vector the same criteria is used for forward and reverse. If provided as a length 2 vector, the
#' first value is used for the forward reads, the second value for the reverse reads.
#' 
#' @param truncQ (Optional). Truncate reads at the first instance of a quality score less than
#'    or equal to truncQ. Default is 2, a special quality score indicating the end of good quality
#'    sequence in Illumina 1.8+. Can provide truncQ as an integer or appropriate ascii encoding. 
#'  
#' @param truncLen (Optional). A \code{numeric(1)} Truncate after truncLen bases, reads shorter than
#'    this are discarded. THIS PARAMETER TRIMS ALL READS TO THE SAME LENGTH WHICH IS NEEDED FOR THE
#'    dada() FUNCTION!
#'  
#' @param trimLeft (Optional). Remove trimLeft nucleotides from the start of each read. If both
#'    truncLen and trimLeft are used, all filtered reads will have length truncLen-trimLeft.
#'  
#' @param maxN (Optional). After truncation, sequences with more than maxN Ns will be discarded.
#'    Default is 0. Currently dada() does not allow Ns.
#'  
#' @param minQ (Optional). After truncation, reads contain a quality score below minQ will be discarded.
#'
#' @param maxEE (Optional). After truncation, reads with higher than maxEE "expected errors" will be discarded.
#'  Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))
#'  
#' @param n (Optional). The number of records (reads) to read in and filter at any one time. 
#'  This controls the peak memory requirement so that very large fastq files are supported. 
#'  Default is \code{1e6}, one-million reads. See \code{\link{FastqStreamer}} for details.
#'
#' @param compress (Optional). A \code{logical(1)} indicating whether the output should be gz compressed.
#' 
#' @param verbose (Optional). A \code{logical(1)}. If TRUE, some status messages are displayed.
#'    
#' @seealso 
#'  \code{\link{fastqFilter}}
#' 
#'  \code{\link[ShortRead]{FastqStreamer}}
#'  
#'  \code{\link[ShortRead]{srFilter}}
#'  
#'  \code{\link[ShortRead]{trimTails}}
#' 
#' @export
#' 
#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead yield
#' @importFrom ShortRead writeFastq
#' @importFrom ShortRead trimTails
#' @importFrom ShortRead nFilter
#' @importFrom Biostrings quality
#' @importFrom Biostrings narrow
#' @importFrom Biostrings width
#' 
fastqPairedFilter <- function(fn, fout, maxN = c(0,0), truncQ = c(2,2), truncLen = c(0,0), trimLeft = c(0,0), minQ = c(0,0), maxEE = c(Inf, Inf), n = 1e6, compress = TRUE, verbose = FALSE){
  # Warning: This assumes that forward/reverse reads are in the same order
  # IT DOES NOT CHECK THE ID LINES
  if(!is.character(fn) || length(fn) != 2) stop("Two paired input file names required.")
  if(!is.character(fout) || length(fout) != 2) stop("Two paired output file names required.")
  
  for(var in c("maxN", "truncQ", "truncLen", "trimLeft", "minQ", "maxEE")) {
    if(length(get(var)) == 1) { # Double the 1 value to be the same for F and R
      assign(var, c(get(var), get(var)))
    }
    if(length(get(var)) != 2) {
      stop(paste("Input variable", var, "must be length 2 (Forward, Reverse)."))
    }
  }

  startF <- max(1, trimLeft[[1]] + 1)
  startR <- max(1, trimLeft[[2]] + 1)

  endF <- truncLen[[1]]
  if(endF < startF) { endF = NA }
  endR <- truncLen[[2]]
  if(endR < startR) { endR = NA }
  
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
  inseqs = 0
  outseqs = 0
  while( length(suppressWarnings(fqF <- yield(fF))) && length(suppressWarnings(fqR <- yield(fR))) ){
    if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files: ", length(fqF), ", ", length(fqR), ".")
    inseqs <- inseqs + length(fqF)
    
    # Trim on truncQ
    # Convert numeric quality score to the corresponding ascii character
    if(is.numeric(truncQ)) {
      encF <- encoding(quality(fqF))
      encR <- encoding(quality(fqR))
      indF <- which(encF==truncQ[[1]])
      indR <- which(encR==truncQ[[2]])
      if(!(length(indF) == 1 && length(indR) == 1)) stop("Encoding for this truncQ value not found.")
      truncQ <- c(names(encF)[[indF]], names(encR)[[indR]])
    }
    rngF <- trimTails(fqF, 1, truncQ[[1]], ranges=TRUE)
    fqF <- narrow(fqF, 1, end(rngF)) # have to do it this way to avoid dropping the zero lengths
    rngR <- trimTails(fqR, 1, truncQ[[2]], ranges=TRUE)
    fqR <- narrow(fqR, 1, end(rngR)) # have to do it this way to avoid dropping the zero lengths
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
    # Trim left and truncate to truncLen
    fqF <- narrow(fqF, start = startF, end = endF)
    fqR <- narrow(fqR, start = startR, end = endR)
    
    # Filter based on minQ and Ns and maxEE
    keep <- rep(TRUE, length(fqF))
    keep <- keep & minQFilter(minQ[[1]])(fqF) & nFilter(maxN[[1]])(fqF)
    if(maxEE[[1]] < Inf) keep <- keep & maxEEFilter(maxEE[[1]])(fqF)
    keep <- keep & minQFilter(minQ[[2]])(fqR) & nFilter(maxN[[2]])(fqR)
    if(maxEE[[2]] < Inf) keep <- keep & maxEEFilter(maxEE[[2]])(fqR)
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    
    if(length(fqF) != length(fqR)) stop("Filtering caused mismatch between forward and reverse sequence lists: ", length(fqF), ", ", length(fqR), ".")
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
  
  if(verbose) {
    cat("Read in", inseqs, "paired-sequences, outputted", outseqs, "filtered paired-sequences.\n")
  }
}

# # @importFrom ShortRead :::.check_type_and_length
#' @importFrom ShortRead srFilter
#' @importFrom Biostrings quality
minQFilter <- function (minQ = 0L, .name = "MinQFilter") 
{
  ShortRead:::.check_type_and_length(minQ, "numeric", 1)
  srFilter(function(x) {
    apply(as(quality(x), "matrix"), 1, function(qs) min(qs, na.rm=TRUE) >= minQ)
  }, name = .name)
}

# # @importFrom ShortRead :::.check_type_and_length
#' @importFrom Biostrings quality
#' @importFrom ShortRead srFilter
maxEEFilter <- function (maxEE = Inf, .name = "MaxEEFilter") 
{
  ShortRead:::.check_type_and_length(maxEE, "numeric", 1)
  srFilter(function(x) {
    apply(as(quality(x), "matrix"), 1, function(qs) sum(10^(-qs[!is.na(qs)]/10.0)) <= maxEE)
  }, name = .name)
}

# # @importFrom ShortRead :::.check_type_and_length
#' @importFrom Biostrings width
minLenFilter <- function(minLen = 0L, .name = "MinLenFilter"){
  ShortRead:::.check_type_and_length(minLen, "numeric", 1)
  srFilter(function(x) {
    width(x) >= minLen
  }, name = .name)
}
