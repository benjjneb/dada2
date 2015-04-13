# TODO: Add EE or AveQ filtering.
# TODO: Added paired filtering.

#' filterFastq implements most of fastq_filter from usearch...
#' 
#' @export
fastqFilter <- function(fn, fout, maxN = 0, truncQ = "#", truncLen = 0, trimLeft = 0, minQ = 0, n = 1e6, verbose = FALSE){
  # See also filterFastq in the ShortRead package
  if(trimLeft <= 0) { 
    start = 1 
  } else { 
    start = trimLeft+1
  } 
  
  if(truncLen <= 0 || truncLen < start) { 
    end = NA
  } else { 
    end = truncLen 
  }
  
  require("ShortRead")
  ## iterating over an entire file using fastq streaming
  f <- FastqStreamer(fn, n = n)
  on.exit(close(f))
  
  first=TRUE
  inseqs = 0
  outseqs = 0
  while( length(suppressWarnings(fq <- yield(f))) ){
    inseqs <- inseqs + length(fq)
    # Trim on truncQ ################## BEST TO MAKE THIS ROBUST TO Q ALPHABET
    fq <- trimTails(fq, 1, truncQ)
    # Filter any with less than required length
    if(!is.na(end)) { fq <- fq[width(fq) >= end] }
    # Trim left and truncate to truncLen
    fq <- narrow(fq, start = start, end = end)
    # Filter based on minQ
    fq <- fq[minQFilter(minQ)(fq)]
    # Filter Ns
    fq <- fq[nFilter(maxN)(fq)]
    
    outseqs <- outseqs + length(fq)
    
    if(first) {
      writeFastq(fq, fout, "w")
      first=FALSE
    } else {
      writeFastq(fq, fout, "a")
    }
  }
  
  if(verbose) {
    cat("Read in", inseqs, "sequences, outputted", outseqs, "filtered sequences.")
  }
}

#' This filters paired read files, keeping only those reads which pass in both forward and reverse files
#' 
#' @export
fastqPairedFilter <- function(fn, fout, maxN = c(0,0), truncQ = c("#","#"), truncLen = c(0,0), trimLeft = c(0,0), minQ = c(0,0), n = 1e6, compress = TRUE, verbose = FALSE){
  # Warning: This assumes that forward/reverse reads are paired by line
  # IT DOES NOT CHECK THE ID LINES
  if(!is.character(fn) || length(fn) != 2) stop("Two paired input file names required.")
  if(!is.character(fout) || length(fout) != 2) stop("Two paired output file names required.")
  
  for(var in c("maxN", "truncQ", "truncLen", "trimLeft", "minQ")) {
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
  if(endR < startR) { endF = NA }
  
  require("ShortRead")
  ## iterating over forward and reverse files using fastq streaming
  fF <- FastqStreamer(fn[[1]], n = n)
  on.exit(close(fF))
  fR <- FastqStreamer(fn[[2]], n = n)
  on.exit(close(fR), add=TRUE)
  
  first=TRUE
  inseqs = 0
  outseqs = 0
  while( length(suppressWarnings(fqF <- yield(fF))) && length(suppressWarnings(fqR <- yield(fR))) ){
    if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files(1): ", length(fqF), ", ", length(fqR), ".")
    inseqs <- inseqs + length(fqF)
    # Trim on truncQ ################## BEST TO MAKE THIS ROBUST TO Q ALPHABET
    rngF <- trimTails(fqF, 1, truncQ[[1]], ranges=TRUE)
    fqF <- narrow(fqF, 1, end(rngF)) # have to do it this way to avoid dropping the zero lengths
    rngR <- trimTails(fqR, 1, truncQ[[2]], ranges=TRUE)
    fqR <- narrow(fqR, 1, end(rngR)) # have to do it this way to avoid dropping the zero lengths
    if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files(2): ", length(fqF), ", ", length(fqR), ".")
    # Filter any with less than required length
    keep <- rep(TRUE, length(fqF))
    if(!is.na(endF)) { keep <- keep & (width(fqF) >= endF) }
    if(!is.na(endR)) { keep <- keep & (width(fqR) >= endR) }
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    # Trim left and truncate to truncLen
    fqF <- narrow(fqF, start = startF, end = endF)
    fqR <- narrow(fqR, start = startR, end = endR)
    # Filter based on minQ and Ns
    keep <- rep(TRUE, length(fqF))
    keep <- keep & minQFilter(minQ[[1]])(fqF) & nFilter(maxN[[1]])(fqF)
    keep <- keep & minQFilter(minQ[[2]])(fqR) & nFilter(maxN[[2]])(fqF)
    fqF <- fqF[keep]
    fqR <- fqR[keep]
    
    outseqs <- outseqs + length(fqF)
    if(length(fqF) != length(fqR)) stop("Mismatched forward and reverse sequence files(3): ", length(fqF), ", ", length(fqR), ".")
    
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
    cat("Read in", inseqs, "sequences, outputted", outseqs, "filtered sequences.")
  }
}


minQFilter <- function (minQ = 0L, .name = "MinQFilter") 
{
  ShortRead:::.check_type_and_length(minQ, "numeric", 1)
  ShortRead:::srFilter(function(x) {
    apply(as(quality(x), "matrix"), 1, function(qs) min(qs) >= minQ)
  }, name = .name)
}

minLenFilter <- function (minLen = 0L, .name = "MinLenFilter") 
{
  ShortRead:::.check_type_and_length(minLen, "numeric", 1)
  ShortRead:::srFilter(function(x) {
    width(x) >= minLen
  }, name = .name)
}

