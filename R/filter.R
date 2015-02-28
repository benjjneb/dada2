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

