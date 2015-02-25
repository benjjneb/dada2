#' filterFastq implements most of fastq_filter from usearch...
#' 
#' @export
filterFastq <- function(fn, fout, maxN = 0, truncQ = "#", truncLen = 0, trimLeft = 0, n = 1e6, verbose = FALSE){
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
  while( length(suppressWarnings(fq <- yield(f))) ){
    # Trim on truncQ ################## BEST TO MAKE THIS ROBUST TO Q ALPHABET
    fq <- trimTails(fq, 1, truncQ)
    # Filter any with less than required length
    if(!is.na(end)) { fq <- fq[width(fq) >= end] }
    # Trim left and truncate to truncLen
    fq <- narrow(fq, start = start, end = end)
    # Filter Ns
    fq <- fq[nFilter(maxN)(fq)]
    
    if(first) {
      writeFastq(fq, fout, "w")
      first=FALSE
    } else {
      writeFastq(fq, fout, "a")
    }
  }
}


minQFilter <- function (minQ = 0L, .name = "MinQFilter") 
{
  ShortRead:::.check_type_and_length(minQ, "numeric", 1)
  ShortRead:::srFilter(function(x) {
    apply(as(quality(x), "matrix"), 1, function(qs) {
      min(qs) > minQ
    })
  }, name = .name)
}

minLenFilter <- function (minLen = 0L, .name = "MinLenFilter") 
{
  ShortRead:::.check_type_and_length(minLen, "numeric", 1)
  ShortRead:::srFilter(function(x) {
    width(x) >= minLen
  }, name = .name)
}

