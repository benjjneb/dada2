
dada_opts <- new.env()
assign("OMEGA_A", 0.01, envir = dada_opts) # NOT YET PASSED IN
assign("USE_KMERS", TRUE, envir = dada_opts)
assign("KDIST_CUTOFF", 0.5, envir = dada_opts)
assign("BAND_SIZE", 50, envir = dada_opts) # NOT YET PASSED IN
assign("MAX_CONSIST", 25, envir = dada_opts)
assign("HOMOPOLYMER_GAPPING", FALSE, envir = dada_opts) # NOT YET IMPLEMENTED
assign("SCORE_MATRIX", matrix(c(5, -4, -4, -4, -4, 5, -4, -4, -4, -4, 5, -4, -4, -4, -4, 5),
                              nrow=4, byrow=TRUE), envir = dada_opts)  # UNUSED. NECESSARY?
assign("GAP_PENALTY", -8, envir = dada_opts) # UNUSED. NECESSARY?

################################################################################
#' Set DADA options
#'
#' @param option (Required). Character.
#'  The DADA options to set.
#' 
#' @param value (Required). Corresponding variables of the appropriate data type.
#'  option and value must have the same length.
#'  WARNING: DATA TYPE CONSISTENCY IS NOT CURRENTLY CHECKED!!
#' 
#' @export
#'
set_dada_opt <- function(option, value) {
  if(!all(option %in% ls(dada_opts))) {
    stop(paste("Invalid DADA option:", option))
  }
  
  if(length(option) != length(value)) {
    stop(paste("Different lengths of name/value vectors:", length(option), ",", length(value)))
  }
  
  for(i in seq(length(option))) {
    assign(option[[i]], value[[i]], envir=dada_opts)
  }
}

################################################################################
#' Get DADA options
#'
#' @param option (Optional). Character.
#'  The DADA options to get.
#' 
#' @return Named list of option/value pairs.
#'  Returns NULL if an invalid option is requested.
#' 
#' @export
#'
get_dada_opt <- function(option = NULL) {
  if(is.null(option)) option <- ls(dada_opts)
  
  if(!all(option %in% ls(dada_opts))) {
    warning("Tried to get an invalid DADA option:", option[!(option %in% ls(dada_opts))])
    option <- option[option %in% ls(dada_opts)]
  }

  ropts <- lapply(option, function(x) get(x, envir=dada_opts))
  names(ropts) <- option
  if(length(ropts) == 1) ropts <- ropts[[1]]  # If just one option requested, return it alone
  return(ropts)
}


#define GAPPEN -8
#define BAND 50 // Size of band in banded alignments. 0 means no banding.
#define USE_KMERS 1
#define KMER_SIZE 6
#define KDIST_CUTOFF 0.5
#define OMEGA_A 0.01



