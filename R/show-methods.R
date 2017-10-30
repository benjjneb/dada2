############################################################################
#' method extensions to show for dada2 objects.
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior. 
#'
#' @seealso \code{\link[methods]{show}}
#' 
#' @inheritParams methods::show
#' @return NULL.
#' @importFrom stats median
#' @rdname show-methods
#' @include allClasses.R
#' # examples
setMethod("show", "derep", function(object){
  cat("derep-class: R object describing dereplicated sequencing reads", fill=TRUE)
  if( length(object$uniques) > 0 ){
    cat("$uniques:", sum(object$uniques, na.rm = TRUE), "reads in",
        length(names(object$uniques)), "unique sequences\n")
    seqlens <- nchar(names(object$uniques))
    cat("  Sequence lengths: min=", min(seqlens), ", median=", median(seqlens),
        ", max=", max(seqlens), sep="", fill=TRUE)
  }
  if( length(object$quals) > 0 & inherits(object$quals, "matrix")){
    cat("$quals: Quality matrix dimension: ", dim(object$quals), fill = TRUE)
    quals <- as.vector(object$quals)
    cat("  Consensus quality scores: min=", min(quals, na.rm=TRUE), ", median=", median(quals, na.rm=TRUE),
        ", max=", max(quals, na.rm=TRUE), sep="", fill=TRUE)
  }
  if( length(object$map) > 0 ){
    cat("$map: Map from reads to unique sequences: ", object$map[1L:5L], "...", fill = TRUE)
  }
})

#' @inheritParams methods::show
#' @rdname show-methods
#' @include allClasses.R
setMethod("show", "dada", function(object){
  cat("dada-class: object describing DADA2 denoising results", fill=TRUE)
  if( length(object$denoised) > 0 && length(object$map) > 0 ) {
    cat(length(object$denoised), "sequence variants were inferred from", length(object$map), "input unique sequences.", fill=TRUE)
  }
  cat("Key parameters: OMEGA_A = ", object$opts$OMEGA_A, ", BAND_SIZE = ", 
      object$opts$BAND_SIZE, ", USE_QUALS = ", object$opts$USE_QUALS, 
      sep="", fill=TRUE)
})

############################################################################
#' Deactivate renaming of derep-class objects.
#'
#' @inheritParams base::`names<-`
#' @return NULL.
#' @include allClasses.R
#' # examples
setMethod("names<-", "derep", function(x, value){
  warning("derep-class objects cannot be renamed.")
  return(x)
})

############################################################################
#' Deactivate renaming of dada-class objects.
#'
#' @inheritParams base::`names<-`
#' @return NULL.
#' @include allClasses.R
#' # examples
setMethod("names<-", "dada", function(x, value){
  warning("dada-class objects cannot be renamed.")
  return(x)
})

############################################################################
#' Change concatenation to list construction.
#'
#' @inheritParams base::c
#' @param x A derep-class object
#' @return list.
#' @include allClasses.R
#' # examples
setMethod("c", signature("derep"), function(x,...,recursive=FALSE){
  list(x,...)
})

############################################################################
#' Change concatenation to list construction.
#'
#' @inheritParams base::c
#' @param x A dada-class object
#' @return list.
#' @include allClasses.R
#' # examples
setMethod("c", signature("dada"), function(x,...,recursive=FALSE){
  list(x,...)
})
