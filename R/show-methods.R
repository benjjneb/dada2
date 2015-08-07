############################################################################
#' method extensions to show for dada2 objects.
#'
#' See the general documentation of \code{\link[methods]{show}} method for
#' expected behavior. 
#'
#' @seealso \code{\link[methods]{show}}
#' 
#' @inheritParams methods::show
#' @export
#' @rdname show-methods
#' @include allClasses.R
#' @examples
#' # examples
setMethod("show", "derep", function(object){
  cat("derep-class: object describing dereplicated sequences", fill=TRUE)
  if( length(object$uniques) > 0 ){
    cat("Sum of counts:", sum(object$uniques, na.rm = TRUE), fill = TRUE)
  }
  cat("Uniques (`$uniques`) length: ", length(object$uniques), fill = TRUE)
  if( length(names(object$uniques)) > 0 ){
    cat(names(object$uniques)[1], "\n...", fill = TRUE)
  }
  if( length(object$quals) > 0 & inherits(object$quals, "matrix")){
    cat("Quality Matrix, (`$quals`) dimension: ", dim(object$quals), fill = TRUE)
  }
  if( length(object$map) > 0 ){
    cat("Dereplication Map (`$map`): ", object$map[1L:5L], "...", fill = TRUE)
  }
})
############################################################################
#' @rdname show-methods
setMethod("show", "dada", function(object){
  cat("dada-class: object describing DADA2 denoising results", fill=TRUE)
})
