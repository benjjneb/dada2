################################################################################
#' Determine if input sequence is an exact bimera of putative parent sequences.
#' 
#' This function attempts to find an exact bimera of the parent sequences that
#' matches the input sequence. A bimera is a two-parent chimera, in which the
#' left side is made up of one parent sequence, and the right-side made up of
#' a second parent sequence. If an exact bimera is found (no mismatches are
#' allowed and positions must line up exactly) TRUE is returned, otherwise FALSE.
#' 
#' @param sq (Required). A \code{character(1)}.
#'  The sequence being evaluated as a possible bimera.
#' 
#' @param parents (Required). Character vector.
#'  A vector of possible "parent" sequence that could form the left and right
#'  sides of the bimera.
#'   
#' @param verbose (Optional). \code{logical(1)}.
#'  If TRUE, some informative output is printed. Default is FALSE.
#'
#' @return \code{logical(1)}.
#'  TRUE if sq is an exact bimera of two of the parents. Otherwise FALSE.
#'
#' @export
#' @import Biostrings 
#' 
isBimera <- function(sq, parents, verbose=FALSE) {
  same <- t(hasLetterAt(DNAStringSet(parents), sq, seq(nchar(sq))))
  # Output is logical matrix, rows:position, columns:parent, value:match
  lhs <- c(max(same[1,]), sapply(seq(2, nrow(same)), function(x) max(apply(same[1:x,], 2, sum))))
  min_ham <- nchar(sq) - max(lhs)
  rhs <- c(0,  max(same[nrow(same),]), sapply( seq(2, nrow(same)-1), function(x) max(apply(same[(nrow(same)-x+1):nrow(same),], 2, sum)) ) )
  # rhs has a (] bound, while lhs has a [] bound on indices
  most_matches = max(lhs+rev(rhs))
  min_bi_ham <- nchar(sq) - most_matches
  if(verbose) {
    cat("Sequence has min hamming distance of", min_ham, "to a parent, and", min_bi_ham, "to a bimera:", (min_bi_ham == 0 && min_ham > 0), "\n")
  }
  return(min_bi_ham == 0 && min_ham > 0)
}

# DEPRECATED: USES SEQINR, BUT A BIT FASTER THAN NEW VERSION
#isBimera <- function(sq, parents, verbose=FALSE) {
#  require(seqinr)
#  sq <- s2c(sq)
#  parents <- sapply(parents, s2c) # returns a character matrix with the parents as columns
#  same <- (parents == sq)
#  lhs <- c(max(same[1,]), sapply(seq(2, nrow(same)), function(x) max(apply(same[1:x,], 2, sum))))
#  min_ham <- length(sq) - max(lhs)
#  rhs <- c(0,  max(same[nrow(same),]), sapply( seq(2, nrow(same)-1), function(x) max(apply(same[(nrow(same)-x+1):nrow(same),], 2, sum)) ) )
#  # rhs has a (] bound, while lhs has a [] bound on indices
#  most_matches = max(lhs+rev(rhs))
#  min_bi_ham <- length(sq) - most_matches
#  if(verbose) {
#    cat("Sequence has min hamming distance of", min_ham, "to a parent, and", min_bi_ham, "to a bimera:", (min_bi_ham == 0 && min_ham > 0), "\n")
#  }
#  return(min_bi_ham == 0 && min_ham > 0)
#}
