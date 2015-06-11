################################################################################
#' DEPRECATED IN FAVOR OF ISBIMERA2: FASTER AND FINDS SHIFT-BIMERAS AS WELL
#' 
#' Determine if input sequence is an exact bimera of putative parent sequences.
#' 
#' This function attempts to find an exact bimera of the parent sequences that
#' matches the input sequence. A bimera is a two-parent chimera, in which the
#' left side is made up of one parent sequence, and the right-side made up of
#' a second parent sequence. If an exact bimera is found (no mismatches are
#' allowed and positions must line up exactly) TRUE is returned, otherwise FALSE.
#' 
#' WARNING: This is a very simplistic bimera checker that does no alignment and
#' allows no mismatches. While it seems to perform well on DADA denoised Illumina
#' data, more caution is warranted for 454 data where alignment may be necessary.
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
#' @seealso \code{\link{isBimeraDenovo}}
#'
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

#' @export
#' 
getOverlaps <- function(parent, sq) {  # parent must be first, sq is being evaluated as a potential bimera
  al <- nwalign(parent, sq)
  lr <- C_get_overlaps(al[1], al[2])
  return(lr)
}

################################################################################
#' DEPRECATED IN FAVOR OF ISBIMERA2: FASTER AND FINDS SHIFT-BIMERAS AS WELL
#' 
#' Determine if input sequence is an exact bimera of putative parent sequences.
#' 
#' This function attempts to find an exact bimera of the parent sequences that
#' matches the input sequence. A bimera is a two-parent chimera, in which the
#' left side is made up of one parent sequence, and the right-side made up of
#' a second parent sequence. If an exact bimera is found (no mismatches are
#' allowed, but strict shifts are allowed) TRUE is returned, otherwise FALSE.
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
#' @seealso \code{\link{isBimeraDenovo}}
#' @export
#' 
isBimera2 <- function(sq, parents, verbose=FALSE) { # Note that order of sq/parents is reversed here from internals
  ov <- t(unname(sapply(parents, function(x) getOverlaps(x, sq))))
  return((max(ov[,1]) + max(ov[,2])) >= nchar(sq))
  # TRUE if the max left/right overlaps are as long as the sequence (or longer)
}

################################################################################
#' Identify bimeras de-novo from collections of unique sequences.
#' 
#' This function is a wrapper around isBimera for collections of DADA denoised
#' sequences. Each sequence is evaluated against a set of "parents" drawn from the
#' sequence collection that are sufficiently more abundant than the sequence
#' being evaluated.
#' 
#' @param unqs (Required). A named integer vector or data.frame with "sequence" and "abundance" cols.
#'   This contains the sequences (and corresponding abundances) to be checked for bimeras.
#'   
#' @param minFoldParentOverAbundance (Optional). A \code{numeric(1)}. Default is 10.
#'   Only sequences at least this much-fold more abundant than a sequence can be its "parents".
#'   
#' @param minParentAbundance (Optional). A \code{numeric(1)}. Default is 100.
#'   Only sequences at least this abundant can be "parents".
#' 
#' @seealso \code{\link{isBimera}}
#' 
#' @export
#' 
isBimeraDenovo <- function(unqs, minFoldParentOverAbundance = 10, minParentAbundance = 100, verbose=FALSE) {
  if(is.data.frame(unqs) && "sequence" %in% colnames(unqs) && "abundance" %in% colnames(unqs)) {
    seqs <- unqs$sequence
    abunds <- unqs$abundance
  } else if(is.integer(unqs) && !is.null(names(unqs))) {
    seqs <- names(unqs)
    abunds <- unname(unqs)
  } else {
    stop("Improper input: Requires named integer vector or $clustering data.frame.")
  }
  
  loopFun <- function(sq, abund) {
    pars <- seqs[(abunds>(minFoldParentOverAbundance*abund) & abunds>minParentAbundance)]
    if(length(pars) == 0) {
      if(verbose) print("No possible parents.")
      return(FALSE)
    } else if (length(pars) == 1) {
      if(verbose) print("Only one possible parent.")
      return(FALSE)
    } else {
      isBimera2(sq, pars, verbose)
    }
  }
  bims <- mapply(loopFun, seqs, abunds)
  return(bims)
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
