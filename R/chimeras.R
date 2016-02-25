################################################################################
#' Determine if input sequence is a bimera of putative parent sequences.
#' 
#' This function attempts to find an exact bimera of the parent sequences that
#' matches the input sequence. A bimera is a two-parent chimera, in which the
#' left side is made up of one parent sequence, and the right-side made up of
#' a second parent sequence. If an exact bimera is found TRUE is returned, 
#' otherwise FALSE. Bimeras that are one-off from exact are also identified if
#' the allowOneOff argument is TRUE.
#' 
#' @param sq (Required). A \code{character(1)}.
#'  The sequence being evaluated as a possible bimera.
#' 
#' @param parents (Required). Character vector.
#'  A vector of possible "parent" sequence that could form the left and right
#'  sides of the bimera.
#'   
#' @param allowOneOff (Optional). A \code{logical(1)}. Default is TRUE.
#'   If TRUE, sq will be identified as a bimera if it is one mismatch or indel away 
#'   from an exact bimera.
#' 
#' @param minOneOffParentDistance (Optional). A \code{numeric(1)}. Default is 4.
#'   Only sequences with at least this many mismatches to sq are considered as possible
#'   "parents" when flagging one-off bimeras. There is no such screen when identifying
#'   exact bimeras.
#'   
#' @param maxShift (Optional). A \code{numeric(1)}. Default is 16.
#'   Maximum shift allowed when aligning sequences to potential "parents".
#' 
#' @return \code{logical(1)}.
#'  TRUE if sq is a bimera of two of the parents. Otherwise FALSE.
#'
#' @seealso \code{\link{isBimeraDenovo}}
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' unqs1 <- getUniques(derep1)
#' isBimera(names(unqs1)[[20]], names(unqs1)[1:10])
#' 
isBimera <- function(sq, parents, allowOneOff=TRUE, minOneOffParentDistance=4, maxShift=16) { # Note that order of sq/parents is reversed here from internals
  # (l0,r0) or (l0,r0,l1,r1,match,mismatch,indel) if allowOneOff
  ov <- t(unname(sapply(parents, function(x) getOverlaps(x, sq, allowOneOff=allowOneOff, maxShift=maxShift))))
  # Remove identical (or strictly shifted) parents
  id_or_fullshift <- apply(ov[,1:2], 1, function(x) max(x) >= nchar(sq))
  ov <- ov[!id_or_fullshift,,drop=FALSE]
  # Return FALSE if now too few parents
  if(nrow(ov) <= 1) return(FALSE)
  # Return TRUE if the max left/right 0-overlaps are as long as the sequence
  if((max(ov[,1]) + max(ov[,2])) >= nchar(sq)) {
    return(TRUE)
  } else if(allowOneOff) {
    too_close <- (apply(ov[,6:7], 1, sum) < minOneOffParentDistance)
    ov <- ov[!too_close,,drop=FALSE]
    # Return FALSE if now too few parents
    if(nrow(ov) <= 1) return(FALSE)
    # Return TRUE if the max left/right 0/1-overlaps are as long as the sequence
    if((max(ov[,1]) + max(ov[,4])) >= nchar(sq) || (max(ov[,2]) + max(ov[,3])) >= nchar(sq)) {
      return(TRUE)
    }
  }
  # If haven't returned TRUE, then not a bimera
  return(FALSE)
}

################################################################################
#' Identify bimeras de-novo from collections of unique sequences.
#' 
#' This function is a wrapper around isBimera for collections of DADA denoised
#' sequences. Each sequence is evaluated against a set of "parents" drawn from the
#' sequence collection that are sufficiently more abundant than the sequence
#' being evaluated.
#' 
#' @param unqs (Required). A "uniques vector" or any object that can be coerced into one with \code{\link{getUniques}}.
#'   This named integer vector is named by the sequences to be checked for bimeras and valued by their abundances.
#'   
#' @param minFoldParentOverAbundance (Optional). A \code{numeric(1)}. Default is 1.
#'   Only sequences more than this-fold abundant than a sequence can be its "parents".
#'   Default is intentionally permissive, as aggressively removing chimeras
#'    is the conservative choice for downstream analysis.
#'   
#' @param minParentAbundance (Optional). A \code{numeric(1)}. Default is 8.
#'   Only sequences at least this abundant can be "parents".
#' 
#' @param allowOneOff (Optional). A \code{logical(1)}. Default is TRUE.
#'   If TRUE, sequences that have one mismatch or indel to an exact bimera are also
#'   flagged.
#' 
#' @param minOneOffParentDistance (Optional). A \code{numeric(1)}. Default is 4.
#'   Only sequences with at least this many mismatches to the potential bimeric sequence
#'   considered as possible "parents" when flagging one-off bimeras. Note that there is
#'   no such screen when considering exact bimeras.
#'   
#' @param maxShift (Optional). A \code{numeric(1)}. Default is 16.
#'   Maximum shift allowed when aligning sequences to potential "parents".
#' 
#' @return \code{logical} of length the number of input unique sequences.
#'  TRUE if sequence is a bimera of more abundant "parent" sequences. Otherwise FALSE.
#'
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#'
#' @seealso \code{\link{isBimera}}
#' 
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' isBimeraDenovo(dada1)
#' isBimeraDenovo(getUniques(dada1), minFoldParentOverAbundance = 2, allowOneOff=FALSE)
#' 
isBimeraDenovo <- function(unqs, minFoldParentOverAbundance = 1, minParentAbundance = 8, allowOneOff=TRUE, minOneOffParentDistance=4, maxShift = 16, verbose=FALSE) {
  unqs <- getUniques(unqs)
  abunds <- unname(unqs)
  seqs <- names(unqs)
  
  loopFun <- function(sq, abund) {
    pars <- seqs[(abunds>(minFoldParentOverAbundance*abund) & abunds>minParentAbundance)]
    if(length(pars) == 0) {
      return(FALSE)
    } else if (length(pars) == 1) {
      return(FALSE)
    } else {
      isBimera(sq, pars, allowOneOff=allowOneOff, minOneOffParentDistance=minOneOffParentDistance, maxShift=maxShift)
    }
  }
  bims <- mapply(loopFun, seqs, abunds)
  if(verbose) message("Identified ", sum(bims), " bimeras out of ", length(bims), " input sequences.")
  return(bims)
}

# Internal function that finds the best overlap between two sequences.
# Uses NW alignment with ends-free gapping.
getOverlaps <- function(parent, sq, allowOneOff=FALSE, maxShift=16) {  # parent must be first, sq is being evaluated as a potential bimera
  al <- nwalign(parent, sq, band=maxShift)
  lr <- C_get_overlaps(al[1], al[2], 0, maxShift)
  if(allowOneOff) {
    lr1 <- C_get_overlaps(al[1], al[2], 1, maxShift)
    diff1 <- C_eval_pair(al[1], al[2])
    lr <- c(lr,lr1,diff1)
  }
  return(lr)
}

################################################################################
#' Identify sequences that are identical to a more abundant sequence up to an
#' overall shift.
#' 
#' This function is a wrapper around isShift for collections of DADA denoised
#' sequences. Each sequence is evaluated against a set of "parents" drawn from the
#' sequence collection that are more abundant than the sequence being evaluated.
#' 
#' @param unqs (Required). A "uniques vector" or any object that can be coerced into one with \code{\link{getUniques}}.
#'   This named integer vector is named by the sequences to be checked for bimeras and valued by their abundances.
#'   
#' @param minOverlap (Optional). A \code{numeric(1)}. Default is 20.
#'   Minimum overlap required to call something a shift.
#'   
#' @return \code{logical} of length the number of input unique sequences.
#'  TRUE if sequence is an exact shift of a more abundant sequence. Otherwise FALSE.
#'
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#'
#' @seealso \code{\link{isBimera}}
#' 
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' isShiftDenovo(dada1)
#' isShiftDenovo(getUniques(dada1), minOverlap=50, verbose=TRUE)
#' 
isShiftDenovo <- function(unqs, minOverlap = 20, verbose=FALSE) {
  unqs <- getUniques(unqs)
  abunds <- unname(unqs)
  seqs <- names(unqs)
  
  loopFun <- function(sq, abund) {
    pars <- seqs[abunds>abund]
    if(length(pars) == 0) {
      if(verbose) print("No possible parents.")
      return(FALSE)
    } else {
      isShift(sq, pars, minOverlap=minOverlap)
    }
  }
  shifts <- mapply(loopFun, seqs, abunds)
  return(shifts)
}


# Internal function that determines if two sequences are identical up to a shift
# Uses NW alignment with ends-free gapping
# 
# @param sq1 A \code{character(1)}. The first DNA sequence.
# 
# @param sq2 A \code{character(1)}. The second DNA sequence.
# 
# @param minOverlap (Optional). A \code{numeric(1)}. Default is 20.
#   Minimum overlap required to call something a shift.
#   
isShiftedPair <- function(sq1, sq2, minOverlap=20) {
  al <- nwalign(sq1, sq2, band=-1)
  foo <- C_eval_pair(al[1], al[2])
  return(foo["match"] < nchar(sq1) && foo["match"] < nchar(sq2) && foo["match"] >= minOverlap && foo["mismatch"]==0 && foo["indel"]==0)
}

isShift <- function(sq, pars, minOverlap=20) {
  return(any(sapply(pars, function(par) isShiftedPair(sq, par))))
}
