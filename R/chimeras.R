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
#' @seealso 
#'  \code{\link{isBimeraDenovo}}, \code{\link{removeBimeraDenovo}}
#'  
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' sqs1 <- getSequences(derep1)
#' isBimera(sqs1[[20]], sqs1[1:10])
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
#' Identify bimeras from collections of unique sequences.
#' 
#' This function is a wrapper around \code{\link{isBimera}} for collections of unique
#' sequences (i.e. sequences with associated abundances). Each sequence is evaluated 
#' against a set of "parents" drawn from the sequence collection that are sufficiently
#' more abundant than the sequence being evaluated. A logical vector is returned, with
#' an entry for each input sequence indicating whether it was (was not) consistent with
#' being a bimera of those more abundant "parents".
#' 
#' @param unqs (Required). A \code{\link{uniques-vector}} or any object that can be coerced
#'  into one with \code{\link{getUniques}}.
#'   
#' @param minFoldParentOverAbundance (Optional). A \code{numeric(1)}. Default is 1.
#'   Only sequences greater than this-fold more abundant than a sequence can be its 
#'   "parents".
#'   
#' @param minParentAbundance (Optional). A \code{numeric(1)}. Default is 8.
#'   Only sequences at least this abundant can be "parents".
#' 
#' @param allowOneOff (Optional). A \code{logical(1)}. Default is TRUE.
#'   If TRUE, sequences that have one mismatch or indel to an exact bimera are also
#'   flagged as bimeric.
#' 
#' @param minOneOffParentDistance (Optional). A \code{numeric(1)}. Default is 4.
#'   Only sequences with at least this many mismatches to the potential bimeric sequence
#'   considered as possible "parents" when flagging one-off bimeras. There is
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
#' @seealso 
#'  \code{\link{isBimera}}, \code{\link{removeBimeraDenovo}}
#' 
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' isBimeraDenovo(dada1)
#' isBimeraDenovo(dada1$denoised, minFoldParentOverAbundance = 2, allowOneOff=FALSE)
#' 
isBimeraDenovo <- function(unqs, minFoldParentOverAbundance = 1, minParentAbundance = 8, allowOneOff=TRUE, minOneOffParentDistance=4, maxShift = 16, verbose=FALSE, test=FALSE) {
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
      if(test) {
        C_is_bimera(sq, pars, maxShift)
      } else {
        isBimera(sq, pars, allowOneOff=allowOneOff, minOneOffParentDistance=minOneOffParentDistance, maxShift=maxShift)
      }
    }
  }
  bims <- mapply(loopFun, seqs, abunds)
  if(verbose) message("Identified ", sum(bims), " bimeras out of ", length(bims), " input sequences.")
  return(bims)
}

################################################################################
#' Remove bimeras from collections of unique sequences.
#' 
#' This function is a wrapper around \code{\link{isBimeraDenovo}}. Bimeras identified by
#' \link{isBimeraDenovo} are removed, and a bimera-free collection of unique sequences is returned.
#' 
#' @param unqs (Required). A \code{\link{uniques-vector}} or any object that can be coerced
#'  into one with \code{\link{getUniques}}. A list of such objects can also be provided.
#'   
#' @param ... (Optional). Arguments to be passed to \code{\link{isBimeraDenovo}}.
#'   
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#'
#' @return A uniques vector, or an object of matching class if a data.frame or sequence table is provided.
#'  A list of such objects is returned if a list of input unqs was provided.
#'
#' @seealso \code{\link{isBimeraDenovo}}
#' 
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' out.nobim <- removeBimeraDenovo(dada1)
#' out.nobim <- removeBimeraDenovo(dada1$clustering, minFoldParentOverAbundance = 2, allowOneOff=FALSE)
#' 
removeBimeraDenovo <- function(unqs, ..., verbose=FALSE) {
  if(class(unqs)!="list") {
    unqs <- list(unqs)
  }
  outs <- list()
  for(i in seq_along(unqs)) {
    bim <- isBimeraDenovo(unqs[[i]], ..., verbose=verbose)
    # The following code is adapted from getUniques
    object <- unqs[[i]]
    if(is.integer(object) && length(names(object)) != 0 && !any(is.na(names(object)))) { # Named integer vector already
      outs[[i]] <- object[!bim]
    } else if(class(object) == "dada") {  # dada return 
      outs[[i]] <- object$denoised[!bim]
    } else if(class(object) == "derep") {
      outs[[i]] <- object$uniques[!bim]
    } else if(is.data.frame(object) && all(c("sequence", "abundance") %in% colnames(object))) {
      outs[[i]] <- object[!bim,]
    } else if(class(object) == "matrix" && !any(is.na(colnames(object)))) { # Tabled sequences
      outs[[i]] <- object[,!bim,drop=FALSE]
    } else {
      stop("Unrecognized format: Requires named integer vector, dada-class, derep-class, sequence matrix, or a data.frame with $sequence and $abundance columns.")
    }
  }
  names(outs) <- names(unqs)
  if(length(outs) == 1) {
    outs <- outs[[1]]
  }
  outs
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
#' This function is a wrapper around isShift for collections of unique
#' sequences. Each unique sequence is evaluated against a set of "parents" drawn from
#' the sequence collection that are more abundant than the sequence being evaluated.
#' 
#' @param unqs (Required). A \code{\link{uniques-vector}} or any object that can be coerced
#'  into one with \code{\link{getUniques}}.
#'   
#' @param minOverlap (Optional). A \code{numeric(1)}. Default is 20.
#'   Minimum overlap required to call something a shift.
#'   
#' @param flagSubseqs (Optional). A \code{logical(1)}. Default is FALSE.
#'   Whether or not to flag strict subsequences as shifts.
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
#' isShiftDenovo(dada1$denoised, minOverlap=50, verbose=TRUE)
#' 
isShiftDenovo <- function(unqs, minOverlap = 20, flagSubseqs=FALSE, verbose=FALSE) {
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
isShiftedPair <- function(sq1, sq2, minOverlap=20, flagSubseqs=FALSE) {
  al <- nwalign(sq1, sq2, band=-1)
  foo <- C_eval_pair(al[1], al[2])
  return((foo["match"] < nchar(sq1) || flagSubseqs) && (foo["match"] < nchar(sq2) || flagSubseqs) &&
           foo["match"] >= minOverlap && foo["mismatch"]==0 && foo["indel"]==0)
}

isShift <- function(sq, pars, minOverlap=20, flagSubseqs=FALSE) {
  return(any(sapply(pars, function(par) isShiftedPair(sq, par, minOverlap=minOverlap, flagSubseqs=flagSubseqs))))
}
