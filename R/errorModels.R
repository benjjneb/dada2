#' Use a loess fit to estimate error rates from transition counts.
#'
#' This function accepts a matrix of observed transitions, with each transition
#' corresponding to a row (eg. row 2 = A->C) and each column to a quality score
#' (eg. col 31 = Q30). It returns a matrix of estimated error
#' rates of the same shape. Error rates are estimates by a \code{\link{loess}} fit
#' of the observed rates of each transition as a function of the quality score.
#' Self-transitions (i.e. A->A) are taken to be the left-over probability.
#'
#' @param trans (Required). A matrix of the observed transition counts. Must be 16 rows,
#' with the rows named "A2A", "A2C", ...
#'
#' @return A numeric matrix with 16 rows and the same number of columns as trans.
#'  The estimated error rates for each transition (row, eg. "A2C") and quality score
#'  (column, eg. 31), as determined by \code{\link{loess}} smoothing over the quality
#'  scores within each transition category.
#'
#' @export
#'
#' @examples
#' derep1 <- derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1)
#' err.new <- loessErrfun(dada1$trans)
#'
loessErrfun <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights=errs)
        #        mod.lo <- loess(rlogp ~ q, df)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))

  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE

  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

#' Inflates an error rate matrix by a specified factor, while accounting for saturation.
#'
#' Error rates are "inflated" by the specified factor, while appropriately saturating so that rates
#' cannot exceed 1. The formula is:
#'   new_err_rate <- err_rate * inflate / (1 + (inflate-1) * err_rate)
#'
#' @param err (Required). A numeric matrix of transition rates (16 rows, named "A2A", "A2C", ...).
#'
#' @param inflation (Required). The fold-factor by which to inflate the transition rates.
#'
#' @param inflateSelfTransitions (Optional). Default FALSE.
#'  If True, self-transitions (eg. A->A) are also inflated.
#'
#' @return An error rate matrix of the same dimensions as the input error rate matrix.
#'
#' @export
#'
#' @examples
#'  tperr2 <- inflateErr(tperr1, 2)
#'  tperr3.all <- inflateErr(tperr1, 3, inflateSelfTransitions=TRUE)
#'
inflateErr <- function(err, inflation, inflateSelfTransitions = FALSE) {
  t_errs <- c("A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G")
  err[t_errs,] <- (err[t_errs,] * inflation)/(1 + (inflation-1) * err[t_errs,])
  if(inflateSelfTransitions) { # Also inflate the non-substitution probabilities
    t_nonsubs <- c("A2A", "C2C", "G2G", "T2T")
    err[t_nonsubs,] <- (err[t_nonsubs,] * inflation)/(1 + (inflation-1) * err[t_nonsubs,])
  }
  return(err)
}

################################################################################
#  --------------------- REQUIRES FURTHER TESTING --------------------------
# Identify False Positive inferred sequences due to bad bases.
#
# Illumina sequencing sometimes produces "bad bases", positions at which
# error rates are significantly higher than expected by the assigned quality
# score. This function identifies the inferred sequences that are likely to
# have been driven by those bad bases.
#
# @param clust (Required). The $clustering data frame from the dada() output.
#   May be subsetted from the original prior to using this function.
#
# @param birth_subs (Required). The $birth_subs data frame from the dada() output.
#
# @param minFraction (Optional). A \code{numeric(1)}. Default is 0.51.
#  The minimum fraction of bad bases among the base positions used to infer the
#  sequence required to call the inferred sequence a false positive.
#
# @param omegaB (Optional). A \code{numeric(1)}. Default is 1e-10.
#  The p-value threshold below which a base is assigned as "bad".
#  The p-value is calculated by the number of repeated occurrences of a particular
#    base position individually driving the formation of a new cluster. Bad bases
#    drive many new "1-away" clusters.
#  The null hypothesis being tested is that real differences are distributed
#    uniformly along the sequence. This is not true, biological differences are
#    non-uniform, so this pvalue threshold should be set conservatively.
#
# @param minOccurence (Optional). A \code{numeric(1)}. Default is 4.
#  The minimum times a single base position must drive the formation of a new cluster
#    before it can be considered a "bad base".
#
# @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#
# @return Logical vector of length the number of inferred sequences.
#  TRUE if inferred sequence a false positive.
#  FALSE otherwise.
#
# @seealso \code{\link{getBadBases}}
#
isBadBaseFP <- function(clust, birth_subs, minFraction = 0.51, omegaB = 1e-10, minOccurence = 4, verbose=FALSE) {
  bb <- getBadBases(clust, birth_subs, omegaB, minOccurence, verbose=verbose)
  fps <- tapply(birth_subs$pos, birth_subs$clust, function(x) mean(x %in% bb) >= minFraction)
  fps <- names(fps)[fps]
  rval <- rownames(clust) %in% fps
  if(verbose) {
    cat(sum(rval), "false positives caused by bad bases identified from", nrow(clust), "input sequences.\n")
  }
  rval
}

################################################################################
#  --------------------- REQUIRES FURTHER TESTING --------------------------
# Identify bad base positions.
#
# Illumina sequencing sometimes produces "bad bases", positions at which
# error rates are significantly higher than expected by the assigned quality
# score. This function identifies those bad bases.
#
# @param clust (Required). The $clustering data frame from the dada() output.
#   May be subsetted from the original prior to using this function.
#
# @param birth_subs (Required). The $birth_subs data frame from the dada() output.
#
# @param omegaB (Optional). A \code{numeric(1)}. Default is 1e-10.
#  The p-value threshold below which a base is assigned as "bad".
#  The p-value is calculated by the number of repeated occurrences of a particular
#    base position individually driving the formation of a new cluster. Bad bases
#    drive many new "1-away" clustes.
#  The null hypothesis being tested is that real differences are distributed
#    uniformly along the sequence. This is not true, biological differences are
#    non-uniform, so this pvalue threshhold should be set conservatively.
#
# @param minOccurence (Optional). A \code{numeric(1)}. Default is 4.
#  The minimum times a single base position must drive the formation of a new cluster
#    before it can be considered a "bad base".
#
# @param verbose (Optional). \code{logical(1)} indicating verbose text output. Defaults FALSE.
#
# @return Integer vector of the bad base positions.
#
# @seealso \code{\link{isBadBaseFP}}
#
getBadBases <- function(clust, birth_subs, omegaB = 1e-20, minOccurence = 4, verbose=FALSE) {
  oos <- which(clust$birth_ham == 1)
  oopos <- birth_subs[birth_subs$clust %in% oos,]
  tab <- table(oopos$pos)
  if(length(unique(nchar(clust$sequence)))>1) stop("Requires same length sequences.")
  seqlen <- nchar(clust$sequence[[1]])
  posp <- ppois(tab, length(oos)/seqlen, lower.tail=FALSE) * seqlen
  bad_bases <- as.integer(names(posp)[posp<omegaB & tab>=minOccurence])
  if(verbose) {
    cat(length(bad_bases), "bad bases identified.\n")
  }
  return(bad_bases)
}

#' An empirical error matrix.
#'
#' A dataset containing the error matrix estimated by fitting a piecewise linear model to
#' the errors observed in the mock community featured in Schirmer 2015 (metaID 35).
#'
#' @format A numerical matrix with 16 rows and 41 columns.
#'  Rows correspond to the 16 transition (eg. A2A, A2C, ...)
#'  Columns correspond to consensus quality scores 0 to 40.
#'
#' @name tperr1
NULL

#' An empirical error matrix.
#'
#' A dataset containing the error matrix estimated by DADA2 from the forward reads of the
#' Illumina Miseq 2x250 sequenced Balanced mock community (see manuscript).
#'
#' @format A numerical matrix with 16 rows and 41 columns.
#'  Rows correspond to the 16 transition (eg. A2A, A2C, ...)
#'  Columns correspond to consensus quality scores 0 to 40.
#'
#' @name errBalancedF
NULL

#' An empirical error matrix.
#'
#' A dataset containing the error matrix estimated by DADA2 from the reverse reads of the
#' Illumina Miseq 2x250 sequenced Balanced mock community (see manuscript).
#'
#' @format A numerical matrix with 16 rows and 41 columns.
#'  Rows correspond to the 16 transition (eg. A2A, A2C, ...)
#'  Columns correspond to consensus quality scores 0 to 40.
#'
#' @name errBalancedR
NULL

#' An empirical error matrix.
#'
#' A dataset containing the error matrix estimated by DADA2 from the forward reads of the
#' Illumina Miseq 2x250 sequenced HMP mock community (see manuscript).
#'
#' @format A numerical matrix with 16 rows and 41 columns.
#'  Rows correspond to the 16 transition (eg. A2A, A2C, ...)
#'  Columns correspond to consensus quality scores 0 to 40.
#'
#' @name errHmpF
NULL

#' An empirical error matrix.
#'
#' A dataset containing the error matrix estimated by DADA2 from the reverse reads of the
#' Illumina Miseq 2x250 sequenced HMP mock community (see manuscript).
#'
#' @format A numerical matrix with 16 rows and 41 columns.
#'  Rows correspond to the 16 transition (eg. A2A, A2C, ...)
#'  Columns correspond to consensus quality scores 0 to 40.
#'
#' @name errHmpR
NULL

#' An empirical error matrix.
#'
#' A dataset containing the error matrix estimated by DADA2 from the forward reads of the
#' Illumina Miseq 2x250 sequenced Extreme mock community (see manuscript).
#'
#' @format A numerical matrix with 16 rows and 41 columns.
#'  Rows correspond to the 16 transition (eg. A2A, A2C, ...)
#'  Columns correspond to consensus quality scores 0 to 40.
#'
#' @name errExtremeF
NULL

#' An empirical error matrix.
#'
#' A dataset containing the error matrix estimated by DADA2 from the reverse reads of the
#' Illumina Miseq 2x250 sequenced Extreme mock community (see manuscript).
#'
#' @format A numerical matrix with 16 rows and 41 columns.
#'  Rows correspond to the 16 transition (eg. A2A, A2C, ...)
#'  Columns correspond to consensus quality scores 0 to 40.
#'
#' @name errExtremeR
NULL
