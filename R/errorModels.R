#' A loess fit for the error rates
#' 
#' @export
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
        mod.lo <- loess(rlogp ~ q, df, weight=errs)
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

#' Make function that accepts the output transitions matrix from DADA ($subqual)
#'  and returns inferred error rates for each transition type (t_ij) and quality
#'  
#' @export
makeErrfun <- function(err_model, init_params, inflation=1.0) {
  errfun <- function(trans, init = init_params) {
    qq <- as.numeric(colnames(trans))
    est <- matrix(0, nrow=0, ncol=length(qq))
    for(nti in c("A","C","G","T")) {
      for(ntj in c("A","C","G","T")) {
        if(nti != ntj) {
          x <- qq
          numer <- trans[paste0(nti,"2",ntj),]
          denom <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
          # Determine which entries to keep, no zeros or NAs or low-outlier denoms
          outs <- 10^(boxplot.stats(log10(denom[denom>0]))$out)
          outs <- outs[outs<median(denom)] # low outliers
          keep <- !(is.na(numer) | is.na(denom) | denom==0 | (denom %in% outs))
          x <- x[keep]
          y <- (numer[keep]+1)/(4+denom[keep]) # Pseudocounts... REVISIT
          mod.tp <- optim(init, makeLoglsqObj(err_model,x,y))
          pred <- err_model(mod.tp$par, qq)
          est <- rbind(est, pred)
        } # if(nti != ntj)
      } # for(ntj in c("A","C","G","T"))
    } # for(nti in c("A","C","G","T"))
    
    # Expand the err matrix with the self-transition probs
    err <- rbind(1-colSums(est[1:3,]), est[1:3,],
                 est[4,], 1-colSums(est[4:6,]), est[5:6,],
                 est[7:8,], 1-colSums(est[7:9,]), est[9,],
                 est[10:12,], 1-colSums(est[10:12,]))
    rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
    colnames(err) <- colnames(trans)
    
    # Inflate
    err <- inflateErr(err, inflation, FALSE)
    # Return
    return(err)
  } # errfun <- function(trans, init = init_params)
}

#' An error model for the dependence of error rates on quality in Illumina data.
#' 
#' A piecewise linear model with three segments. The first and third segments are flat.
#' 
#' @export
threepiece <- function(parms, qave) {
  if(!(is.numeric(parms) && length(parms) == 4)) {
    stop("Invalid parms argument to threepiece.")
  }
  e0 <- parms[1]
  knot1 <- parms[2]
  knot2 <- parms[3]
  e40 <- parms[4]
  pred <- rep(e0, length(qave))
  pred[qave>knot1] <- e0 + (qave[qave>knot1]-knot1)*(e40-e0)/(knot2-knot1)
  pred[qave>knot2] <- e40
  pred <- 10^pred
  pred[pred > 0.25] <- 0.25 ### Hacky
  return(pred)
}

# Suggested parameters to initialize threepiece with when using optim().
tp_init_parms <- c(-1, 20, 35, -4)

#' Creates an objective function suitable for optim().
#'   Sums the square errors after taking the log transform.
makeLoglsqObj <- function(err_model, qave, obs) {
  objfun <- function(parms) {
    pred <- err_model(parms, qave)
    return(sum((log(pred)-log(obs))^2))
  }
  return(objfun)
}

#' @export
#' @import ggplot2
#' @import gridExtra
showErrors <- function(dq, nti, ntj="all", erri=TRUE, erro=TRUE, ...) {
  ACGT <- c("A", "C", "G", "T")
  if(ntj == "all") {
    err_plots <- mapply(function(x,y) .showErrors(dq, x, y, erri, erro, title=paste(x, "->", y)), nti, ACGT, SIMPLIFY=FALSE)
    do.call(grid.arrange, c(err_plots, ncol=2))
  } else {
    .showErrors(dq, nti, ntj, erri, erro)
  }
}

#' @import ggplot2
.showErrors <- function(dq, nti, ntj, erri=TRUE, erro=TRUE, ...) {
  ACGT <- c("A", "C", "G", "T")
  tij <- 4*(which(ACGT==nti)-1) + which(ACGT==ntj)
  nij <- paste0(nti,"2",ntj)
  subdf <- as.data.frame(t(dq$trans))
  subdf$qave <- as.numeric(rownames(subdf))
  subdf[,nij] <- subdf[,nij]/(subdf[,paste0(nti,"2A")]+subdf[,paste0(nti,"2C")]+subdf[,paste0(nti,"2G")]+subdf[,paste0(nti,"2T")])

  if(erri) {
    erridf <- dq$err_in
    if(is.list(erridf)) erridf <- erridf[[1]]
    erridf <- as.data.frame(t(erridf))
    colnames(erridf) <- paste0(rep(ACGT, each=4), "2", ACGT)
    rownames(erridf) <- seq(dq$opts$QMIN, dq$opts$QMAX, length.out=nrow(erridf))
    erridf$qave <- as.numeric(rownames(erridf))
  }
  
  if(!("err_out" %in% ls(dq))) erro <- FALSE
  if(erro) {
    errodf <- as.data.frame(t(dq$err_out))
    errodf$qave <- as.numeric(rownames(errodf))
  }
  
  p <- ggplot(data=subdf, aes_string(x="qave", y=paste0("log10(",nij,")")), ...)
  p <- p + geom_point()
  if(erri) p <- p + geom_line(data=erridf)
  if(erro) p <- p + geom_line(data=errodf, linetype="dashed")
  return(p)
}

#' Function that inflates an error rate matrix by a specified factor. Accounts for saturation.
#'   new_err_rate <- err_rate * inflate / (1 + (inflate-1) * err_rate)
#' @export
inflateErr <- function(err, inflation=1.0, inflateNonSubs = FALSE) {
  t_errs <- c("A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G")
  err[t_errs,] <- (err[t_errs,] * inflation)/(1 + (inflation-1) * err[t_errs,])
  if(inflateNonSubs) { # Also inflate the non-substitution probabilities
    t_nonsubs <- c("A2A", "C2C", "G2G", "T2T")
    err[t_nonsubs,] <- (err[t_nonsubs,] * inflation)/(1 + (inflation-1) * err[t_nonsubs,])
  }
  return(err)
}

################################################################################
#' Identify False Positive inferred sequences due to bad bases.
#' 
#' Illumina sequencing sometimes produces "bad bases", positions at which 
#' error rates are significantly higher than expected by the assigned quality
#' score. This function identifies the inferred sequences that are likely to
#' have been driven by those bad bases.
#' 
#' @param dadaO (Required). Output of dada() function.
#' 
#' @param minFraction (Optional). A \code{numeric(1)}. Default is 0.51.
#'  The minimum fraction of bad bases among the base positions used to infer the
#'  sequence required to call the inferred sequence a false positive.
#'   
#' @param omegaB (Optional). A \code{numeric(1)}. Default is 1e-10.
#'  The p-value threshhold below which a base is assigned as "bad".
#'  The p-value is calculated by the number of repeated occurences of a particular
#'    base position individually driving the formation of a new cluster. Bad bases 
#'    drive many new "1-away" clustes.
#'  The null hypothesis being tested is that real differences are distributed
#'    uniformly along the sequence. This is not true, biological differences are
#'    non-uniform, so this pvalue threshhold should be set conservatively.
#' 
#' @param minOccurence (Optional). A \code{numeric(1)}. Default is 4.
#'  The minimum times a single base position must drive the formation of a new cluster
#'    before it can be considered a "bad base".
#'
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Defaults FALSE.
#'
#' @return Logical vector of length the number of inferred sequences. 
#'  TRUE if inferred sequence a false positive.
#'  FALSE otherwise.
#'
#' @seealso \code{\link{getBadBases}}
#'
#' @export
#' 
getBadBaseFPs <- function(dadaO, minFraction = 0.51, omegaB = 1e-10, minOccurence = 4, verbose=FALSE) {
  bb <- getBadBases(dadaO, omegaB, minOccurence)
  fps <- c("1"=FALSE, tapply(dadaO$subpos$pos, dadaO$subpos$clust, function(x) mean(x %in% bb) >= minFraction))
  if(verbose) {
    cat(length(fps), "false positives caused by bad bases identified.\n")
  }
  unname(fps)
}

################################################################################
#' Identify bad base positions.
#' 
#' Illumina sequencing sometimes produces "bad bases", positions at which 
#' error rates are significantly higher than expected by the assigned quality
#' score. This function identifies those bad bases.
#' 
#' @param dadaO (Required). Output of dada() function.
#' 
#' @param omegaB (Optional). A \code{numeric(1)}. Default is 1e-10.
#'  The p-value threshhold below which a base is assigned as "bad".
#'  The p-value is calculated by the number of repeated occurences of a particular
#'    base position individually driving the formation of a new cluster. Bad bases 
#'    drive many new "1-away" clustes.
#'  The null hypothesis being tested is that real differences are distributed
#'    uniformly along the sequence. This is not true, biological differences are
#'    non-uniform, so this pvalue threshhold should be set conservatively.
#' 
#' @param minOccurence (Optional). A \code{numeric(1)}. Default is 4.
#'  The minimum times a single base position must drive the formation of a new cluster
#'    before it can be considered a "bad base".
#'
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Defaults FALSE.
#'
#' @return Integer vector of the bad base positions. 
#'
#' @seealso \code{\link{getBadBases}}
#'
#' @export
#' 
getBadBases <- function(dadaO, omegaB = 1e-10, minOccurence = 4, verbose=FALSE) {
  oos <- which(dadaO$clustering$birth_ham == 1)
  oopos <- dadaO$subpos[dadaO$subpos$clust %in% oos,]
  tab <- table(oopos$pos)
  seqlen <- nchar(dadaO$clustering$sequence[[1]]) ##
  posp <- ppois(tab, length(oos)/seqlen, lower.tail=FALSE) * seqlen
  bad_bases <- as.integer(names(posp)[posp<omegaB & tab>=minOccurence])
  if(verbose) {
    cat(length(bad_bases), "bad bases identified.\n")
  }
  return(bad_bases)
}

