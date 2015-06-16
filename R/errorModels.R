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
showErrors <- function(dq, nti, ntj="all", erri=TRUE, erro=TRUE, ...) {
  require(ggplot2)
  require(gridExtra)
  ACGT <- c("A", "C", "G", "T")
  if(ntj == "all") {
    err_plots <- mapply(function(x,y) .showErrors(dq, x, y, erri, erro, title=paste(x, "->", y)), nti, ACGT, SIMPLIFY=FALSE)
    do.call(grid.arrange, c(err_plots, ncol=2))
  } else {
    .showErrors(dq, nti, ntj, erri, erro)
  }
}

.showErrors <- function(dq, nti, ntj, erri=TRUE, erro=TRUE, ...) {
  require(ggplot2)
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
