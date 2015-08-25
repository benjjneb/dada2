#' Plot Substitution Pairs from DADA Result
#' 
#' This is similar to original DADA article, Figure 6.
#' 
#' @param dadaOut (Required). The object returned by \code{\link{dada}}.
#' 
#' @param facetByGrp (Optional). Logical(1).
#'  Whether to plot all substitution groups together in one panel
#'  or separately on a grid of panels with a linear model fit.
#' 
#' @return A \code{\link{ggplot}2} object that will be rendered
#'  to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  See \code{\link{ggsave}} for additional options.
#' 
#' @importFrom reshape2 melt
#' @importFrom reshape2 dcast
#' @importFrom data.table data.table
#' @importFrom data.table setnames
#' @importFrom data.table dcast.data.table
#' @importFrom data.table :=
#' @import ggplot2
#' 
#' @export
#' 
#' @examples 
#' # Examples here.
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"), verbose = TRUE)
#' res1 <- dada(derep1,
#'              err = inflateErr(tperr1, 2), 
#'              errorEstimationFunction = loessErrfun, 
#'              selfConsist = TRUE) 
#' plotSubstitutions(res1)
#' 
plotSubstitutions = function(dadaOut, facetByGrp = TRUE){
  transdt = data.table(melt(dadaOut$trans))
  setnames(transdt, c("Substitution", "Quality", "Count"))
  transdt[, Sub1 := substr(Substitution, 1, 1)]
  transdt[, Sub2 := substr(Substitution, 3, 3)]
  # For a plot like from the DADA-1 paper, Figure 6,
  # Map color to complementary pairs of errors
  # red = (A→G,T→C) cyan=(C→T,G→A) green=(A→T,T→A) black=(C→A,G→T) blue=(A→C,T→G) purple=(C→G,G→C).
  CompSubGroups = c(A2G = "A2GT2C", T2C = "A2GT2C",
                    C2T = "C2TG2A", G2A = "C2TG2A",
                    A2T = "A2TT2A", T2A = "A2TT2A",
                    C2A = "C2AG2T", G2T = "C2AG2T",
                    A2C = "A2CT2G", T2G = "A2CT2G",
                    C2G = "C2GG2C", G2C = "C2GG2C") 
  transdt[, SubGrp := CompSubGroups[as.character(Substitution)]]
  # Redefine the Forward/Reverse pairings
  transdt[(substr(SubGrp, 1, 3) == Substitution), Direction := "Forward"]
  transdt[(substr(SubGrp, 4, 6) == Substitution), Direction := "Reverse"]
  # Define plot as in DADA1 article Figure 6
  # (but better because it is not matlab)
  tCast = dcast.data.table(data = transdt[(Direction != "NoChange")][Count > 0],
                           formula = SubGrp + Quality ~ Direction,
                           value.var = "Count",
                           drop = TRUE)
  # Remove missing
  tCast <- tCast[!is.na(Forward)][!is.na(Reverse)]
  # Define ggplot2 plot
  p1 = ggplot(tCast, aes(Forward, Reverse, color = SubGrp, size = Quality)) + 
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept=0, slope=1, alpha = 0.5, linetype = 2) +
    geom_point(alpha = 0.7) + 
    ggtitle("Substitution Pairs Count Comparison")
  if(facetByGrp){
    # If TRUE, facet and add simple linear model fit
    p1 <- p1 + 
      stat_smooth(method = "lm") +
      facet_wrap(~SubGrp, drop = TRUE)
  }
  return(p1)
}

#' Plot error rates from after dada processing.
#' 
#' This function plots the observed freuency of substitutions for each transition
#' (eg. A->C) as a function of the associated quality score. It also plots the error
#' rates passed in to this instantiation of the dada algorithm, and the output error
#' rates if they exist.
#' 
#' @param dq Required. A dada return object.
#' 
#' @param nti Optional. Default "all". Choice of "A", "C", "G", "T".
#' @param ntj Optional. Default "all". Choice of "A", "C", "G", "T".
#'  The error rate from nti->ntj will be plotted. If "all" is selected,
#'  the error rate from all 4 nti, or to all 4 ntj, will be plotted.
#'  
#' @param err_out Optional. Default TRUE.
#'  A \code{logical(1)} determining whether to plot the output error rates (solid).
#' 
#' @param err_in Optional. Default FALSE.
#'  A \code{logical(1)} determining whether to plot the input error rates (dashed).
#'
#' @param nominalQ Optional. Default FALSE.
#'  A \code{logical(1)} determining whether to plot the expected error rate (red) if the
#'  quality score exactly matched its nominal definition: Q = -10 log10(p_err).
#'  
#' @param ...
#'  Further arguments that are passed to the ggplot() function.
#'
#' @return A \code{\link{ggplot2}} object that will be rendered
#'  to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  
#' @export
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' 
plotErrors <- function(dq, nti="all", ntj="all", err_out=TRUE, err_in=FALSE, nominalQ=FALSE, ...) {
  ACGT <- c("A", "C", "G", "T")
  if(!(nti %in% c(ACGT, "all") && (ntj %in% c(ACGT, "all")))) {
    stop("nti and ntj must be a nucleotide (A/C/G/T) or all.")
  }
  if(nti == "all" && ntj == "all") {
    err_plots <- mapply(function(x,y) .plotErrors(dq, x, y, err_out, err_in, nominalQ, ...), rep(ACGT, each=4), rep(ACGT,4), SIMPLIFY=FALSE)
    do.call(grid.arrange, c(err_plots, ncol=4))
  }
  else if(nti == "all") {
    err_plots <- mapply(function(x,y) .plotErrors(dq, x, y, err_out, err_in, nominalQ, ...), ACGT, ntj, SIMPLIFY=FALSE)
    do.call(grid.arrange, c(err_plots, ncol=2))  
  }
  else if(ntj == "all") {
    err_plots <- mapply(function(x,y) .plotErrors(dq, x, y, err_out, err_in, nominalQ, ...), nti, ACGT, SIMPLIFY=FALSE)
    do.call(grid.arrange, c(err_plots, ncol=2))
  } else {
    .plotErrors(dq, nti, ntj, err_out, err_in, nominalQ, ...)
  }
}

#' @import ggplot2
.plotErrors <- function(dq, nti, ntj, err_out=TRUE, err_in=FALSE, nominalQ=FALSE, ...) {
  ACGT <- c("A", "C", "G", "T")
  tij <- 4*(which(ACGT==nti)-1) + which(ACGT==ntj)
  nij <- paste0(nti,"2",ntj)
  subdf <- as.data.frame(t(dq$trans))
  subdf$qave <- as.numeric(rownames(subdf))
  subdf[,nij] <- subdf[,nij]/(subdf[,paste0(nti,"2A")]+subdf[,paste0(nti,"2C")]+subdf[,paste0(nti,"2G")]+subdf[,paste0(nti,"2T")])
  
  p <- ggplot(data=subdf, aes_string(x="qave", y=paste0("log10(",nij,")")), ...)
  p <- p + geom_point(na.rm=TRUE)

  if(err_in) {
    erridf <- dq$err_in
    if(is.list(erridf)) erridf <- erridf[[1]]
    erridf <- as.data.frame(t(erridf))
    colnames(erridf) <- paste0(rep(ACGT, each=4), "2", ACGT)
    rownames(erridf) <- seq(dq$opts$QMIN, dq$opts$QMAX, length.out=nrow(erridf))
    erridf$qave <- as.numeric(rownames(erridf))
    p <- p + geom_line(data=erridf, linetype="dashed")
  }
  
  if(!("err_out" %in% ls(dq))) err_out <- FALSE
  if(err_out) {
    errodf <- as.data.frame(t(dq$err_out))
    errodf$qave <- as.numeric(rownames(errodf))
    p <- p + geom_line(data=errodf)
  }
  
  if(nominalQ) {
    qq <- subdf$qave
    if(nti != ntj) { 
      nom_err <- (1/3) * 10^(-qq/10)
    } else {
      nom_err <- 1. - 10^(-qq/10)
    }
    nomqdf <- data.frame(qave=qq, nom=nom_err)
    colnames(nomqdf)[[2]] <- nij
    p <- p + geom_line(data=nomqdf, color="red")
  }
  
  return(p)
}

#' @import ggplot2
#' @importFrom gridExtra grid.arrange
plotSubPos <- function(subpos, ...) {
  subpos$pos <- seq(nrow(subpos))
  subpos <- subpos[1:match(0,subpos$nts)-1,]
  p <- ggplot(data=subpos, aes(x=pos))
  pA <- p + geom_line(aes(y=A2C/(1+A)), color="red") + geom_line(aes(y=A2G/(1+A)), color="orange") + geom_line(aes(y=A2T/(1+A)), color="blue") + ylab("Subs at As")
  pC <- p + geom_line(aes(y=C2A/(1+C)), color="grey") + geom_line(aes(y=C2G/(1+C)), color="orange") + geom_line(aes(y=C2T/(1+C)), color="blue") + ylab("Subs at Cs")
  pG <- p + geom_line(aes(y=G2A/(1+G)), color="grey") + geom_line(aes(y=G2C/(1+G)), color="red") + geom_line(aes(y=G2T/(1+G)), color="blue") + ylab("Subs at Gs")
  pT <- p + geom_line(aes(y=T2A/(1+T)), color="grey") + geom_line(aes(y=T2C/(1+T)), color="red") + geom_line(aes(y=T2G/(1+T)), color="orange") + ylab("Subs at Ts")
  pAll <- p + geom_line(aes(y=subs/nts)) + ylab("Sub rate (all nts)")
  grid.arrange(pAll, pAll, pA, pC, pG, pT, nrow=3, ...)
}

