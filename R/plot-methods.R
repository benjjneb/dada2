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
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"), verbose = TRUE)
#' dada1 <- dada(derep1, err = inflateErr(tperr1, 2), selfConsist = TRUE) 
#' plotComplementarySubstitutions(dada1)
#' 
plotComplementarySubstitutions = function(dadaOut, facetByGrp = TRUE){
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

#' Plot error rates after dada processing.
#' 
#' This function plots the observed frequency of each transition
#' (eg. A->C) as a function of the associated quality score. It also plots the final
#' estimated error rates (if they exist). The initial input rates and the expected error
#' rates under the nominal definition of quality scores can also be included.
#' 
#' @param dq Required. A \code{\link{dada-class}} object.
#' 
#' @param nti Optional. Default c("A","C","G","T"). Some combination of the 4 DNA nucleotides.
#' @param ntj Optional. Default c("A","C","G","T"). Some combination of the 4 DNA nucleotides.
#'  The error rates from nti->ntj will be plotted. If multiple nti or ntj are chosen,
#'   error rates from each-to-each will be plotted in a grid.
#'  
#' @param err_out Optional. Default TRUE.
#'  A \code{logical(1)} determining whether to plot the output error rates (solid).
#' 
#' @param err_in Optional. Default FALSE.
#'  A \code{logical(1)} determining whether to plot the initial input error rates (dashed).
#'
#' @param nominalQ Optional. Default FALSE.
#'  A \code{logical(1)} determining whether to plot the expected error rate (red) if the
#'  quality score exactly matched its nominal definition: Q = -10 log10(p_err).
#'  
#' @return A \code{\link{ggplot}} object that will be rendered
#'  to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  
#' @importFrom reshape2 melt
#' @import ggplot2
# @importFrom ggplot2 ggplot
# @importFrom ggplot2 aes
# @importFrom ggplot2 geom_point
# @importFrom ggplot2 geom_line
# @importFrom ggplot2 facet_wrap
# @importFrom ggplot2 xlab
# @importFrom ggplot2 ylab
#' 
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"), verbose = TRUE)
#' dada1 <- dada(derep1, err = inflateErr(tperr1, 2), errorEstimationFunction = loessErrfun) 
#' plotErrors(dada1)
#' plotErrors(dada1, "A", "C")
#' plotErrors(dada1, nti="A", ntj=c("A","C","G","T"), err_in=TRUE, nominalQ=TRUE)
#' 
plotErrors <- function(dq, nti=c("A","C","G","T"), ntj=c("A","C","G","T"), err_out=TRUE, err_in=FALSE, nominalQ=FALSE) {
  ACGT <- c("A", "C", "G", "T")
  if(!(all(nti %in% ACGT) && all(ntj %in% ACGT)) || any(duplicated(nti)) || any(duplicated(ntj))) {
    stop("nti and ntj must be nucleotide(s): A/C/G/T.")
  }
  transdf = melt(dq$trans, factorsAsStrings=TRUE)
  colnames(transdf) <- c("Transition", "Qual", "count")
  transdf$from <- substr(transdf$Transition, 1, 1)
  transdf$to <- substr(transdf$Transition, 3, 3)
  tot.count <- tapply(transdf$count, list(transdf$from, transdf$Qual), sum)
  transdf$tot <- mapply(function(x,y) tot.count[x,y], transdf$from, as.character(transdf$Qual))
  transdf$Observed <- transdf$count/transdf$tot
  if(!is.null(dq$err_out)) {
    transdf$Estimated <- mapply(function(x,y) dq$err_out[x,y], transdf$Transition, as.character(transdf$Qual))
  } else {
    transdf$Estimated <- NA
  }
  # If selfConsist, then err_in is a list. Use the first err_in, the initial error rates provided.
  ei <- dq$err_in; if(is.list(ei)) ei <- ei[[1]]
  transdf$Input <- mapply(function(x,y) ei[x,y], transdf$Transition, as.character(transdf$Qual))
  transdf$Nominal <- (1/3)*10^-(transdf$Qual/10)
  transdf$Nominal[transdf$Transition %in% c("A2A", "C2C", "G2G", "T2T")] <- 1 - 10^-(transdf$Qual[transdf$Transition %in% c("A2A", "C2C", "G2G", "T2T")]/10)
  
  p <- ggplot(data=transdf[transdf$from %in% nti & transdf$to %in% ntj,], aes(x=Qual, y=log10(Observed))) + geom_point(na.rm=TRUE)
  if(err_out && !is.null(dq$err_out))  p <- p + geom_line(aes(y=log10(Estimated)))
  if(err_in)   p <- p + geom_line(aes(y=log10(Input)), linetype="dashed")
  if(nominalQ) p <- p + geom_line(aes(y=log10(Nominal)), color="red")
  p <- p + facet_wrap(~Transition, nrow=length(nti))
  p <- p + xlab("Consensus quality score") + ylab("Error frequency (log10)")
  p
}
