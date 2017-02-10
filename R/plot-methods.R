#' Plot Substitution Pairs from DADA Result
#' 
#' This is similar to original DADA article, Figure 6.
#' 
#' @param dadaOut (Required). A \code{\link{dada-class}} object.
#' 
#' @param facetByGrp (Optional). Default TRUE.
#'  Whether to plot all substitution groups together in one panel
#'  or separately on a grid of panels with a linear model fit.
#' 
#' @return A \code{\link{ggplot}2} object.
#'  Will be rendered to default device if \code{\link{print}ed},
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

#' Plot observed and estimated error rates.
#' 
#' This function plots the observed frequency of each transition
#' (eg. A->C) as a function of the associated quality score. It also plots the final
#' estimated error rates (if they exist). The initial input rates and the expected error
#' rates under the nominal definition of quality scores can also be shown.
#' 
#' @param dq (Required). An object from which error rates can be extracted. Valid inputs are
#'  coercible by \code{\link{getErrors}}. This includes the output of the \code{\link{dada}}
#'  and \code{\link{learnErrors}} functions.
#' 
#' @param nti (Optional). Default c("A","C","G","T"). 
#'  Some combination of the 4 DNA nucleotides.
#'  
#' @param ntj (Optional). Default c("A","C","G","T"). 
#'  Some combination of the 4 DNA nucleotides.
#'  
#'  The error rates from nti->ntj will be plotted. If multiple nti or ntj are chosen,
#'   error rates from each-to-each will be plotted in a grid.
#'  
#' @param obs (Optional). Default TRUE.
#'  If TRUE, the observed error rates are plotted as points.  
#'  
#' @param err_out (Optional). Default TRUE.
#'  If TRUE, plot the output error rates (solid line).
#' 
#' @param err_in (Optional). Default FALSE.
#'  If TRUE, plot the input error rates (dashed line).
#'
#' @param nominalQ (Optional). Default FALSE.
#'  If TRUE, plot the expected error rates (red line) if quality scores
#'  exactly matched their nominal definition: Q = -10 log10(p_err).
#'  
#' @return A \code{\link{ggplot}2} object.
#'  Will be rendered to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  See \code{\link{ggsave}} for additional options.
#'  
#' @seealso 
#'  \code{\link{learnErrors}}, \code{\link{getErrors}}
#'
#' @importFrom reshape2 melt
#' @importFrom methods is
#' @import ggplot2
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
plotErrors <- function(dq, nti=c("A","C","G","T"), ntj=c("A","C","G","T"), obs=TRUE, err_out=TRUE, err_in=FALSE, nominalQ=FALSE) {
  ACGT <- c("A", "C", "G", "T")
  if(!(all(nti %in% ACGT) && all(ntj %in% ACGT)) || any(duplicated(nti)) || any(duplicated(ntj))) {
    stop("nti and ntj must be nucleotide(s): A/C/G/T.")
  }
  
  dq <- getErrors(dq, detailed=TRUE, enforce=FALSE)
  
  if(!is.null(dq$trans)) {
    if(ncol(dq$trans) <= 1) {
      stop("plotErrors only supported when using quality scores in the error model (i.e. USE_QUALS=TRUE).")
    }
    transdf = melt(dq$trans, factorsAsStrings=TRUE)
    colnames(transdf) <- c("Transition", "Qual", "count")
  } else if(!is.null(dq$err_out)) {
    if(ncol(dq$err_out) <= 1) {
      stop("plotErrors only supported when using quality scores in the error model (i.e. USE_QUALS=TRUE).")
    }
    transdf = melt(dq$err_out, factorsAsStrings=TRUE)
    colnames(transdf) <- c("Transition", "Qual", "est")
  } else {
    stop("Non-null observed and/or estimated error rates (dq$trans or dq$err_out) must be provided.")
  }
  transdf$from <- substr(transdf$Transition, 1, 1)
  transdf$to <- substr(transdf$Transition, 3, 3)
  
  if(!is.null(dq$trans)) {
    tot.count <- tapply(transdf$count, list(transdf$from, transdf$Qual), sum)
    transdf$tot <- mapply(function(x,y) tot.count[x,y], transdf$from, as.character(transdf$Qual))
    transdf$Observed <- transdf$count/transdf$tot
  } else {
    transdf$Observed <- NA
    obs <- FALSE
  }
  if(!is.null(dq$err_out)) {
    transdf$Estimated <- mapply(function(x,y) dq$err_out[x,y], transdf$Transition, as.character(transdf$Qual))
  } else {
    transdf$Estimated <- NA
    err_out <- FALSE
  }
  if(!is.null(dq$err_in)) {
    # If selfConsist, then err_in is a list. Use the first err_in, the initial error rates provided.
    ei <- dq$err_in; if(is.list(ei)) ei <- ei[[1]]
    transdf$Input <- mapply(function(x,y) ei[x,y], transdf$Transition, as.character(transdf$Qual))
  } else {
    transdf$Input <- NA
    err_in <- FALSE
  }
  transdf$Nominal <- (1/3)*10^-(transdf$Qual/10)
  transdf$Nominal[transdf$Transition %in% c("A2A", "C2C", "G2G", "T2T")] <- 1 - 10^-(transdf$Qual[transdf$Transition %in% c("A2A", "C2C", "G2G", "T2T")]/10)
  
  p <- ggplot(data=transdf[transdf$from %in% nti & transdf$to %in% ntj,], aes(x=Qual))
  if(obs) p <- p + geom_point(aes(y=Observed), na.rm=TRUE)
  if(err_out)  p <- p + geom_line(aes(y=Estimated))
  if(err_in)   p <- p + geom_line(aes(y=Input), linetype="dashed")
  if(nominalQ) p <- p + geom_line(aes(y=Nominal), color="red")
  p <- p + scale_y_log10()
  p <- p + facet_wrap(~Transition, nrow=length(nti))
  p <- p + xlab("Consensus quality score") + ylab("Error frequency (log10)")
  p
}

#' Plot quality profile of a fastq file.
#' 
#' This function plots a visual summary of the distribution of quality scores
#' as a function of sequence position for the input fastq file.
#' 
#' The distribution of quality scores at each position is shown as a grey-scale
#' heat map, with dark colors corresponding to higher frequency. The plotted lines
#' show positional summary statistics: green is the mean, orange is the median, and
#' the dashed orange lines are the 25th and 75th quantiles.
#' 
#' @param fl (Required). \code{character(1)}.
#'  The file path to the fastq or fastq.gz file.
#' 
#' @param n (Optional). Default 1,000,000.
#'  The number of records to sample from the fastq file.
#' 
#' @return A \code{\link{ggplot}2} object.
#'  Will be rendered to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  See \code{\link{ggsave}} for additional options.
#'  
#' @importFrom ShortRead qa
#' @import ggplot2
#' 
#' @export
#' 
#' @examples
#' plotQualityProfile(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' 
# This code is adapted from ShortRead:::.plotCycleQuality
#
plotQualityProfile <- function(fl, n=1000000) {
  srqa <- qa(fl, n=n)
  df <- srqa[["perCycle"]]$quality
  rc <- srqa[["readCounts"]]$read
  if (rc >= n){
    rclabel <- paste("Readcounts >= ", n)
  } else {
    rclabel <- paste("Readcounts: ", rc)
  }
  # Calculate summary statistics at each position
  means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
  get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >=q)][[1]] }
  q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
  q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
  q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
  if(!all(sapply(list(names(q25s), names(q50s), names(q75s)), identical, rownames(means)))) {
    stop("Calculated quantiles/means weren't compatible.")
  }
  statdf <- data.frame(Cycle=as.integer(rownames(means)),
                       Mean=means, 
                       Q25=as.vector(q25s), Q50=as.vector(q50s), Q75=as.vector(q75s))
  # Create plot
  ggplot(data=df, aes(x=Cycle, y=Score)) + geom_tile(aes(fill=Count)) + 
    scale_fill_gradient(low="#F5F5F5", high="black") + 
    geom_line(data=statdf, aes(y=Mean), color="#66C2A5") +
    geom_line(data=statdf, aes(y=Q25), color="#FC8D62", size=0.25, linetype="dashed") +
    geom_line(data=statdf, aes(y=Q50), color="#FC8D62", size=0.25) +
    geom_line(data=statdf, aes(y=Q75), color="#FC8D62", size=0.25, linetype="dashed") +
    ylab("Quality Score") + xlab("Cycle") +
    theme_bw() + theme(panel.grid=element_blank()) + guides(fill=FALSE) +
    annotate("text", 1, min(df$Score), label = basename(fl), hjust=0, vjust=0) +
    annotate("text", 1, min(df$Score), label = rclabel, hjust=0, vjust=2)
}
