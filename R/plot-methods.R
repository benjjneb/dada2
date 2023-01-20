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
  if(obs) p <- p + geom_point(aes(y=Observed), color="gray40", na.rm=TRUE)
  if(err_in)   p <- p + geom_line(aes(y=Input), linetype="dashed")
  if(err_out)  p <- p + geom_line(aes(y=Estimated))
  if(nominalQ) p <- p + geom_line(aes(y=Nominal), color="red")
  p <- p + scale_y_log10()
  p <- p + facet_wrap(~Transition, nrow=length(nti))
  p <- p + xlab("Consensus quality score") + ylab("Error frequency (log10)")
  p <- p + theme_bw()
  p
}

#' Plot quality profile of a fastq file.
#' 
#' This function plots a visual summary of the distribution of quality scores
#' as a function of sequence position for the input fastq file(s).
#' 
#' The distribution of quality scores at each position is shown as a grey-scale
#' heat map, with dark colors corresponding to higher frequency. The plotted lines
#' show positional summary statistics: green is the mean, orange is the median, and
#' the dashed orange lines are the 25th and 75th quantiles.
#' 
#' If the sequences vary in length, a red line will be plotted showing the percentage
#' of reads that extend to at least that position.
#' 
#' @param fl (Required). \code{character}.
#'  File path(s) to fastq or fastq.gz file(s).
#' 
#' @param n (Optional). Default 500,000.
#'  The number of records to sample from the fastq file.
#' 
#' @param aggregate (Optional). Default FALSE.
#'  If TRUE, compute an aggregate quality profile for all fastq files provided.
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
plotQualityProfile <- function(fl, n=500000, aggregate=FALSE) {
  statdf <- data.frame(Cycle=integer(0), Mean=numeric(0), Q25=numeric(0), Q50=numeric(0), Q75=numeric(0), Cum=numeric(0), file=character(0))
  anndf <- data.frame(minScore=numeric(0), label=character(0), rclabel=character(0), rc=numeric(0), file=character(0))

  FIRST <- TRUE
  for(f in fl[!is.na(fl)]) {
    srqa <- qa(f, n=n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read) # Handle aggregate form from qa of a directory
    if (rc >= n) { 
      rclabel <- paste("Reads >= ", n)
    } else {
      rclabel <- paste("Reads: ", rc)
    }
    # Calculate summary statistics at each position
    means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
    get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >=q)][[1]] }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), simplify=TRUE)
    if(!all(sapply(list(names(q25s), names(q50s), names(q75s), names(cums)), identical, rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if(FIRST) {
      plotdf <- cbind(df, file=basename(f))
      FIRST <- FALSE
    } else { plotdf <- rbind(plotdf, cbind(df, file=basename(f))) }
    statdf <- rbind(statdf, data.frame(Cycle=as.integer(rownames(means)), Mean=means, 
                                       Q25=as.vector(q25s), Q50=as.vector(q50s), Q75=as.vector(q75s), Cum=10*as.vector(cums)/min(rc, n), file=basename(f)))
    anndf <- rbind(anndf, data.frame(minScore=min(df$Score), label=basename(f), rclabel=rclabel, rc=rc, file=basename(f)))
  }
  anndf$minScore <- min(anndf$minScore)
  # Create plot
  if (aggregate) {
  	plotdf.summary <- aggregate(Count ~ Cycle + Score, plotdf, sum)
  	plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
  	means <- rowsum(plotdf.summary$Score*plotdf.summary$Count, plotdf.summary$Cycle)/rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
    q25s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
    q50s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
    q75s <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
    cums <- by(plotdf.summary, plotdf.summary$Cycle, function(foo) sum(foo$Count), simplify=TRUE)
    statdf.summary <- data.frame(Cycle=as.integer(rownames(means)), Mean=means, Q25=as.vector(q25s), Q50=as.vector(q50s), Q75=as.vector(q75s), Cum=10*as.vector(cums)/sum(pmin(anndf$rc, n)))
    p <- ggplot(data=plotdf.summary, aes(x=Cycle, y=Score)) + geom_tile(aes(fill=Count)) + 
		  scale_fill_gradient(low="#F5F5F5", high="black") + 
		  geom_line(data=statdf.summary, aes(y=Mean), color="#66C2A5") +
		  geom_line(data=statdf.summary, aes(y=Q25), color="#FC8D62", linewidth=0.25, linetype="dashed") +
		  geom_line(data=statdf.summary, aes(y=Q50), color="#FC8D62", linewidth=0.25) +
		  geom_line(data=statdf.summary, aes(y=Q75), color="#FC8D62", linewidth=0.25, linetype="dashed") +
		  ylab("Quality Score") + xlab("Cycle") + 
		  annotate("text", x=0, y=0, label=sprintf("Total reads: %d", sum(anndf$rc)), color="red", hjust=0) + 
		  theme_bw() + theme(panel.grid=element_blank()) + guides(fill="none") + 
      facet_wrap(~label)
    if(length(unique(statdf$Cum))>1) {
      p <- p + geom_line(data=statdf.summary, aes(y=Cum), color="red", linewidth=0.25, linetype="solid") +
        scale_y_continuous(limits = c(0,NA), sec.axis=sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) + 
        theme(axis.text.y.right = element_text(color = "red"), axis.title.y.right = element_text(color = "red"))
    } else {
      p <- p + ylim(c(0,NA))
    }
  } else {
  	p <- ggplot(data=plotdf, aes(x=Cycle, y=Score)) + geom_tile(aes(fill=Count)) + 
		  scale_fill_gradient(low="#F5F5F5", high="black") + 
		  geom_line(data=statdf, aes(y=Mean), color="#66C2A5") +
		  geom_line(data=statdf, aes(y=Q25), color="#FC8D62", linewidth=0.25, linetype="dashed") +
		  geom_line(data=statdf, aes(y=Q50), color="#FC8D62", linewidth=0.25) +
      geom_line(data=statdf, aes(y=Q75), color="#FC8D62", linewidth=0.25, linetype="dashed") +
      ylab("Quality Score") + xlab("Cycle") +
		  theme_bw() + theme(panel.grid=element_blank()) + guides(fill="none") +
		  geom_text(data=anndf, aes(x=0, label=rclabel, y=0), color="red", hjust=0) + 
      facet_wrap(~file)
    if(length(unique(statdf$Cum))>1) {
      p <- p + geom_line(data=statdf, aes(y=Cum), color="red", linewidth=0.25, linetype="solid") +
        scale_y_continuous(limits = c(0,NA), sec.axis=sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) + 
        theme(axis.text.y.right = element_text(color = "red"), axis.title.y.right = element_text(color = "red"))
    } else {
      p <- p + ylim(c(0,NA))
    }
  }
  p
}

#' Plot sequence complexity profile of a fastq file.
#' 
#' This function plots a histogram of the distribution of sequence complexities
#' in the form of effective numbers of kmers as determined by \code{\link{seqComplexity}}.
#' By default, kmers of size 2 are used, in which case a perfectly random sequences
#' will approach an effective kmer number of 16 = 4 (nucleotides) ^ 2 (kmer size).
#' 
#' @param fl (Required). \code{character}.
#'  File path(s) to fastq or fastq.gz file(s).
#' 
#' @param kmerSize (Optional). Default 2.
#'  The size of the kmers (or "oligonucleotides" or "words") to use.
#'
#' @param window (Optional). Default NULL.
#' The width in nucleotides of the moving window. If NULL the whole sequence is used.
#'
#' @param by (Optional). Default 5.
#' The step size in nucleotides between each moving window tested.
#'
#' @param n (Optional). Default 100,000.
#'  The number of records to sample from the fastq file.
#' 
#' @param bins (Optional). Default 100.
#'  The number of bins to use for the histogram.
#'
#' @param aggregate (Optional). Default FALSE.
#'  If TRUE, compute an aggregate quality profile for all fastq files provided.
#'  
#' @param ... (Optional). Arguments passed on to \code{\link{geom_histogram}}.
#' 
#' @return A \code{\link{ggplot}2} object.
#'  Will be rendered to default device if \code{\link{print}ed},
#'  or can be stored and further modified.
#'  See \code{\link{ggsave}} for additional options.
#'  
#' @importFrom ShortRead FastqSampler
#' @importFrom ShortRead yield
#' @import ggplot2
#' 
#' @seealso
#'  \code{\link{seqComplexity}}
#'  \code{\link[Biostrings]{oligonucleotideFrequency}}
#'
#' @export
#' 
#' @examples
#' plotComplexity(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' 
plotComplexity <- function(fl, kmerSize=2, window=NULL, by=5, n=100000, bins=100, aggregate=FALSE, ...) {
  cmplx <- lapply(fl, function(fli) {
    f <- FastqSampler(fli, n)
    srq <- yield(f)
    close(f)
    seqComplexity(sread(srq), kmerSize=kmerSize, window=window, by=by) # Also seqlens for warning?
  })
  df <- data.frame( complexity=unlist(cmplx), 
                    file=rep(basename(fl), times=sapply(cmplx, length)) )
  if(aggregate) { df$file <- paste(length(fl), "files (aggregated)") }
  p <- ggplot(data=df, aes(x=complexity)) + geom_histogram(bins=bins, na.rm=TRUE, ...) +
    ylab("Count") + xlab("Effective Oligonucleotide Number") +
    theme_bw() +
    facet_wrap(~file) + 
    scale_x_continuous(limits=c(0, 4^kmerSize), breaks=seq(0, 4^kmerSize, (4^kmerSize)/4))
  p
}


