# Run the auto-installer/updater code:
source("install.R", local = TRUE)


library("shiny")
library("shinyFiles")

# Number of available cores
message("Number of available cores:\n")
(NCores <- max(1L, RcppParallel::defaultNumThreads()))

################################################################################
#' Tabulate quality by cycle for a random subset of reads in a fastq file
#' 
#' @param fastqFile Valid path to a fastq (or fastq.gz) file.
#' @param nReads The number of reads to sample from this file.
#' 
#' @return A tidy data.table of position, quality, frequency (counts).
#' 
#' @import ShortRead
#' @import data.table
tabulate_quality = function(fastqFile, nReads = 1e4){
  require("ShortRead")
  require("data.table")
  message("Reading from file: ", fastqFile, "\n")
  # fastqFile = fls[1]
  FQS = FastqSampler(con = fastqFile, n = nReads)
  fq = yield(FQS)
  # Fun borrowed from ShortRead internals
  .qa_perCycleQuality = function(abc, quality){
    if (missing(abc) || dim(abc)[[3]] == 0) {
      df <- data.frame(Cycle=integer(0), Quality=numeric(0),
                       Score=numeric(0), Count=integer(0),
                       lane=character(0))
      return(df)
    }
    abc <- apply(abc, 2:3, sum)
    q <- factor(rownames(abc)[row(abc)], levels=rownames(abc))
    q0 <- as(do.call(class(quality), list(rownames(abc))), "matrix")
    df <- data.frame(Cycle=as.integer(colnames(abc)[col(abc)]),
                     Quality=q, Score=as.integer(q0)[q],
                     Count=as.vector(abc),
                     row.names=NULL)
    df[df$Count != 0, ]
  }
  abc <- alphabetByCycle(fq)
  perCycleQuality <- data.table(.qa_perCycleQuality(abc, quality(fq)))
  perCycleQuality[, fastqFile := fastqFile]
  return(perCycleQuality)
}
################################################################################
#' Wrapper for running DADA2 algorithm from sequence file to sequence result.
#' 
#' Check that this isn't redundant with recent additions. Migrate if so.
#'
wrap_dada2_workflow = function(seqFiles,
                               dadaOutFiles = c("DADA2-Forward.RDS", "DADA2-Reverse.RDS"),
                               err = NULL,
                               selfConsist = TRUE,
                               minOverlap = 20,
                               maxMismatch = 0,
                               # Performance params
                               nReads = 1e6,
                               multithread = TRUE){
  require("dada2")
  merged = NULL
  
  stopifnot(all(file.exists(seqFiles)))
  seqFileF = seqFiles[1]
  seqFileR = seqFiles[2]
  
  stopifnot(
    dir.exists(
      dirname(
        c(dadaOutFiles[1],
          dadaOutFiles[2]))))
  
  dadaOutFileF = dadaOutFiles[1]
  dadaOutFileR = dadaOutFiles[2]
  
  if(length(err) > 2){
    warning("Provided more than two error matrices. Most likely something is wrong.")
  }
  
  if(length(err) == 2){
    # Assume in forward-then-reverse order
    errF = err[[1]]
    errR = err[[2]]
  } else {
    errF = errR = err[[1]]
  }
  dadaF = dadaR = derepF = derepR = NULL
  
  message("Dereplicating forward reads:\n", seqFileF, "\n")
  derepF <- derepFastq(seqFileF, n = nReads)
  message("DADA2-ing:\n", seqFileF, "\n")
  dadaF <- dada(derepF, err=errF, selfConsist = selfConsist, multithread = multithread)
  saveRDS(dadaF, dadaOutFileF)
  
  message("Dereplicating reverse reads:\n", seqFileR, "\n")
  derepR <- derepFastq(seqFileR, n = nReads)
  message("DADA2-ing:\n", seqFileR, "\n")
  dadaR <- dada(derepR, err=errR, selfConsist = selfConsist, multithread = multithread)
  saveRDS(dadaR, dadaOutFileR)
  
  message("Merging DADA2 results for read-pairs:\n",
          paste0(seqFiles, collapse = "\n"), "\n")
  # merger <- mergePairs(ddF, derepF, ddR, derepR)
  trash = try(expr = {
    ## Merge seq directions by ID, return an abundance data.table
    merged = dada2:::mergePairsByID(
      # Forward
      dadaF = dadaF,
      derepF = derepF,
      srF = seqFileF,
      # Reverse
      dadaR = dadaR,
      derepR = derepR,
      srR = seqFileR,
      # Additional params
      minOverlap = minOverlap,
      maxMismatch = maxMismatch,
      returnRejects = FALSE,
      verbose = TRUE)
    message("Sum of read-pairs properly merged after denoising:\n",
            round(100 * merged[, sum(abundance[(accept)])/sum(abundance)], digits = 1), "%")
    merged <- merged[(accept & !is.na(sequence)), .(sequence, abundance)]
    saveRDS(merged, file = file.path(dirname(dadaOutFileF), 
                                     gsub("\\.RDS$", "merge.RDS", basename(dadaOutFileF))))
  }, silent = TRUE)
  return(merged)
}
################################################################################
#' Plot sequence quality by cycle (position in sequence)
#' 
#' 
plot_quality_by_cycle = function(CycleStats, CycleCounts, TrimTable = NULL){
  p = ggplot(data = CycleStats,
             mapping = aes(Cycle, y = Quality, color = Statistic)) +
    ylim(0, 40) +
    geom_raster(data = CycleCounts,
                mapping = aes(x = Cycle,
                              y = Score,
                              fill = log10(Proportion)),
                inherit.aes = FALSE) +
    # grey0 (black) to grey100 (white)
    # scale_fill_gradient(
    scale_fill_gradient2(
      midpoint = -1,
      low = "grey99",
      mid = "grey90",
      high = "grey10",
      space = "white",
      na.value = "white",
      guide = "colourbar") +
    geom_path(size = 0.25) +
    geom_path(mapping = aes(y = Smooth), size = 1, alpha = 0.35) +
    geom_text(mapping = aes(label = Cycle,
                            y = 1,
                            hjust = ifelse(Side == "Right", yes = 1.1, no = -0.1)),
              data = TrimTable,
              vjust = 0.5,
              color = "black", size = 3) +
    facet_wrap(~Direction, nrow = 2)
  if( !is.null(TrimTable) ){
    p <- p + geom_vline(mapping = aes(xintercept = Cycle), data = TrimTable, size = 0.25)
  }
  return(p)
}
################################################################################
#' Write tab delimited table
#' 
#' @inheritParams utils::write.table
#' @param ... additional args passed to [write.table]
#' 
write_table_tab = function(x, file, ...){
  write.table(
    x = x, 
    file = file,
    sep = "\t", 
    row.names = FALSE, 
    col.names = TRUE,
    append = FALSE,
    quote = FALSE,
    ...) 
}
################################################################################

