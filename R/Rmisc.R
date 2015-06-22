#' @export
isHit100 <- function(clust, fn) {
  bb <- read.table(fn, comment.char="#", col.names=c("seqid", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score"))
  bbHit100 <- bb[bb$identity == 100 & bb$coverage == nchar(clust[match(bb$seqid,clust$id),"sequence"]),]
  return(clust$id %in% bbHit100$seqid)
}

#' @export
isOneOff <- function(clust, fn) {
  bb <- read.table(fn, comment.char="#", col.names=c("seqid", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score"))
  bb <- bb[bb$coverage == nchar(clust[match(bb$seqid,clust$id),"sequence"]),] # Only full length hits
  tab <- tapply(bb$identity, bb$seqid, max)
  tab <- tab[match(clust$id, names(tab))]
  seqlens <- nchar(clust$sequence)
  oneOff <- tab<100 & (abs(tab - 100.0*(seqlens-1)/seqlens) < 0.01)
  oneOff[is.na(oneOff)] <- FALSE # happens if no hits were full coverage
  names(oneOff) <- clust$id # Also drop the name to NA so fix here
  return(oneOff)
}

#' @export
checkConvergence <- function(dadaO) {
  sapply(dadaO$err_in, function(x) sum(abs(dadaO$err_out-x)))
}

#' @export
nwalign <- function(s1, s2, score=getDadaOpt("SCORE_MATRIX"), gap=getDadaOpt("GAP_PENALTY"), band=getDadaOpt("BAND_SIZE")) {
  if(nchar(s1) != nchar(s2)) {
    if(band != -1) message("Sequences of unequal length must use unbanded alignment.")
    band = -1
  }
  C_nwalign(s1, s2, score, gap, band)
}

#' @export 
nwhamming <- Vectorize(function(s1, s2, ...) {
  al <- nwalign(s1, s2, ...)
  out <- dadac:::C_eval_pair(al[1], al[2])
  return(out["mismatch"])
})

#' @export 
nweval <- function(s1, s2, ...) {
  al <- nwalign(s1, s2, ...)
  C_eval_pair(al[1], al[2])
}

#' @export
strdiff <- function(s1, s2) {
  xx = unlist(strsplit(s1,""))
  yy = unlist(strsplit(s2,""))
  dd <- which(xx != yy)
  data.frame(pos=dd,nt0=xx[dd],nt1=yy[dd])
}

#' @export
rc <- function(x) as(reverseComplement(DNAString(x)), "character")

#' @export
hamming <- Vectorize(function(x, y) nrow(strdiff(x, y)))

#' @export
#' @import ggplot2
#' @import gridExtra
showSubPos <- function(subpos, ...) {
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

#' @export
getIll <- function(fn, remove_singletons = FALSE) {
  if(is.numeric(fn)) {
    fn <- paste0("~/Desktop/Illumina/metaID-", fn, "_R1.fastq.gz")
  }
  ill <- derepFastq(fn)
  ill <- filterNs(ill)
  if(remove_singletons) { ill <- ill[ill>1] }
  ill
}

#' @export
filterNs <- function(unqs) { # For now all Ns must be removed
  acgts <- sapply(names(unqs), function(x) nchar(gsub("[ACGT]", "", x))==0)
  unqs[acgts]  
}

#' @export
subseqUniques <- function(unqs, start, end) {
  subnms <- subseq(names(unqs), start, end)
  newNames <- unique(subnms)
  newUniques <- as.integer(rep(0,length(newNames)))
  names(newUniques) <- newNames
  for(i in seq(length(unqs))) {
    newnm <- subnms[[i]]
    newUniques[[newnm]] <- newUniques[[newnm]] + unqs[[i]]
  }
  newUniques[sapply(names(newUniques), function(nm) nchar(nm) == (end-start+1))]
}

#' @export
mergeUniques <- function(unqsList, ...) {
  if(!is.list(unqsList) && length(list(...))>=1) {
    unqsList = list(unqsList, unlist(unname(list(...))))
  }
  concat <- c(unlist(unqsList))
  unqs <- unique(names(concat))
  mrg <- as.integer(rep(0, length(unqs)))
  names(mrg) <- unqs
  # Probably a better way to do this than for loop...
  for(i in seq(length(concat))) {
    unq <- names(concat)[[i]]
    mrg[[unq]] <- mrg[[unq]] + concat[[i]]
  }
  mrg
}

#' @export
as.uniques <- function(df) {
  if(is.integer(df) && length(names(df)) != 0 && !any(is.na(names(df)))) { # Named integer vector already
    return(df)
  } else if(all(c("genotypes", "clustering") %in% names(df))) {  # dada return 
    return(df$genotypes)
  } else if(is.data.frame(df) && all(c("sequence", "abundance") %in% colnames(df))) {
    unqs <- as.integer(df$abundance)
    names(unqs) <- df$sequence
    return(unqs)
  } else {
    stop("Unrecognized format: Requires named integer vector, dada object or $clustering data.frame.")
  }
}
