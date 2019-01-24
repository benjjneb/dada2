################################################################################
#' Get the uniques-vector from the input object.
#' 
#' This function extracts the \code{\link{uniques-vector}} from several different data objects, 
#'  including \code{\link{dada-class}} and \code{\link{derep-class}} objects, as well as 
#'  \code{data.frame} objects that have both $sequence and $abundance columns.
#'  The return value is an integer vector named by sequence and valued by abundance. If the input is
#'  already in \code{\link{uniques-vector}} format, that same vector will be returned.
#' 
#' @param object (Required). The object from which to extract the \code{\link{uniques-vector}}.
#' 
#' @param collapse (Optional). Default TRUE.
#'  Should duplicate sequences detected in \code{object} be collapsed together, thereby
#'   imposing uniqueness on non-unique input.
#'  
#' @param silence (Optional). Default FALSE.
#'  Suppress reporting of the detection and merger of duplicated input sequences.
#' 
#' @return \code{integer}.
#'  An integer vector named by unique sequence and valued by abundance.
#' 
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1)
#' getUniques(derep1)
#' getUniques(dada1)
#' getUniques(dada1$clustering)
#' 
getUniques <- function(object, collapse=TRUE, silence=FALSE) {
  if(is.integer(object) && length(names(object)) != 0 && !any(is.na(names(object)))) { # Named integer vector already
    unqs <- object
  } else if(class(object) == "dada") {  # dada return 
    unqs <- object$denoised
  } else if(class(object) == "derep") {
    unqs <- object$uniques
  } else if(is.data.frame(object) && all(c("sequence", "abundance") %in% colnames(object))) {
    unqs <- as.integer(object$abundance)
    names(unqs) <- object$sequence
  } else if(class(object) == "matrix" && is.numeric(object) && !any(is.na(colnames(object)))) { # Tabled sequences
    unqs <- as.integer(colSums(object))
    names(unqs) <- colnames(object)
  }
  else {
    stop("Unrecognized format: Requires named integer vector, dada-class, derep-class, sequence matrix, or a data.frame with $sequence and $abundance columns.")
  }
  #### ENFORCE UNIQUENESS HERE!!!
  if(any(duplicated(names(unqs)))) {
    if(collapse) {
      unqs <- tapply(unqs, names(unqs), sum)
      if(!silence) message("Duplicate sequences detected and merged.")
    } else if(!silence) {
      message("Duplicate sequences detected.")
    }
  }
  return(unqs)
}

################################################################################
#' Get vector of sequences from input object.
#' 
#' This function extracts the sequences from several different data objects, including
#'  including \code{\link{dada-class}} and \code{\link{derep-class}} objects, as well as 
#'  \code{data.frame} objects that have both $sequence and $abundance columns. This function 
#'  wraps the \code{\link{getUniques}} function, but return only the names (i.e. the sequences).
#'  Can also be provided the file path to a fasta or fastq file, a taxonomy table, or a
#'  DNAStringSet object. 
#' 
#' @param object (Required). The object from which to extract the sequences.
#' 
#' @param collapse (Optional). Default FALSE.
#'  Should duplicate sequences detected in \code{object} be collapsed together, thereby
#'  imposing uniqueness on non-unique input.
#'  
#' @param silence (Optional). Default TRUE.
#'  Suppress reporting of the detection and merger of duplicated input sequences.
#' 
#' @return \code{character}. A character vector of the sequences.
#' 
#' @importFrom methods is
#' @importFrom methods as
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead readFastq
#' @importFrom ShortRead sread
#' @importFrom ShortRead id
#' 
#' @export
#' 
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' dada1 <- dada(derep1, err=tperr1)
#' getSequences(derep1)
#' getSequences(dada1)
#' getSequences(dada1$clustering)
#' 
getSequences <- function(object, collapse=FALSE, silence=TRUE) {
  if(is(object, "character")) {
    if(length(object)==1 && file.exists(object)) {
      sr <- tryCatch(readFasta(object), error=function(err) { readFastq(object) })
      seqs <- toupper(as.character(sread(sr)))
      names(seqs) <- id(sr)
      return(seqs)
    } else if(collapse) {
      if(any(duplicated(object)) && !silence) message("Duplicate sequences detected and merged.")
      return(unique(object))
    } else {
      return(object)
    }
  } else if(class(object) == "DNAStringSet") {
    return(as.character(object))
  } else if(class(object) == "matrix" && is.character(object) && !any(is.na(rownames(object)))) { # Taxonomy table
    seqs <- rownames(object)
    if(any(duplicated(seqs))) {
      if(collapse) seqs <- unique(seqs)
      if(collapse && !silence) message("Duplicate sequences detected and merged.")
      if(!collapse && !silence) message("Duplicate sequences detected.")
    }
    return(seqs)
  } else {
    return(names(getUniques(object, collapse=collapse, silence=silence)))
  }
}

getAbund <- function(object) {
  return(sum(getUniques(object)))
}

getNseq <- function(object) {
  return(length(getUniques(object)))
}

################################################################################
#' Needleman-Wunsch alignment.
#' 
#' This function performs a Needleman-Wunsch alignment between two sequences.
#' 
#' @param s1 (Required). \code{character(1)}. The first sequence to align. A/C/G/T only.
#' 
#' @param s2 (Required). \code{character(1)}. The second sequence to align. A/C/G/T only.
#' 
#' @param match (Optional). \code{numeric(1)}. Default is getDadaOpt("MATCH").
#'  The score of a match in the alignment.
#' 
#' @param mismatch (Optional). \code{numeric(1)}. Default is getDadaOpt("MISMATCH").
#'  The score of a mismatch in the alignment.
#' 
#' @param gap (Optional). \code{numeric(1)}. Default is getDadaOpt("GAP_PENALTY").
#'  The alignment gap penalty. Should be negative.
#'  
#' @param homo_gap (Optional). \code{numeric(1)}. Default NULL (no special homopolymer penalty).
#'  The alignment gap penalty within homopolymer regions. Should be negative.
#'  
#' @param band (Optional). \code{numeric(1)}.  Default -1 (no banding).
#'  The Needleman-Wunsch alignment can be banded. This value specifies the radius of that band.
#'  Set \code{band = -1} to turn off banding.
#'  
#' @param endsfree (Optional). \code{logical(1)}. Default TRUE.
#'  Allow unpenalized gaps at the ends of the alignment.
#'  
#' @param vec (Optional). \code{logical(1)}. Default FALSE.
#'  Use DADA2's vectorized aligner instead of standard DP matrix. Not intended for long sequences (>1kb).
#'  
#' @return \code{character(2)}. The aligned sequences.
#' 
#' @export
#' 
#' @examples
#'  sq1 <- "CTAATACATGCAAGTCGAGCGAGTCTGCCTTGAAGATCGGAGTGCTTGCACTCTGTGAAACAAGATA"
#'  sq2 <- "TTAACACATGCAAGTCGAACGGAAAGGCCAGTGCTTGCACTGGTACTCGAGTGGCGAACGGGTGAGT"
#'  nwalign(sq1, sq2)
#'  nwalign(sq1, sq2, band=16)
#' 
nwalign <- function(s1, s2, match=getDadaOpt("MATCH"), mismatch=getDadaOpt("MISMATCH"), gap=getDadaOpt("GAP_PENALTY"), homo_gap=NULL, band=-1, endsfree=TRUE, vec=FALSE) {
  if(!is.character(s1) || !is.character(s2)) stop("Can only align character sequences.")
  if(is.null(homo_gap)) { homo_gap <- gap }
  if(vec) {
    if(homo_gap != gap) stop("Homopolymer gap penalties are not implemented in the vectorized aligner.")
    return(C_nwvec(s1, s2, match, mismatch, gap, band, endsfree))
  } else {
    if(!C_isACGT(s1) || !C_isACGT(s2)) {
      stop("Sequences must contain only A/C/G/T characters.")
    }
    return(C_nwalign(s1, s2, match, mismatch, gap, homo_gap, band, endsfree))
  }
}

################################################################################
#' Hamming distance after Needlman-Wunsch alignment.
#' 
#' This function performs a Needleman-Wunsch alignment between two sequences, and then counts
#' the number of mismatches and indels in that alignment. End gaps are not included in this count.
#' 
#' @param s1 (Required). \code{character(1)}. The first sequence to align. A/C/G/T only.
#' 
#' @param s2 (Required). \code{character(1)}. The second sequence to align. A/C/G/T only.
#' 
#' @param ... (Optional). Further arguments to pass on to \code{\link{nwalign}}.
#' 
#' @return \code{integer(1)}. The total number of mismatches and gaps, excluding gaps at the beginning
#'  and end of the alignment.
#' 
#' @export
#' 
#' @examples
#'  sq1 <- "CTAATACATGCAAGTCGAGCGAGTCTGCCTTGAAGATCGGAGTGCTTGCACTCTGTGAAACAAGATA"
#'  sq2 <- "TTAACACATGCAAGTCGAACGGAAAGGCCAGTGCTTGCACTGGTACTCGAGTGGCGAACGGGTGAGT"
#' nwhamming(sq1, sq2)
#' nwhamming(sq1, sq2, band=16)
#' 
nwhamming <- Vectorize(function(s1, s2, ...) {
  al <- nwalign(s1, s2, ...)
  out <- C_eval_pair(al[1], al[2])
  return(unname(out["mismatch"]+out["indel"]))
}, USE.NAMES=FALSE)

nweval <- Vectorize(function(s1, s2, ...) {
  al <- nwalign(s1, s2, ...)
  C_eval_pair(al[1], al[2])
}, USE.NAMES=FALSE)

nwextract <- function(query, ref, ...) {
  al <- nwalign(query, ref, ...)
  ntq <- gregexpr("[ACGT]", al[[1]])
  rval <- substr(al[[2]], min(ntq[[1]]), max(ntq[[1]]))
  rval <- gsub("-", "", rval)
  rval
}

strdiff <- function(s1, s2) {
  xx = unlist(strsplit(s1,""))
  yy = unlist(strsplit(s2,""))
  dd <- which(xx != yy)
  data.frame(pos=dd,nt0=xx[dd],nt1=yy[dd])
}

#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
#' @importFrom methods as
rc <- function(sqs) {
  if(length(sqs) < 1) {
    return(character(0))
  } else if(length(sqs) == 1) {
    as(reverseComplement(DNAString(sqs)), "character")
  } else {
    as(reverseComplement(DNAStringSet(sqs)), "character")
  }
}

checkConvergence <- function(dadaO) {
  sapply(dadaO$err_in, function(x) sum(abs(dadaO$err_out-x)))
}

pfasta <- function(seqs, ids=seq(length(seqs))) {
  seqs <- getSequences(seqs, collapse=FALSE)
  cat(paste(">", ids, "\n", seqs, sep="", collapse="\n"))
}

#' @importFrom methods is
#' @keywords internal
is.list.of <- function(x, ctype) {
  if(!is.list(x)) return(FALSE)
  else return(all(sapply(x, is, ctype)))
}

#' @keywords internal
seqtab_to_qiime <- function(st, fout) {
  st <- t(st) # QIIME has OTUs as rows
  col.names <- colnames(st)
  col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
  write.table(st, fout, sep="\t",
              row.names=TRUE, col.names=col.names, quote=FALSE)
}

#' @keywords internal
seqtab_to_mothur <- function(st, fout) {
  # mothur has OTUs as columns, and a couple required columns
  df.shared <- data.frame(label=rep("DADA2", nrow(st)), Group=rownames(st), numOtus=ncol(st))
  df.shared <- cbind(df.shared, st)
  write.table(df.shared, four, row.names=FALSE, col.names=TRUE, quote=FALSE)
}

#' @keywords internal
samdf_to_qiime2 <- function(df, fout) {
  col.names <- colnames(df)
  col.names[[1]] <- paste0("#SampleID\t", col.names[[1]])
  write.table(df, fout, sep="\t",
              row.names=TRUE, col.names=col.names, quote=FALSE)
}

#' @keywords internal
bs1ham <- function(dd, ham=1) {
  is.1ham <- which(dd$clustering$birth_ham %in% ham)
  dd$birth_subs[dd$birth_subs$clust %in% is.1ham,]
}

#' @keywords internal
getSRR <- Vectorize(function(run, outdir="sra", verbose=TRUE, ...) {
  if(!grepl("^SRR[0-9]{6+}$", run)) stop("Requires SRA Run accessions in format: SRR1234567")
  if(!dir.exists(outdir)) dir.create(outdir)
  loc <- paste0("ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/",
                substr(run, 1, 6), "/", run, "/")
  loc <- paste0(loc, run, ".sra")
  download.file(loc, file.path(outdir, paste0(run, ".sra")), ...)
  if(verbose) cat(run, "\n")
})



