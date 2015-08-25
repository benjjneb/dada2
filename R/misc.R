################################################################################
#' Get the "uniques vector": Integer vector named by sequence and valued by abundance.
#' 
#' This function extracts the "uniques vector" from several different data objects, including
#'  \code{\link{dada-class}}, \code{\link{derep-class}} and \code{data.frame} objects that have both
#'  $sequence and $abundance columns.
#'  The return value is an integer vector named by sequence and valued by abundance. If input is already
#'  a "uniques vector", that same vector will be returned.
#'  The uniques format is used by several functions within the dada2 package.
#' 
#' @param object (Required). The object from which to extract the uniques vector.
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
getUniques <- function(object) {
  if(is.integer(object) && length(names(object)) != 0 && !any(is.na(names(object)))) { # Named integer vector already
    unqs <- object
  } else if(class(object) == "dada") {  # dada return 
    unqs <- object$denoised
  } else if(class(object) == "derep") {
    unqs <- object$uniques
  } else if(is.data.frame(object) && all(c("sequence", "abundance") %in% colnames(object))) {
    unqs <- as.integer(object$abundance)
    names(unqs) <- object$sequence
  } else {
    stop("Unrecognized format: Requires named integer vector, dada-class, derep-class, or a data.frame with $sequence and $abundance columns.")
  }
  if(any(duplicated(names(unqs)))) message("Duplicate sequences detected.")
  return(unqs)
}

################################################################################
#' Get vector of sequences.
#' 
#' This function extracts the unique sequences from several different data objects, including
#'  \code{\link{dada-class}}, \code{\link{derep-class}} and \code{data.frame} objects that have both
#'  $sequence and $abundance columns. This function wraps the \code{\link{getUniques}} function, but
#'  return only the names (i.e. the sequences).
#' 
#' @param object (Required). The object from which to extract the sequences.
#' 
#' @return \code{character}. A character vector of the sequences.
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
getSequences <- function(object) {
  return(names(getUniques(object)))
}

################################################################################
#' Needlman-Wunsch alignment with ends-free gapping.
#' 
#' This function performs a Needleman-Wunsch alignment with free gaps at the ends of the sequences.
#' 
#' @param s1 (Required). \code{character(1)}. The first sequence to align. A/C/G/T only.
#' 
#' @param s2 (Required). \code{character(1)}. The second sequence to align. A/C/G/T only.
#' 
#' @param score (Optional). A 4x4 numeric matrix.
#'  The transition scores to use for the alignment. Default is getDadaOpt("SCORE_MATRIX").
#' 
#' @param gap (Optional). A \code{numeric(1)}. The gap penalty.
#'  Default is getDadaOpt("GAP_PENALTY").
#'  
#' @param band (Optional). A \code{numeric(1)}. The band size.
#'  This Needleman-Wunsch alignment is banded. This value specifies the radius of that band.
#'  Default is getDadaOpt("BAND_SIZE").
#'  
#' @return \code{character(2)}. The aligned sequences.
#' 
#' @export
#' 
#' @examples
#'  sq1 <- "CTAATACATGCAAGTCGAGCGAGTCTGCCTTGAAGATCGGAGTGCTTGCACTCTGTGAAACAAGATACAGGCTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCCAAGAGATCGGGATAACACCTGGAAACAGATGCTAATACCGGATAACAACAGATGATGCCTATCAACTGTTTAAAAGATGGTTCTGCTATCACTCTTGGATGGACCTGCG"
#'  sq2 <- "TTAACACATGCAAGTCGAACGGAAAGGCCAGTGCTTGCACTGGTACTCGAGTGGCGAACGGGTGAGTAACACGTGGGTGATCTGCCCTGTACTTCGGGATAAGCTTGGGAAACTGGGTCTAATACCGGATAGGACAACTTTTTGGATATTGTTGTGGAAAGCTTTTGCGGTATGGGATGAGCTCGCGGCCTATCAGCTTGTTGGTGGGGTAATGGCCTAC"
#'  nwalign(sq1, sq2)
#'  nwalign(sq1, sq2, band=-1)
#' 
nwalign <- function(s1, s2, score=getDadaOpt("SCORE_MATRIX"), gap=getDadaOpt("GAP_PENALTY"), band=getDadaOpt("BAND_SIZE")) {
  if(!is.character(s1) || !is.character(s2)) stop("Can only align character sequences.")
  if(nchar(s1) >= 1000 || nchar(s2) >= 1000) stop("Can only align strings up to 999 nts in length.")
  if(nchar(s1) != nchar(s2)) {
    if(band != -1) message("Sequences of unequal length must use unbanded alignment.")
    band = -1
  }
  C_nwalign(s1, s2, score, gap, band)
}

################################################################################
#' Hamming distance after Needlman-Wunsch alignment with ends-free gapping.
#' 
#' This function performs a Needleman-Wunsch alignment with free gaps at the ends of the sequences, and then counts
#' the number of mismatches and indels in that alignment. Gaps at the beginning and end are ignored.
#' 
#' @param s1 (Required). \code{character(1)}. The first sequence to align. A/C/G/T only.
#' 
#' @param s2 (Required). \code{character(1)}. The second sequence to align. A/C/G/T only.
#' 
#' @param ... (Optional). Further arguments to pass on to \code{\link{nwalign}}.
#' 
#' @return \code{integer(1)}. The total number of mismatches and gaps, excluding gaps at the beginning and end of the alignment.
#' 
#' @export
#' 
#' @examples
#' sq1 <- "CTAATACATGCAAGTCGAGCGAGTCTGCCTTGAAGATCGGAGTGCTTGCACTCTGTGAAACAAGATACAGGCTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCCAAGAGATCGGGATAACACCTGGAAACAGATGCTAATACCGGATAACAACAGATGATGCCTATCAACTGTTTAAAAGATGGTTCTGCTATCACTCTTGGATGGACCTGCG"
#' sq2 <- "TTAACACATGCAAGTCGAACGGAAAGGCCAGTGCTTGCACTGGTACTCGAGTGGCGAACGGGTGAGTAACACGTGGGTGATCTGCCCTGTACTTCGGGATAAGCTTGGGAAACTGGGTCTAATACCGGATAGGACAACTTTTTGGATATTGTTGTGGAAAGCTTTTGCGGTATGGGATGAGCTCGCGGCCTATCAGCTTGTTGGTGGGGTAATGGCCTAC"
#' nwhamming(sq1, sq2)
#' nwhamming(sq1, sq2, band=-1)
#' 
nwhamming <- Vectorize(function(s1, s2, ...) {
  al <- nwalign(s1, s2, ...)
  out <- C_eval_pair(al[1], al[2])
  return(out["mismatch"]+out["indel"])
})

nweval <- Vectorize(function(s1, s2, ...) {
  al <- nwalign(s1, s2, ...)
  C_eval_pair(al[1], al[2])
})

strdiff <- function(s1, s2) {
  xx = unlist(strsplit(s1,""))
  yy = unlist(strsplit(s2,""))
  dd <- which(xx != yy)
  data.frame(pos=dd,nt0=xx[dd],nt1=yy[dd])
}

hamming <- Vectorize(function(x, y) nrow(strdiff(x, y)))

#' @importFrom Biostrings DNAString
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings reverseComplement
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

