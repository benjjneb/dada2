################################################################################
#' Merge paired forward and reverse reads after DADA denoising.
#' 
#' This function attempts each denoised pair of forward and reverse reads, rejecting any 
#' which do not sufficiently overlap or which contain too many (>0 by default) mismatches in the
#' overlap region. Note: This function assumes that the fastq files for the 
#' forward and reverse reads were in the same order. Use of the concatenate option will
#' result in concatenating forward and reverse reads without attempting a merge/alignment step.
#' 
#' @param dadaF (Required). A \code{\link{dada-class}} object.
#'  The output of dada() function on the forward reads.
#' 
#' @param derepF (Required). A \code{\link{derep-class}} object.
#'  The derep-class object returned by derepFastq() that was used as the input to the
#'  dada-class object passed to the dadaF argument.
#'  
#'  Alternatively the map itself (derepFastq()$map) can be provided in place of the
#'  derep-class object.
#'   
#' @param dadaR (Required). A \code{\link{dada-class}} object.
#'  The output of dada() function on the reverse reads.
#' 
#' @param derepR (Required). A \code{\link{derep-class}} object.
#'  See derepF description, but for the reverse reads.
#'
#' @param minOverlap (Optional). A \code{numeric(1)} of the minimum length of the overlap (in nucleotides)
#'  required for merging the forward and reverse reads. Default is 20.
#'
#' @param maxMismatch (Optional). A \code{numeric(1)} of the maximum mismatches allowed in the overlap region.
#'  Default is 0 (i.e. only exact matches in the overlap region are accepted).
#'  
#' @param returnRejects (Optional). A \code{logical(1)}. Default is False.
#'  If true, the pairs that that were rejected based on mismatches in the overlap
#'  region are retained in the return data.frame.
#'
#' @param propagateCol (Optional). \code{character}. Default is \code{character(0)}.
#'  The mergePairs return data.frame will include copies of columns with names specified
#'  in the dada-class$clustering data.frame.
#'
#' @param justConcatenate (Optional). \code{logical(1)}, Default FALSE.
#'  If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged,
#'    with a NNNNNNNNNN (10 Ns) spacer inserted between them.
#' 
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#'
#' @return Dataframe. One row for each unique pairing of forward/reverse denoised sequences.
#' \itemize{
#'  \item{$abundance: }{Number of reads corresponding to this forward/reverse combination.}
#'  \item{$sequence: The merged sequence.}
#'  \item{$forward: The index of the forward denoised sequence.}
#'  \item{$reverse: The index of the reverse denoised sequence.}
#'  \item{$nmatch: Number of matching nts in the overlap region.}
#'  \item{$nmismatch: Number of mismatching nts in the overlap region.}
#'  \item{$nindel: Number of indels in the overlap region.}
#'  \item{$prefer: The sequence used for the overlap region. 1=forward; 2=reverse.}
#'  \item{$accept: TRUE if overlap between forward and reverse denoised sequences was at least minOverlap and had at most maxMismatch differences.
#'          FALSE otherwise.}
#'  \item{$...: Additional columns specified in the propagateCol argument.}
#' }
#' 
#' @seealso \code{\link{derepFastq}}, \code{\link{dada}}
#' @export
#' 
#' @examples
#' derepF = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' derepR = derepFastq(system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#' dadaF <- dada(derepF, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' dadaR <- dada(derepR, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' mergePairs(dadaF, derepF, dadaR, derepR)
#' mergePairs(dadaF, derepF, dadaR, derepR, returnRejects=TRUE, propagateCol=c("n0", "birth_ham"))
#' mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE)
#' 
mergePairs <- function(dadaF, derepF, dadaR, derepR, minOverlap = 20, maxMismatch=0, returnRejects=FALSE, propagateCol=character(0), justConcatenate=FALSE, verbose=FALSE) {
  if(class(derepF) == "derep") mapF <- derepF$map
  else mapF <- derepF
  if(class(derepR) == "derep") mapR <- derepR$map
  else mapR <- derepR
  if(!(is.integer(mapF) && is.integer(mapR))) stop("Incorrect format of derep arguments.")
  rF <- dadaF$map[mapF]
  rR <- dadaR$map[mapR]
  if(any(is.na(rF)) || any(is.na(rR))) stop("Non-corresponding maps and dada-outputs.")
  
  pairdf <- data.frame(sequence = "", abundance=0, forward=rF, reverse=rR)
  ups <- unique(pairdf) # The unique forward/reverse pairs of denoised sequences
  Funqseq <- unname(as.character(dadaF$clustering$sequence[ups$forward]))
  Runqseq <- rc(unname(as.character(dadaR$clustering$sequence[ups$reverse])))
  
  if (justConcatenate == TRUE) {
    # Simply concatenate the sequences together
    ups$sequence <- mapply(function(x,y) paste0(x,"NNNNNNNNNN", y), Funqseq, Runqseq, SIMPLIFY=FALSE);  
    ups$nmatch <- 0
    ups$nmismatch <- 0
    ups$nindel <- 0
    ups$prefer <- NA
    ups$accept <- TRUE
  } else {
    # Align forward and reverse reads.
    # Use unbanded N-W align to compare forward/reverse
    # May want to adjust align params here, but for now just using dadaOpt
    alvecs <- mapply(function(x,y) nwalign(x,y,band=-1), Funqseq, Runqseq, SIMPLIFY=FALSE)
    outs <- t(sapply(alvecs, function(x) C_eval_pair(x[1], x[2])))
    ups$nmatch <- outs[,1]
    ups$nmismatch <- outs[,2]
    ups$nindel <- outs[,3]
    ups$prefer <- 1 + (dadaR$clustering$n0[ups$reverse] > dadaR$clustering$n0[ups$forward])
    ups$accept <- (ups$nmatch > minOverlap) & ((ups$nmismatch + ups$nindel) <= maxMismatch)
    # Make the sequence
    ups$sequence <- mapply(C_pair_consensus, sapply(alvecs,`[`,1), sapply(alvecs,`[`,2), ups$prefer);
    # Additional param to indicate whether 1:forward or 2:reverse takes precedence
    # Must also strip out any indels in the return
    # This function is only used here.
  }
  
  # Add abundance and sequence to the output data.frame
  tab <- table(pairdf$forward, pairdf$reverse)
  ups$abundance <- tab[cbind(ups$forward, ups$reverse)]
  ups$sequence[!ups$accept] <- ""
  # Add columns from forward/reverse clustering
  propagateCol <- propagateCol[propagateCol %in% colnames(dadaF$clustering)]
  for(col in propagateCol) {
    ups[,paste0("F.",col)] <- dadaF$clustering[ups$forward,col]
    ups[,paste0("R.",col)] <- dadaR$clustering[ups$reverse,col]
  }
  # Sort output by abundance and name
  ups <- ups[order(ups$abundance, decreasing=TRUE),]
  rownames(ups) <- paste0("s", ups$forward, "_", ups$reverse)
  if(verbose) {
    message(sum(ups$abundance[ups$accept]), " paired-reads (in ", sum(ups$accept), " unique pairings) successfully merged out of ", sum(ups$abundance), " (in ", nrow(ups), " pairings) input.")
  }
  if(!returnRejects) { ups <- ups[ups$accept,] }
  
  if(any(duplicated(ups$sequence))) {
    message("Duplicate sequences in merged output.")
  }
  
  ups
}

#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead id
#' @importFrom ShortRead yield
sameOrder <- function(fnF, fnR) {
  matched <- TRUE
  fF <- FastqStreamer(fnF)
  on.exit(close(fF))
  fR <- FastqStreamer(fnR)
  on.exit(close(fR), add=TRUE)
  
  while( length(suppressWarnings(fqF <- yield(fF))) && length(suppressWarnings(fqR <- yield(fR))) ) {
    idF <- trimTails(id(fqF), 1, " ")
    idR <- trimTails(id(fqR), 1, " ")
    matched <- matched && all(idF == idR)
  }
  return(matched)
}

isMatch <- function(al, minOverlap, verbose=FALSE) {
  out <- C_eval_pair(al[1], al[2]) # match, mismatch, indel
  if(verbose) { cat("Match/mismatch/indel:", out, "\n") }
  if(out[1] >= minOverlap && out[2] == 0 && out[3] == 0) {
    return(TRUE);
  } else {
    return(FALSE);
  }
}
