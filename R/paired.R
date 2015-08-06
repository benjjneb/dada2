################################################################################
#' Merge paired forward and reverse reads after DADA denoising.
#' 
#' This function takes the output of dada() and the map from read to unique index
#' returned by derepFastq()$map for both the forward and reverse data set. It attempts
#' to merge each pair of reads, rejecting any which do not perfectly overlap.
#' Note: This function assumes that the fastq files for the forward and reverse reads
#' were in the same order.
#' 
#' @param dadaF (Required). Output of dada() function.
#'  The output list returned by the dada() function on the forward reads.
#' 
#' @param mapF (Required). An integer vector map from read index to unique index.
#'  The map returned by derepFastq() for the forward reads.
#'   
#' @param dadaR (Required). Output of dada() function.
#'  The output list returned by the dada() function on the reverse reads.
#' 
#' @param mapR (Required). An integer vector map from read index to unique index.
#'  The map returned by derepFastq() for the reverse reads.
#'
#' @param keepMismatch(Optional). A \code{logical(1)}. Default is False.
#'  If true, the pairs that did not match are retained in the return data.frame.
#'
#' @param minOverlap (Optional). A \code{numeric(1)} of the minimum overlap
#'  required for merging the forward and reverse reads. Default is 20.
#'
#' @param propagateCol (Optional). Character vector. Default is empty.
#'  The mergePairs return data.frame will include copies of columns with names specified
#'  in the dada()$clustering data.frame.
#'
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#'
#' @return Dataframe. One row for each unique pairing of forward/reverse denoised sequences.
#' \itemize{
#'  \item{$abundance: }{Number of reads corresponding to this forward/reverse combination.}
#'  \item{$sequence: The merged sequence if match=TRUE. Otherwise an empty sequence, i.e. "".}
#'  \item{$forward: The index of the forward denoised genotypes.}
#'  \item{$reverse: The index of the reverse denoised genotype.}
#'  \item{$nmatch: Number of matching nts in the overlap region.}
#'  \item{$nmismatch: Number of mismatching nts in the overlap region.}
#'  \item{$nindel: Number of indels in the overlap region.}
#'  \item{$match: TRUE if a perfect match between the forward and reverse denoised sequences of at least MIN_OVERLAP.
#'          FALSE otherwise.}
#'  \item{$...: Additional columns specified in the propagateCol argument.}
#' }
#' 
#' @seealso \code{\link{derepFastq}}, \code{\link{dada}}
#'
#' @export
#' @import Biostrings 
#' 
mergePairs <- function(dadaF, mapF, dadaR, mapR, keepMismatch=FALSE, minOverlap = 20, propagateCol=character(0), verbose=TRUE, align=TRUE) {
  rF <- dadaF$map[mapF]
  rR <- dadaR$map[mapR]
  if(any(is.na(rF)) || any(is.na(rR))) stop("Non-corresponding maps and dada-outputs.")
  
  pairdf <- data.frame(sequence = "", abundance=0, forward=rF, reverse=rR)
  ups <- unique(pairdf) # The unique forward/reverse pairs of denoised sequences
  Funqseq <- unname(dadaF$clustering$sequence[ups$forward])
  Runqseq <- as(reverseComplement(DNAStringSet(unname(dadaR$clustering$sequence[ups$reverse]))), "character")
  
  # Use unbanded N-W align to compare forward/reverse
  # May want to adjust align params here, but for now just using dadaOpt
  alvecs <- mapply(function(x,y) nwalign(x,y,band=-1), Funqseq, Runqseq, SIMPLIFY=FALSE)
  outs <- t(sapply(alvecs, function(x) C_eval_pair(x[1], x[2])))
  ups$nmatch <- outs[,1]
  ups$nmismatch <- outs[,2]
  ups$nindel <- outs[,3]
  ups$match <- (ups$nmatch > minOverlap) & (ups$nmismatch==0) & (ups$nindel==0)
  # Make the sequence
  ups$sequence <- sapply(alvecs, function(x) C_pair_consensus(x[[1]], x[[2]]));

  # Add abundance and sequence to the output data.frame
  tab <- table(pairdf$forward, pairdf$reverse)
  ups$abundance <- tab[cbind(ups$forward, ups$reverse)]
  ups$sequence[!ups$match] <- ""
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
    cat(sum(ups$abundance[ups$match]), "paired-reads (in", sum(ups$match), "unique pairings) successfully merged out of", sum(ups$abundance), "(in", nrow(ups), "pairings) input.\n")
  }
  
  if(!keepMismatch) { ups <- ups[ups$match,] }
  ups
}

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


