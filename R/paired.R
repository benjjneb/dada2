#' Function to merge DADA-denoised paired end reads
#' @export
#' 
mergePairs <- function(dqF, mapF, dqR, mapR, minOverlap = 20) {
  rF <- dqF$map[mapF]
  rR <- dqR$map[mapR]
  if(any(is.na(rF)) || any(is.na(rR))) stop("Non-corresponding maps and dada-outputs.")
  
  pairdf <- data.frame(forward=rF, reverse=rR)
  ups <- unique(pairdf) # The unique forward/reverse pairs of denoised sequences
  Funqseq <- dqF$clustering$sequence[ups$forward]
  Runqseq <- sapply(dqR$clustering$sequence[ups$reverse], rc)
    
  Fstart <- unname(mapply(regexpr, subseq(Runqseq,1,minOverlap), Funqseq))
  Fstart[Fstart==-1] <- 1 # Should cause matching step to come up unequal
  Fend <- nchar(Funqseq)
  Rstart <- rep(1, length(Fstart))
  Rend <- Rstart + Fend - Fstart
  Rend[Rend > nchar(Runqseq)] <- nchar(Runqseq)[Rend > nchar(Runqseq)]
  # Should come up unequal in next step
  
  ups$match <- mapply(function(x,y) x==y, subseq(Funqseq, Fstart, Fend), subseq(Runqseq, Rstart, Rend))

  tab <- table(pairdf$forward, pairdf$reverse)
  ups$abundance <- tab[cbind(ups$forward, ups$reverse)]
  ups$sequence <- paste0(Funqseq, subseq(Runqseq,Rend+1,nchar(Runqseq)))
  rownames(ups) <- paste0("s", ups$forward, "_", ups$reverse)
  ups
}




