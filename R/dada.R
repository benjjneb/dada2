#' The dada function takes as input the unique sequences in an amplicon sequencing sample paired
#'  with their abundances, and returns the sample genotypes paired with their abundances
#'  inferred by the Divisive Amplicon Denoising Algorithm (Rosen, Callahan, Fisher, Holmes 2012).
#'
#' @param uniques (Required). Named integer vector.
#'  This is a vector of abundances of unique sequences, with the sequence specified by the name.
#'  Sequences are only allowed to contain A/C/G/T/N/- characters.
#'  Warning: N and - need more testing!
#'  
#' @param err (Optional). 4x4 numeric matrix.
#'  The matrix of estimated error rates from one nucleotide to another.
#'  err[i,j] = Prob(j in sequence | i in sample genotype).
#'  A=1, C=2, G=3, T=4 for indexing. Rows required to sum to 1.
#'  
#' @param score (Optional). 4x4 numeric matrix.
#'  The score matrix used by the Needleman-Wunsch alignment algorithm. Defaults to nuc44.
#'  
#' @param gap_penalty (Optional). \code{numeric(1)}
#'  The gap penalty for the Needleman-Wunsch alignemnt algoirithm. Defaults to -8.
#'  
#' @param self_consist (Optional). \code{logical(1)}
#'  When true, DADA will re-estimate error rates after inferring sample genotypes, and then repeat
#'  the algorithm using the newly estimated error rates. This continues until convergence.
#'
#' @return List.
#'  $genotypes: named integer vector of the inferred sample genotypes.
#'  $trans: 4x4 integer matrix of the observed transitions between nts. 
#'   
#' @export
#'
dada <- function(uniques,
                 err = matrix(c(0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991), nrow=4, byrow=T),
                 score = matrix(c(5, -4, -4, -4, -4, 5, -4, -4, -4, -4, 5, -4, -4, -4, -4, 5), nrow=4, byrow=T),
                 gap_penalty = -8,
                 self_consist = FALSE) {
  # Validate inputs
  if(!(is.integer(uniques))) {
    if(is.numeric(uniques) && all.equal(uniques, as.integer(uniques))) {
      nms <- names(uniques)
      uniques <- as.integer(uniques)
      names(uniques) <- nms
    } else {
      stop("dada: Invalid uniques vector. Must be integer valued.")
    }
  }
  if(!(all(sapply(names(bar), function(x) nchar(gsub("[ACGTN-]", "", x))==0, USE.NAMES=FALSE)))) {
    stop("dada: Invalid uniques vector. Names must be sequences made up of A/C/G/T/N/-")
  }
  
  if(!( is.numeric(err) && dim(err) == c(4,4) && all(err>=0) && all.equal(rowSums(err), c(1,1,1,1)) )) {
    stop("dada: Invalid error matrix.")
  }
  if(any(err==0)) warning("dada: Zero in error matrix.")

  if(!(is.numeric(score) && dim(err) == c(4,4))) {
    stop("dada: Invalid score matrix.")
  }
  
  if(!(is.numeric(gap_penalty) && gap_penalty <=0)) {
    stop("dada: Invalid gap penalty.")
  }
  if(gap_penalty > -1) warning("dada: Very small gap penalty.")
  
  prev <- NULL
  nn <- 1
  cur <- dada_uniques(names(uniques), unname(uniques), err, score, gap_penalty,
                      get("USE_KMERS", envir=dada_opts), get("KDIST_CUTOFF", envir=dada_opts))
  while(self_consist && !identical(cur, prev) && nn < get("MAX_CONSIST", envir=dada_opts)) {
    cat(".")
    prev <- cur
    err <- cur$trans + 1   # ADD ONE PSEUDOCOUNT TO EACH TRANSITION
    err <- t(apply(err, 1, function(x) x/sum(x)))  # apply returns a transposed result
    cur <- dada_uniques(names(uniques), unname(uniques), err, score, gap_penalty,
                        get("USE_KMERS", envir=dada_opts), get("KDIST_CUTOFF", envir=dada_opts))
    nn <- nn+1
  }
  cat("\n")
  if(self_consist && nn == get("MAX_CONSIST", envir=dada_opts)) {
    warning("dada: Self-consistency loop terminated before convergence.")
  }
  # Convert cur$genotypes to the named integer vector being used as the uniques format
  foo <- as.integer(cur$genotypes$abundance)
  names(foo) <- cur$genotypes$sequence
  cur$genotypes <- sort(foo, decreasing=TRUE)
  cur
}
