#' The dada function takes as input the unique sequences in an amplicon sequencing sample paired
#'  with their abundances, and returns the sample genotypes paired with their abundances
#'  inferred by the Divisive Amplicon Denoising Algorithm (Rosen, Callahan, Fisher, Holmes 2012).
#'
#' @param uniques (Required). Named integer vector or list of named integer vectors.
#'  Each uniques vector is an integer vector of abundances, named by the unique DNA
#'    sequence to which that abundance corresponds.
#'  Sequences are only allowed to contain A/C/G/T characters.
#'  A list of such vectors can be provided, in which case the error model will be 
#'    shared across these samples when denoising.
#'  
#' @param err (Optional). 4x4 numeric matrix.
#'  The matrix of estimated error rates from one nucleotide to another.
#'  err[i,j] = Prob(j in sequence | i in sample genotype).
#'  A=1, C=2, G=3, T=4 for indexing. Rows required to sum to 1.
#'  
#' @param self_consist (Optional). \code{logical(1)}
#'  When true, DADA will re-estimate error rates after inferring sample genotypes, and then repeat
#'  the algorithm using the newly estimated error rates. This continues until convergence.
#'
#' @return List.
#'  $genotypes: named integer vector of the denoised sample genotypes.
#'  $trans: 4x4 integer matrix of the inferred substutions ("errors") between nts. 
#'  $opts: A list of the dada_opts used for this function call.
#'   
#' @export
#'
dada <- function(uniques,
                 err = matrix(c(0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991), nrow=4, byrow=T),
                 self_consist = FALSE) {
  
  call <- sys.call(1)
  if(!is.list(uniques)) { # If a single vector, make into a length 1 list
    uniques <- list(uniques)
  }
  
  # Validate uniques vector(s)
  for(i in seq(length(uniques))) {
    if(!(is.integer(uniques[[i]]))) {
      if(is.numeric(uniques[[i]]) && all.equal(uniques[[i]], as.integer(uniques[[i]]))) {
        nms <- names(uniques[[i]])
        uniques[[i]] <- as.integer(uniques[[i]])
        names(uniques[[i]]) <- nms
      } else {
        stop("dada: Invalid uniques vector. Must be integer valued.")
      }
    }
    if(!(all(sapply(names(uniques[[i]]), function(x) nchar(gsub("[ACGT]", "", x))==0, USE.NAMES=FALSE)))) {
      stop("dada: Invalid uniques vector. Names must be sequences made up only of A/C/G/T.")
    }
  }
  
  if(!( is.numeric(err) && dim(err) == c(4,4) && all(err>=0) && all.equal(rowSums(err), c(1,1,1,1)) )) {
    stop("dada: Invalid error matrix.")
  }
  if(any(err==0)) warning("dada: Zero in error matrix.")

  prev <- NULL
  clustering <- list()
  trans <- matrix(0, nrow=4, ncol=4)
  nconsist <- 1
  for(i in seq(length(uniques))) {
    res <- dada_uniques(names(uniques[[i]]), unname(uniques[[i]]), err, 
                        get("SCORE_MATRIX", envir=dada_opts), get("GAP_PENALTY", envir=dada_opts),
                        get("USE_KMERS", envir=dada_opts), get("KDIST_CUTOFF", envir=dada_opts),
                        get("BAND_SIZE", envir=dada_opts),
                        get("OMEGA_A", envir=dada_opts), 
                        get("USE_SINGLETONS", envir=dada_opts), get("OMEGA_S", envir=dada_opts))
    clustering[[i]] <- res$clustering
    trans <- trans + res$trans
  }
  cur = list(clustering=clustering, trans=trans)
    
  while(self_consist && !identical(cur, prev) && nconsist < get("MAX_CONSIST", envir=dada_opts)) {
    cat(".")
    prev <- cur
    clustering <- list()
    trans <- matrix(0, nrow=4, ncol=4)
    err <- cur$trans + 1   # ADD ONE PSEUDOCOUNT TO EACH TRANSITION
    err <- t(apply(err, 1, function(x) x/sum(x)))  # apply returns a transposed result
    for(i in seq(length(uniques))) {
      res <- dada_uniques(names(uniques[[i]]), unname(uniques[[i]]), err, 
                          get("SCORE_MATRIX", envir=dada_opts), get("GAP_PENALTY", envir=dada_opts),
                          get("USE_KMERS", envir=dada_opts), get("KDIST_CUTOFF", envir=dada_opts),
                          get("BAND_SIZE", envir=dada_opts),
                          get("OMEGA_A", envir=dada_opts), 
                          get("USE_SINGLETONS", envir=dada_opts), get("OMEGA_S", envir=dada_opts))
      clustering[[i]] <- res$clustering
      trans <- trans + res$trans  
    }
    cur = list(clustering=clustering, trans=trans)
    nconsist <- nconsist+1
  }
  cat("\n")
  if(self_consist && nconsist == get("MAX_CONSIST", envir=dada_opts)) {
    warning("dada: Self-consistency loop terminated before convergence.")
  }
  
  # Construct dada return object
  # Convert cur$genotypes to the named integer vector being used as the uniques format
  rval = list()
  rval$genotypes <- list()
  rval$clustering <- list()
  for(i in seq(length(uniques))) {
    rval$genotypes[[i]] <- as.integer(cur$clustering[[i]]$abundance)
    names(rval$genotypes[[i]]) <- cur$clustering[[i]]$sequence
    rval$clustering <- cur$clustering
  }
  if(length(rval$genotypes)==1) { # one sample, return a naked uniques vector
    rval$genotypes <- rval$genotypes[[1]]
    rval$clustering <- rval$clustering[[1]]
  } else { # keep names if it is a list
    names(rval$genotypes) <- names(uniques)
    names(rval$clustering) <- names(uniques)
  }

  rval$trans <- cur$trans
  rval$err <- err
  
  # Store all the options in the return object
  opts <- ls(dada_opts)
  ropts <- lapply(opts, function(x) get(x, envir=dada_opts))
  names(ropts) <- opts
  rval$opts <- ropts
  
  rval$call <- call
  
  rval
}
