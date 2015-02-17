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
                 self_consist = FALSE, ...) {
  
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
        stop("Invalid uniques vector. Must be integer valued.")
      }
    }
    if(!(all(sapply(names(uniques[[i]]), function(x) nchar(gsub("[ACGT]", "", x))==0, USE.NAMES=FALSE)))) {
      stop("Invalid uniques vector. Names must be sequences made up only of A/C/G/T.")
    }
  }
  
  if(!( is.numeric(err) && dim(err) == c(4,4) && all(err>=0) && all.equal(rowSums(err), c(1,1,1,1)) )) {
    stop("Invalid error matrix.")
  }
  if(any(err==0)) warning("Zero in error matrix.")
  
  # Read in default opts and then replace with any that were passed in to the function
  opts <- get_dada_opt()
  args <- list(...)
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) {
      opts[[opnm]] <- args[[opnm]]
    } else {
      warning(opnm, " is not a valid DADA option.")
    }
  }

  # Initialize
  cur <- NULL
  nconsist <- 1
  errs <- list()
  # The main loop, run once, or repeat until error rate repeats if self_consist=T
  repeat{
    clustering <- list()
    substitutions <- list()
    trans <- matrix(0, nrow=4, ncol=4)
    prev <- cur
    errs[[nconsist]] <- err

    for(i in seq(length(uniques))) {
      cat("Sample", i, "-", sum(uniques[[i]]), "reads in", length(uniques[[i]]), "unique sequences.\n")
      res <- dada_uniques(names(uniques[[i]]), unname(uniques[[i]]), err, 
                          opts[["SCORE_MATRIX"]], opts[["GAP_PENALTY"]],
                          opts[["USE_KMERS"]], opts[["KDIST_CUTOFF"]],
                          opts[["BAND_SIZE"]],
                          opts[["OMEGA_A"]], 
                          opts[["USE_SINGLETONS"]], opts[["OMEGA_S"]])
      clustering[[i]] <- res$clustering
      substitutions[[i]] <- res$substitutions
      trans <- trans + res$trans
    }
    cur = trans # The only thing that changes is err which is set by trans, so this is sufficient
    
    if((!self_consist) || identical(cur, prev) || (nconsist >= get("MAX_CONSIST", envir=dada_opts))) {
      break
    }
    
    err <- trans + 1   # ADD ONE PSEUDOCOUNT TO EACH TRANSITION
    err <- t(apply(err, 1, function(x) x/sum(x)))  # apply returns a transposed result
    nconsist <- nconsist+1
    cat(".")
  } # repeat

  cat("\n")
  if(self_consist) {
    if(nconsist == get("MAX_CONSIST", envir=dada_opts)) {
      warning("dada: Self-consistency loop terminated before convergence.")
    } else {
      cat("\nConvergence after ", nconsist, " rounds.\n")
    }
  }
  
  # Construct dada return object
  # Convert cur$genotypes to the named integer vector being used as the uniques format
  rval = list()
  rval$genotypes <- list()
  for(i in seq(length(uniques))) {
    rval$genotypes[[i]] <- as.integer(clustering[[i]]$abundance)
    names(rval$genotypes[[i]]) <- clustering[[i]]$sequence
  }
  rval$clustering <- clustering
  rval$substitutions <- substitutions
  
  if(length(rval$genotypes)==1) { # one sample, return a naked uniques vector
    rval$genotypes <- rval$genotypes[[1]]
    rval$clustering <- rval$clustering[[1]]
    rval$substitutions <- rval$substitutions[[1]]
  } else { # keep names if it is a list
    names(rval$genotypes) <- names(uniques)
    names(rval$clustering) <- names(uniques)
    names(rval$substitutions) <- names(uniques)
  }

  # Return the error rate(s) used as well as the final sub matrix and estimated error matrix
  if(self_consist) { # Did a self-consist loop
    rval$err <- errs
  } else {
    rval$err <- err
  }
  rval$subs <- trans
  rval$err_out <- trans + 1
  rval$err_out <- t(apply(rval$err_out, 1, function(x) x/sum(x)))
  
  # Store the call and the options that were used in the return object
  rval$opts <- opts
  rval$call <- call
  
  rval
}
