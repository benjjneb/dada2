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
#' @param quals (Optional). Numeric matrix containing average quality scores for each unique at each position.
#'  The quals matrix has a row for each corresponding unique in the uniques vector, and a column for each
#'  sequence position. So nrow(quals) == length(uniques) and ncol(quals) == nchar(uniques[[foo]]).
#'  
#' @param err (Optional). 16xN numeric matrix.
#'  The matrix of estimated rates for each possible nucleotide transition (from sample nucleotide to read nucleotide).
#'  Rows correspond to the 16 possible transitions (t_ij) indexed as so... 
#'    1:A->A,  2:A->C,  3:A->G,  4:A->T,  5:C->A,  6:C->C,  7:C->G,  8:C->T,
#'    9:G->A, 10:G->C, 11:G->G, 12:G->T, 13:T->A, 14:T->C, 15:T->G, 16:T->T
#'    
#'  If USE_QUALS = TRUE, the columns correspond to a linear interpolation of quality score from QMIN and QMAX, eg. seq(QMIN, QMAX, by=QSTEP).
#'    err[t_ij, round((q_ave-QMIN)/QSTEP) + 1] = Prob(j in sequence | i in sample genotype and q_ave).
#'    
#'  If USE_QUALS = FALSE, the matrix must have only one column, which corresponds to the estimated error rate for that transition.
#'    err[t_ij, 1] = Prob(j in sequence | i in sample genotype).
#'  
#' @param err_model (Optional). A parameterized function that predicts error rates from a given transition type and quality score.
#'   If USE_QUALS = TRUE, 
#'   
#'   If USE_QUALS = FALSE, this argument is ignored, and transition rates are estimated by maximum likelihood (t_ij = n_ij/n_i).
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
dada <- function(uniques, quals=NULL,
                 err = matrix(rep(c(0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991, 0.003, 0.003, 0.003, 0.003, 0.991), each=41), nrow=16, byrow=T),
                 err_model = NULL,
                 self_consist = FALSE, ...) {
  
  call <- sys.call(1)
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
  
  # If a single vector, make into a length 1 list
  if(!is.list(uniques)) { uniques <- list(uniques) }
  if(opts$USE_QUALS && is.null(quals)) { stop("Must provide quals if USE_QUALS is TRUE.") }
  if(!is.null(quals) && !is.list(quals)) { quals <- list(quals) }
  
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

  # Validate quals matrix(es)
  if(!is.null(quals)) {
    if(length(uniques) != length(quals)) { stop("Must be a qual matrix for each uniques vector.") }
    for(i in seq(length(uniques))) {
      if(nrow(quals[[i]]) != length(uniques[[i]])) {
        stop("Qual matrices must have one row for each unique.")
      }
      if(any(sapply(names(uniques), nchar) > ncol(quals[[i]]))) {
        stop("Qual matrices must have at least as many columns as the sequence length.")
      }
      if(min(quals[[i]]) < opts$QMIN || max(quals[[i]] > opts$QMAX)) {
        stop("Invalid quality matrix. Quality values must be between QMIN and QMAX.")
      }
    }
  }
  
  # Validate err matrix
  if(!( is.numeric(err) && nrow(err) == 16) && all(err>=0))
  { stop("Invalid error matrix.") }
  if(any(err==0)) warning("Zero in error matrix.")
  # MAY WANT TO ADD BACK IN CHECKING FOR SUMS TO 1? Nah, it may be useful to have "super probabilistic" error matrices
  
  # Validate err_model
  if(!opts$USE_QUALS) {
    if(!is.null(err_model)) warning("The err_model argument is ignored when USE_QUALS is FALSE.")
    err_model = NULL  # NULL error model has different meaning depending on USE_QUALS
  } else {
    if(is.null(err_model)) { 
      if(self_consist) {
        stop("Must provide an error model if USE_QUALS and self_consist=TRUE.")
      } else {
        message("No error model provided, no error estimates will be inferred.") 
      }
    } else {
      if(!is.function(err_model)) stop("Must provide a function for err_model if USE_QUALS is TRUE.")
      if(any(names(formals(err_model)) != c("parms", "qave"))) stop("err_model must be a function of two named arguments: err_model(parms, qave).")
    }
  }
  
  # Initialize
  cur <- NULL
  nconsist <- 1
  errs <- list()
  # The main loop, run once, or repeat until err repeats if self_consist=T

  repeat{
    clustering <- list()
    clusterquals <- list()
    subpos <- list()
    subqual <- list()
    map <- list()
    prev <- cur
    errs[[nconsist]] <- err

    for(i in seq(length(uniques))) {
      if(is.null(quals)) { qi <- matrix(0, nrow=0, ncol=0) }
      else { qi <- unname(t(quals[[i]])) } # Need transpose so that sequences are columns
      cat("Sample", i, "-", sum(uniques[[i]]), "reads in", length(uniques[[i]]), "unique sequences.\n")
      res <- dada_uniques(names(uniques[[i]]), unname(uniques[[i]]), err, qi, 
                          opts[["SCORE_MATRIX"]], opts[["GAP_PENALTY"]],
                          opts[["USE_KMERS"]], opts[["KDIST_CUTOFF"]],
                          opts[["BAND_SIZE"]],
                          opts[["OMEGA_A"]], 
                          opts[["USE_SINGLETONS"]], opts[["OMEGA_S"]],
                          opts[["MAX_CLUST"]],
                          opts[["MIN_FOLD"]], opts[["MIN_HAMMING"]],
                          opts[["USE_QUALS"]],
                          opts[["QMIN"]], opts[["QMAX"]])
      
      # Augment the returns
      # res$clustering$ham <- sapply(res$clustering$sequence, function(x) nrow(strdiff(res$clustering$sequence[[1]], x)))
      
      # List the returns
      clustering[[i]] <- res$clustering
      clusterquals[[i]] <- t(res$clusterquals) # make sequences rows and positions columns
      subpos[[i]] <- res$subpos
      subqual[[i]] <- res$subqual
      map[[i]] <- res$map
      rownames(subqual[[i]]) <- c("A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T", "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T")
      if(!is.null(quals)) colnames(subqual[[i]]) <- seq(opts$QMIN, opts$QMAX)
    }
    # Accumulate the sub matrix
    trans <- Reduce("+", subqual)
    cur = trans # The only thing that changes is err which is set by trans, so this is sufficient
    
    # Estimate the new error model (if applicable)
    if(opts$USE_QUALS) {
      if(is.null(err_model)) {
        err <- NULL
      } else { 
        qq <- as.numeric(colnames(trans))
        est <- matrix(0, nrow=0, ncol=length(qq))
        for(nti in c("A","C","G","T")) {
          for(ntj in c("A","C","G","T")) {
            if(nti != ntj) {
              x <- qq
              numer <- trans[paste0(nti,"2",ntj),]
              denom <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
              # Determine which entries to keep, no zeros or NAs or low-outlier denoms
              outs <- 10^(boxplot.stats(log10(denom[denom>0]))$out)
              outs <- outs[outs<median(denom)] # low outliers
              keep <- !(is.na(numer) | is.na(denom) | denom==0 | (denom %in% outs))
              x <- x[keep]
              y <- (numer[keep]+1)/(4+denom[keep]) # Pseudocounts... REVISIT
              mod.tp <- optim(c(-1, 20, 35, -4), make_log10_lsq_obj(err_model,x,y))
              ##### NEED TO DO SOMETHING ABOUT THE INIT PARAMS HERE!!
              pred <- err_model(mod.tp$par, qq)
              est <- rbind(est, pred)
            } # if(nti != ntj)
          } # for(ntj in c("A","C","G","T"))
        } # for(nti in c("A","C","G","T"))
        
        # Expand the err matrix with the self-transition probs
        err <- rbind(1-colSums(est[1:3,]), est[1:3,],
                     est[4,], 1-colSums(est[4:6,]), est[5:6,],
                     est[7:8,], 1-colSums(est[7:9,]), est[9,],
                     est[10:12,], 1-colSums(est[10:12,]))
        rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
        colnames(err) <- seq(opts$QMIN, opts$QMAX)
      }
    } else { # Not using quals
      err <- trans + 1   # ADD ONE PSEUDOCOUNT TO EACH TRANSITION
      err[1:4,1] <- err[1:4,1]/sum(err[1:4,1])
      err[5:8,1] <- err[5:8,1]/sum(err[5:8,1])
      err[9:12,1] <- err[9:12,1]/sum(err[9:12,1])
      err[13:16,1] <- err[13:16,1]/sum(err[13:16,1])
    }

    # Termination condition for self_consist loop
    if((!self_consist) || identical(cur, prev) || (nconsist >= get("MAX_CONSIST", envir=dada_opts))) {
      break
    } 
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
  rval$clusterquals <- clusterquals
  rval$subpos <- subpos
  rval$subqual <- subqual
  rval$map <- map
  
  if(length(rval$genotypes)==1) { # one sample, return a naked uniques vector
    rval$genotypes <- rval$genotypes[[1]]
    rval$clustering <- rval$clustering[[1]]
    rval$clusterquals <- rval$clusterquals[[1]]
    rval$subpos <- rval$subpos[[1]]
    rval$subqual <- rval$subqual[[1]]
    rval$map <- rval$map[[1]]
  } else { # keep names if it is a list
    names(rval$genotypes) <- names(uniques)
    names(rval$clustering) <- names(uniques)
    names(rval$clusterquals) <- names(uniques)
    names(rval$subpos) <- names(uniques)
    names(rval$subqual) <- names(uniques)
    names(rval$map) <- names(uniques)
  }

  # Return the error rate(s) used as well as the final sub matrix and estimated error matrix
  if(self_consist) { # Did a self-consist loop
    rval$err_in <- errs
  } else {
    rval$err_in <- errs[[1]]
  }
  rval$subs <- trans
  rval$err_out <- err
  
  # Store the call and the options that were used in the return object
  rval$opts <- opts
  rval$call <- call
  
  rval
}
