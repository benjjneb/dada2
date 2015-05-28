#' dada infers the sample sequence in amplicon data.
#' 
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
#' @param err (Required). 16xN numeric matrix.
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
#' @param err_function (Optional). Function.
#'   If USE_QUALS = TRUE, err_function(dada()$trans_out) is computed and taken to be the new err matrix.
#'    If self_consist = TRUE, the next iteration if the dada() algorithm will use this new err, and this
#'    will continue until the loop terminates due to convergence (or by hitting MAX_CONSIST). If
#'    self_consist=FALSE, err_function is only used to calculate return value $err_out, but otherwise
#'    has no effect on the algorithm.
#'   
#'   If USE_QUALS = FALSE, this argument is ignored, and transition rates are estimated by maximum likelihood (t_ij = n_ij/n_i).
#'  
#' @param self_consist (Optional). \code{logical(1)}
#'  When true, DADA will re-estimate error rates after inferring sample genotypes, and then repeat
#'  the algorithm using the newly estimated error rates. This continues until convergence.
#'
#' @param ... (Optional). All dada_opts can be passed in as arguments to the dada() function.
#'    eg. dada(unq, err=err_in, OMEGA_A=1e-50, MAX_CLUST=50) 
#'
#' @return A multi-item List with the following named values...
#' \itemize{
#'  \item{$genotypes: }{Integer vector, named by sequence valued by abundance, of the denoised genotypes.}
#'  \item{$clustering: }{An informative data.frame containing information on each cluster.}
#'  \item{$quality: }{The average quality scores for each cluster (row) by position (col).}
#'  \item{$map: }{Integer vector that maps the unique (index) to the cluster/genotype (value).}
#'  \item{$birth_subs: }{A data.frame containing the substitutions at the birth of each new cluster.}
#'  \item{$trans: }{The matrix of transitions by type (row), eg. A2A, A2C..., and quality score (col)
#'          observed in the final output of the dada algorithm.}
#'  \item{$trans_out: }{The matrix of transitions by type and quality, summed over all denoised samples.
#'          Used to estimate $err_out. Identical to $trans if just one sample denoised.}
#'  \item{$err_out: }{The err matrix estimated from the output of dada. NULL if err_function not provided.}
#'  \item{$err_in: }{The err matrix used for this invocation of dada.}
#'  \item{$opts: }{A list of the dada_opts used for this invocation of dada.}
#'  \item{$uniques: }{The uniques vector(s) used for this invocation of dada.}
#'  \item{$call: }{The function call used for this invocation of dada.}
#' }
#'   
#'  If a list of uniques vectors was provided (i.e. multiple samples) then $genotypes, $clustering,
#'    $quality, $map and $birth_subs are lists with entries corresponding to each provided uniques vector.
#'  
#'  If self_consist=TRUE, $err_in is a list of length the number of times through the self_consist
#'    loop corresponding to the err used in each iteration. $err_out is the final estimated error rate.
#'   
#' @export
#'
dada <- function(uniques, quals=NULL,
                 err,
                 err_function = NULL,
                 self_consist = FALSE, ...) {
  
  call <- sys.call(1)
  # Read in default opts and then replace with any that were passed in to the function
  opts <- getDadaOpt()
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
  
  # Validate err_model
  if(!opts$USE_QUALS) {
    if(!is.null(err_function)) warning("The err_function argument is ignored when USE_QUALS is FALSE.")
    err_function = NULL  # NULL error function has different meaning depending on USE_QUALS
  } else {
    if(is.null(err_function)) { 
      if(self_consist) {
        stop("Must provide an error function if USE_QUALS and self_consist=TRUE.")
      } else {
        message("No error function provided, no post-dada error estimates ($err_out) will be inferred.") 
      }
    } else {
      if(!is.function(err_function)) stop("Must provide a function for err_function.")
#      if(any(names(formals(err_model)) != c("parms", "qave"))) stop("err_model must be a function of two named arguments: err_model(parms, qave).")
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
    trans <- list()
    map <- list()
    prev <- cur
    errs[[nconsist]] <- err

    for(i in seq(length(uniques))) {
      if(is.null(quals)) { qi <- matrix(0, nrow=0, ncol=0) }
      else { qi <- unname(t(quals[[i]])) } # Need transpose so that sequences are columns
      if(nconsist == 1) {
        cat("Sample", i, "-", sum(uniques[[i]]), "reads in", length(uniques[[i]]), "unique sequences.\n")
      } else if(i==1) {
        cat("   Consist step", nconsist, "\n")
      }
      res <- dada_uniques(names(uniques[[i]]), unname(uniques[[i]]), err, qi, 
                          opts[["SCORE_MATRIX"]], opts[["GAP_PENALTY"]],
                          opts[["USE_KMERS"]], opts[["KDIST_CUTOFF"]],
                          opts[["BAND_SIZE"]],
                          opts[["OMEGA_A"]], 
                          opts[["USE_SINGLETONS"]], opts[["OMEGA_S"]],
                          opts[["MAX_CLUST"]],
                          opts[["MIN_FOLD"]], opts[["MIN_HAMMING"]],
                          opts[["USE_QUALS"]],
                          opts[["QMIN"]], opts[["QMAX"]],
                          opts[["FINAL_CONSENSUS"]],
                          opts[["VERBOSE"]])
      
      # Augment the returns
      # res$clustering$ham_nln <- 
      
      # List the returns
      clustering[[i]] <- res$clustering
      clusterquals[[i]] <- t(res$clusterquals) # make sequences rows and positions columns
      subpos[[i]] <- res$subpos
      trans[[i]] <- res$subqual
      map[[i]] <- res$map
      rownames(trans[[i]]) <- c("A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T", "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T")
      if(!is.null(quals)) colnames(trans[[i]]) <- seq(opts$QMIN, opts$QMAX)  # Assumes C sides is returning one col for each integer from QMIN to QMAX
    }
    # Accumulate the sub matrix
    cur <- Reduce("+", trans) # The only thing that changes is err(trans), so this is sufficient
    
    # Estimate the new error model (if applicable)
    if(opts$USE_QUALS) {
      if(is.null(err_function)) {
        err <- NULL
      } else {
        err <- err_function(cur)
      }
    } else { # Not using quals, MLE estimate for each transition type
      err <- cur + 1   # ADD ONE PSEUDOCOUNT TO EACH TRANSITION
      err[1:4,1] <- err[1:4,1]/sum(err[1:4,1])
      err[5:8,1] <- err[5:8,1]/sum(err[5:8,1])
      err[9:12,1] <- err[9:12,1]/sum(err[9:12,1])
      err[13:16,1] <- err[13:16,1]/sum(err[13:16,1])
    }

    # Termination condition for self_consist loop
    if((!self_consist) || identical(cur, prev) || (nconsist >= opts$MAX_CONSIST)) {
      break
    } 
    nconsist <- nconsist+1
  } # repeat

  cat("\n")
  if(self_consist) {
    if(nconsist >= opts$MAX_CONSIST) {
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
  rval$quality <- clusterquals
  rval$subpos <- subpos
  rval$trans <- trans
  rval$map <- map
  rval$uniques <- uniques

  if(length(rval$genotypes)==1) { # one sample, return a naked uniques vector
    rval$genotypes <- rval$genotypes[[1]]
    rval$clustering <- rval$clustering[[1]]
    rval$quality <- rval$quality[[1]]
    rval$subpos <- rval$subpos[[1]]
    rval$trans <- rval$trans[[1]]
    rval$map <- rval$map[[1]]
    rval$uniques <- rval$uniques[[1]]
  } else { # keep names if it is a list
    names(rval$genotypes) <- names(uniques)
    names(rval$clustering) <- names(uniques)
    names(rval$quality) <- names(uniques)
    names(rval$subpos) <- names(uniques)
    names(rval$trans) <- names(uniques)
    names(rval$map) <- names(uniques)
  }

  # Return the error rate(s) used as well as the final sub matrix and estimated error matrix
  if(self_consist) { # Did a self-consist loop
    rval$err_in <- errs
  } else {
    rval$err_in <- errs[[1]]
  }
  rval$trans_out <- cur
  rval$err_out <- err
  
  # Store the call and the options that were used in the return object
  rval$opts <- opts
  rval$call <- call
  
  rval
}
