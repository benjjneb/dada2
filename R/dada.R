#' dada infers the sample sequence in amplicon data.
#' 
#' The dada function takes as input the unique sequences in an amplicon sequencing sample paired
#'  with their abundances, and returns the sample genotypes paired with their abundances
#'  inferred by the Divisive Amplicon Denoising Algorithm:
#'  (Rosen, Callahan, Fisher, Holmes 2012; Callahan, McMurdie, Rosen, Han, Johnson, Holmes 2015).
#'
#' @param derep (Required). A derep-class object, the output of derepFastq.
#'  The derep object contains the $uniques named integer vector of abundances, named by the unique
#'  DNA sequence to which that abundance corresponds. It also contains $quals, a numeric matrix of 
#'  the associated consensus quality scores (by averaging) with a row for each corresponding unique
#'  in the uniques vector, and a column for each sequence position. If USE_QUALS = TRUE then the 
#'  dada() error model uses the associated consensus quality scores.
#'  
#'  A list of derep objects can be provided, in which case each will be independently denoised but
#'  the error model will be shared across these samples when denoising.
#'  
#'  The dada(...) function requires that all input sequences have been trimmed to the same length.
#'  
#' @param err (Required). 16xN numeric matrix. Each entry must be between 0 and 1.
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
#' @param errorEstimationFunction (Optional). Function. Default Null.
#'   If USE_QUALS = TRUE, errorEstimationFunction(dada()$trans_out) is computed and taken to be the new err matrix.
#'    If selfConsist = TRUE, the next iteration of the dada algorithm will use these new error rates, and this
#'    will continue until the loop terminates due to convergence (or by hitting MAX_CONSIST). If
#'    selfConsist=FALSE, errorEstimationFunction is only used to calculate return value $err_out, but otherwise
#'    has no effect on the algorithm.
#'   
#'   If USE_QUALS = FALSE, this argument is ignored, and transition rates are estimated by maximum likelihood (t_ij = n_ij/n_i).
#'  
#' @param selfConsist (Optional). \code{logical(1)}. Default FALSE.
#'  When true, DADA will re-estimate error rates after inferring sample genotypes, and then repeat
#'  the algorithm using the newly estimated error rates. This continues until convergence.
#'
#' @param ... (Optional). All dada_opts can be passed in as arguments to the dada() function.
#'    eg. dada(unq, err=err_in, USE_QUALS=TRUE, OMEGA_A=1e-50, MAX_CLUST=50). See \code{\link{setDadaOpt}}
#'    for a discussion of the various dada options. 
#'
#' @return A \code{\link{dada-class}} object. 
#'   
#'  Not that if \code{selfConsist=TRUE}, 
#'  then \code{$err_in} is a list of length the number of times 
#'  through the selfConsist loop corresponding to the err used in each iteration. 
#'  \code{$err_out} is the final estimated error rate.
#'   
#'  If a list of uniques vectors was provided (i.e. multiple samples)
#'  then a list of the \code{\link{dada-class}} objects is returned corresponding
#'  to each input sample.
#'  
#' @seealso \code{\link{derepFastq}}, \code{\link{setDadaOpt}}
#'
#' @export
#'
dada <- function(derep, #!!!!!
                 err,
                 errorEstimationFunction = NULL,
                 selfConsist = FALSE, ...) {
  
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
  
  # If a single derep object, make into a length 1 list
  if(class(derep) == "derep") { derep <- list(derep) }
  if(opts$USE_QUALS && any(is.null(lapply(derep, function(x) x$quals)))) { stop("The input derep object(s) must include quals if USE_QUALS is TRUE.") }
  
  # Validate derep object(s)
  for(i in seq_along(derep)) {
    if(!class(derep[[i]]) == "derep") stop("The derep argument must be a derep-class object or list of derep-class objects.")
    if(!(is.integer(derep[[i]]$uniques))) {
      stop("Invalid derep$uniques vector. Must be integer valued.")
    }
    if(!(all(sapply(names(derep[[i]]$uniques), function(x) nchar(gsub("[ACGT]", "", x))==0, USE.NAMES=FALSE)))) {
      stop("Invalid derep$uniques vector. Names must be sequences made up only of A/C/G/T.")
    }
    if(sum(tabulate(nchar(names(derep[[i]]$uniques)))>0) > 1) {
      stop("Invalid derep$uniques vector. All sequences must be the same length.")
    }
  }

  # Validate quals matrix(es)
  if(opts$USE_QUALS) {
    for(i in seq_along(derep)) {
      if(nrow(derep[[i]]$quals) != length(derep[[i]]$uniques)) {
        stop("derep$qual matrices must have one row for each derep$unique sequence.")
      }
      if(any(sapply(names(derep[[i]]$uniques), nchar) > ncol(derep[[i]]$quals))) {
        stop("derep$qual matrices must have as many columns as the length of the derep$unique sequences.")
      }
      if(any(is.na(derep[[i]]$quals))) {
        stop("NAs in derep$qual matrix. Check that all input sequences were the same length.")
      }
      if(min(derep[[i]]$quals) < opts$QMIN || max(derep[[i]]$quals > opts$QMAX)) {
        stop("Invalid derep$qual matrix. Quality values must be between QMIN and QMAX.")
      }
    }
  }
  
  # Validate err matrix
  if(!is.numeric(err)) stop("Error matrix must be numeric.")
  if(!(nrow(err)==16)) stop("Error matrix must have 16 rows (A2A, A2C, ...).")
  if(!all(err>=0)) stop("All error matrix entries must be >= 0.")
  if(!all(err<=1)) stop("All error matrix entries must be <=1.")
  if(any(err==0)) warning("Zero in error matrix.")
  # Might want to check for summed transitions from NT < 1 also.
  
  # Validate err_model
  if(!opts$USE_QUALS) {
    if(!is.null(errorEstimationFunction)) warning("The errorEstimationFunction argument is ignored when USE_QUALS is FALSE.")
    errorEstimationFunction = NULL  # NULL error function has different meaning depending on USE_QUALS
  } else {
    if(is.null(errorEstimationFunction)) { 
      if(selfConsist) {
        stop("Must provide an error function if USE_QUALS and selfConsist=TRUE.")
      } else {
        message("No error function provided, no post-dada error estimates ($err_out) will be inferred.") 
      }
    } else {
      if(!is.function(errorEstimationFunction)) stop("Must provide a function for errorEstimationFunction.")
#      if(any(names(formals(err_model)) != c("parms", "qave"))) stop("err_model must be a function of two named arguments: err_model(parms, qave).")
    }
  }
  
  # Initialize
  cur <- NULL
  nconsist <- 1
  errs <- list()
  # The main loop, run once, or repeat until err repeats if selfConsist=T

  repeat{
    clustering <- list()
    clusterquals <- list()
    subpos <- list()
    trans <- list()
    map <- list()
    exp <- list()
    prev <- cur
    errs[[nconsist]] <- err

    for(i in seq_along(derep)) {
      if(!opts$USE_QUALS) { qi <- matrix(0, nrow=0, ncol=0) }
      else { qi <- unname(t(derep[[i]]$quals)) } # Need transpose so that sequences are columns

      if(nconsist == 1) {
        cat("Sample", i, "-", sum(derep[[i]]$uniques), "reads in", length(derep[[i]]$uniques), "unique sequences.\n")
      } else if(i==1) {
        cat("   Consist step", nconsist, "\n")
      }
      res <- dada_uniques(names(derep[[i]]$uniques), unname(derep[[i]]$uniques), err, qi, 
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
      res$clustering$sequence <- as.character(res$clustering$sequence)
      # ... nothing here for now
      
      # List the returns
      clustering[[i]] <- res$clustering
      clusterquals[[i]] <- t(res$clusterquals) # make sequences rows and positions columns
      subpos[[i]] <- res$subpos
      trans[[i]] <- res$subqual
      map[[i]] <- res$map
      exp[[i]] <- res$exp
      rownames(trans[[i]]) <- c("A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T", "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T")
      if(opts$USE_QUALS) colnames(trans[[i]]) <- seq(opts$QMIN, opts$QMAX)  # Assumes C sides is returning one col for each integer from QMIN to QMAX
    }
    # Accumulate the sub matrix
    cur <- Reduce("+", trans) # The only thing that changes is err(trans), so this is sufficient
    
    # Estimate the new error model (if applicable)
    if(opts$USE_QUALS) {
      if(is.null(errorEstimationFunction)) {
        err <- NULL
      } else {
        err <- errorEstimationFunction(cur)
      }
    } else { # Not using quals, MLE estimate for each transition type
      err <- cur + 1   # ADD ONE PSEUDOCOUNT TO EACH TRANSITION
      err[1:4,1] <- err[1:4,1]/sum(err[1:4,1])
      err[5:8,1] <- err[5:8,1]/sum(err[5:8,1])
      err[9:12,1] <- err[9:12,1]/sum(err[9:12,1])
      err[13:16,1] <- err[13:16,1]/sum(err[13:16,1])
    }

    if(selfConsist) { # Validate err matrix
      if(!is.numeric(err)) stop("Error matrix returned by errorEstimationFunction not numeric.")
      if(!(nrow(err)==16)) stop("Error matrix returned by errorEstimationFunction does not have 16 rows.")
      if(!all(err>=0)) stop("Error matrix returned by errorEstimationFunction has entries <0.")
      if(!all(err<=1)) stop("Error matrix returned by errorEstimationFunction has entries >1.")
      if(any(err==0)) warning("Error matrix returned by errorEstimationFunction has 0 entries.")      
    }
    
    # Termination condition for selfConsist loop
    if((!selfConsist) || identical(cur, prev) || (nconsist >= opts$MAX_CONSIST)) {
      break
    } 
    nconsist <- nconsist+1
  } # repeat

  cat("\n")
  if(selfConsist) {
    if(nconsist >= opts$MAX_CONSIST) {
      warning("dada: Self-consistency loop terminated before convergence.")
    } else {
      cat("\nConvergence after ", nconsist, " rounds.\n")
    }
  }
  
  # Construct return object
  # A single dada-class object if one derep object provided.
  # A list of dada-class objects if multiple derep objects provided.
  rval2 = replicate(length(derep), list(genotypes=NULL, clustering=NULL, quality=NULL, subpos=NULL, trans=NULL, map=NULL, uniques_in=NULL,
                                          err_in=NULL, err_out=NULL, opts=NULL, call=NULL), simplify=FALSE)
  for(i in seq_along(derep)) {
    # Convert cur$genotypes to the named integer vector being used as the uniques format
    rval2[[i]]$genotypes <- as.uniques(clustering[[i]])
##!    names(rval2[[i]]$genotypes) <- clustering[[i]]$sequence
    rval2[[i]]$clustering <- clustering[[i]]
    rval2[[i]]$quality <- clusterquals[[i]]
    rval2[[i]]$subpos <- subpos[[i]]
    rval2[[i]]$trans <- trans[[i]]
    rval2[[i]]$map <- map[[i]]
    rval2[[i]]$uniques_in <- derep[[i]]$uniques
    rval2[[i]]$exp <- exp[[i]]
    # Return the error rate(s) used as well as the final estimated error matrix
    if(selfConsist) { # Did a self-consist loop
      rval2[[i]]$err_in <- errs
    } else {
      rval2[[i]]$err_in <- errs[[1]]
    }
    rval2[[i]]$err_out <- err           # maybe better as _final? Just the last one
    
    # Store the call and the options that were used in the return object
    rval2[[i]]$opts <- opts
    rval2[[i]]$call <- call
  }
  names(rval2) <- names(derep)
  if(length(rval2) == 1) {  # Unlist if just a single derep object provided
    rval2 <- rval2[[1]]
    rval2 <- as(rval2, "dada")
  } else {
    for(i in seq_along(rval2)) {
      rval2[[i]] <- as(rval2[[i]], "dada")
    }
  }

  return(rval2)
}
