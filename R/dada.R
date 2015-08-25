dada_opts <- new.env()
assign("OMEGA_A", 1e-40, envir = dada_opts)
assign("USE_SINGLETONS", FALSE, envir=dada_opts)
assign("OMEGA_S", 1e-3, envir = dada_opts)
assign("USE_KMERS", TRUE, envir = dada_opts)
assign("KDIST_CUTOFF", 0.42, envir = dada_opts)
assign("MAX_CONSIST", 10, envir = dada_opts)
assign("SCORE_MATRIX", matrix(c(5L, -4L, -4L, -4L, -4L, 5L, -4L, -4L, -4L, -4L, 5L, -4L, -4L, -4L, -4L, 5L),
                              nrow=4, byrow=TRUE), envir = dada_opts)
assign("GAP_PENALTY", -8L, envir = dada_opts)
assign("BAND_SIZE", 16, envir = dada_opts)
assign("MAX_CLUST", 0, envir=dada_opts)
assign("MIN_FOLD", 1, envir=dada_opts)
assign("MIN_HAMMING", 1, envir=dada_opts)
assign("USE_QUALS", TRUE, envir=dada_opts)
assign("QMIN", 0, envir=dada_opts) # NON-FUNCTIONAL
assign("QMAX", 40, envir=dada_opts) # NON-FUNCTIONAL
assign("FINAL_CONSENSUS", FALSE, envir=dada_opts)
assign("VERBOSE", FALSE, envir=dada_opts)
# assign("HOMOPOLYMER_GAPPING", FALSE, envir = dada_opts) # NOT YET IMPLEMENTED

#' High resolution sample inference from amplicon data.
#' 
#' The dada function takes as input dereplicated amplicon sequencing reads and returns the inferred composition
#'  of the sample (or samples). Put another way, dada removes all sequencing errors to reveal the members of the
#'  sequenced community.
#'  
#' If dada is run in selfConsist=TRUE mode, the algorithm will infer both the sample composition and
#'  the parameters of its error model from the data.
#'  
#' @param derep (Required). A derep-class object, the output of \code{\link{derepFastq}}.
#' 
#'  A list of derep objects can be provided, in which case each will be independently denoised with
#'  a shared error model.
#'  
#' @param err (Required). 16xN numeric matrix. Each entry must be between 0 and 1.
#' 
#'  The matrix of estimated rates for each possible nucleotide transition (from sample nucleotide to read nucleotide).
#'  
#'  Rows correspond to the 16 possible transitions (t_ij) indexed as so... 
#'    1:A->A,  2:A->C,  3:A->G,  4:A->T,  5:C->A,  6:C->C,  7:C->G,  8:C->T,
#'    9:G->A, 10:G->C, 11:G->G, 12:G->T, 13:T->A, 14:T->C, 15:T->G, 16:T->T
#'    
#'  Columns correspond to consensus quality scores. Typically there are 41 columns for the quality scores 0-40.
#'  However, if USE_QUALS=FALSE, the matrix must have only one column.
#'    
#' @param errorEstimationFunction (Optional). Function. Default Null.
#' 
#'  If USE_QUALS = TRUE, \code{errorEstimationFunction(dada()$trans_out)} is computed after sample inference finishes 
#'    and the return value is used as the new estimate of the err matrix.
#'    
#'  If USE_QUALS = FALSE, this argument is ignored, and transition rates are estimated by maximum likelihood (t_ij = n_ij/n_i).
#'  
#' @param selfConsist (Optional). \code{logical(1)}. Default FALSE.
#' 
#'  If selfConsist = TRUE, the algorithm will alternate between sample inference and error rate estimation until convergence.
#'    Error rate estimation is performed by the errorEstimationFunction, which is required for selfConsist mode. If dada is
#'    run in selfConsist mode without specifying this function, the default loessErrfun will be used.
#'    
#'  If selfConsist=FALSE the algorithm performs one round of sample inference based on the provided err matrix.
#'   
#' @param ... (Optional). All dada_opts can be passed in as arguments to the dada() function.
#' 
#'  See \code{\link{setDadaOpt}} for a discussion of the various dada options. 
#'
#' @return A \code{\link{dada-class}} object or list of such objects of a list of derep objects was provided. 
#'   
#' @details
#' 
#' dada() implements the Divisive Amplicon Denoising Algorithm as described in:
#' 
#' Briefly, DADA implements a statiscal test for the notion that a specific sequence was seen too many times
#'  to have been caused by amplicon errors from currently inferred sample sequences. Overly-abundant
#'  sequences are used as the seeds of new clusters of sequencing reads, and the final set of clusters
#'  is taken to represent the denoised composition of the sample. A more detailed explanation of the algorithm
#'  is found in two publications:
#' 
#' \itemize{ 
#'  \item{Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J., & Holmes, S. P. (2015). DADA2: High resolution sample inference from amplicon data. bioRxiv, 024034.}
#'  \item{Rosen, M. J., Callahan, B. J., Fisher, D. S., & Holmes, S. P. (2012). Denoising PCR-amplified metagenome data. BMC bioinformatics, 13(1), 283.}
#' }
#'  
#' DADA depends on a parametric error model of substitutions. Thus the quality of its sample inference is affected
#'  by the accuracy of the estimated error rates. DADA's selfConsist mode allows these error rates to be inferred 
#'  from the data.
#'  
#' All of DADA's comparisons between sequences depend on pairwise alignments. This step is the most computationally
#'  intensive part of the algorithm, and two alignment heuristics have been implemented for speed: A kmer-distance
#'  screen and a banded Needleman-Wunsch alignmemt. See \code{\link{setDadaOpt}}.
#'  
#' @seealso 
#'  \code{\link{derepFastq}}
#' 
#'  \code{\link{setDadaOpt}}
#'
#' @export
#'
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' derep2 = derepFastq(system.file("extdata", "sam2F.fastq.gz", package="dada2"))
#' dada(derep1, err=tperr1)
#' dada(list(sam1=derep1, sam2=derep2), err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' dada(derep1, err=inflateErr(tperr1,2), BAND_SIZE=32, OMEGA_A=1e-20)
#'
dada <- function(derep,
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
        warning("Did not provide an error function for selfConsist mode, using the default loessErrfun.")
        errorEstimationFunction <- loessErrfun
      } else {
        message("No error function provided, no post-dada error estimates ($err_out) will be inferred.") 
      }
    } else {
      if(!is.function(errorEstimationFunction)) stop("Must provide a function for errorEstimationFunction.")
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
  rval2 = replicate(length(derep), list(denoised=NULL, clustering=NULL, quality=NULL, subpos=NULL, trans=NULL, map=NULL, uniques_in=NULL,
                                          err_in=NULL, err_out=NULL, opts=NULL, call=NULL), simplify=FALSE)
  for(i in seq_along(derep)) {
    rval2[[i]]$denoised <- getUniques(clustering[[i]])
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

################################################################################
#' Set DADA options
#'
#' setDadaOpt sets the default options used by the dada(...) function for your current session, much
#'  like \code{par} sets the session default plotting parameters. However, all dada options can be set as
#'  part of the dada(...) function call itself by including a DADA_OPTION_NAME=VALUE argument.
#'
#' @param ... (Required). The DADA options to set, along with their new value.
#'  
#' @details The various dada options...
#' 
#' OMEGA_A: This parameter sets the threshold for when DADA2 calls unique sequences significantly overabundant, and therefore creates a
#'  new cluster with that sequence as the center. The default value is 1e-40, which is a conservative setting to avoid making false
#'  positive inferences, but which comes at the cost of reducing the ability to identify some rare variants.
#' 
#' USE_QUALS: If TRUE, the dada(...) error model takes into account the consensus quality score of the dereplicated unique sequences.
#'  If FALSE, quality scores are ignored. The default is TRUE, however if applying DADA2 to pyrosequenced data it is recommended to set
#'  USE_QUALS to FALSE, as quality scores are not informative about substitution error rates in pyrosequencing.
#' 
#' USE_KMERS: If TRUE, a 5-mer distance screen is performed prior to performing each pairwise alignment, and if the 5mer-distance
#'  is greater than KDIST_CUTOFF, no alignment is performed. TRUE by default.
#' 
#' KDIST_CUTOFF: The default value of 0.42 was chosen to screen pairs of sequences that differ by >10\%, and was
#'  calibrated on Illumina sequenced 16S amplicon data. The assumption is that sequences that differ by such a large
#'  amount cannot be linked by amplicon errors (i.e. if you sequence one, you won't get a read of other) and so
#'  careful (and costly) alignment is unnecessary.
#' 
#' BAND_SIZE: When set, banded Needleman-Wunsch alignments are performed. Banding restricts the net cumulative number of insertion
#'  of one sequence relative to the other. The default value of BAND_SIZE is 16. If DADA is applied to marker genes with high rates
#'  of indels, such as the ITS region in fungi, the BAND_SIZE parameter should be increased. Setting BAND_SIZE to a negative number
#'  turns off banding (i.e. full Needleman-Wunsch).
#' 
#' SCORE_MATRIX: The score matrix for the Needleman-Wunsch alignment. This is a 4x4 matrix as no ambiguous nucleotides
#'  are allowed. Default is nuc44: -4 for mismatches, +5 for matchces.
#'  
#' GAP_PENALTY: The cost of gaps in the Needlman-Wunsch alignment. Default is -8.
#'  
#' MIN_FOLD: The minimum fold-overabundance for sequences to form new clusters. Default value is 1, which means this
#'  criteria is ignored.
#'  
#' MIN_HAMMING: The minimum hamming-separation for sequences to form new clusters. Default value is 1, which means this
#'  criteria is ignored.
#'
#' MAX_CLUST: The maximum number of clusters. Once this many clusters have been created, the algorithm terminates regardless
#'  of whether the statistical model suggests more sample sequences exist. If set to 0 this argument is ignored. Default
#'  value is 0.
#'  
#' VERBOSE: If TRUE progress messages from the algorithm are printed. Warning: There is a lot of output. Default is FALSE.
#' 
#' USE_SINGLETONS: CURRENTLY BROKEN. Default is FALSE.
#' 
#' OMEGA_S: CURRENTLY BROKEN. Default is 1e-3.
#'
#' @seealso 
#'  \code{\link{getDadaOpt}}
#'
#' @export
#'
#' @examples
#' setDadaOpt(OMEGA_A = 1e-20)
#' setDadaOpt(OMEGA_A = 1e-20, VERBOSE = TRUE)
#'   
setDadaOpt <- function(...) {
  opts <- getDadaOpt()
  args <- list(...)
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) {
      if(class(getDadaOpt(opnm)) != class(args[[opnm]])) {
        warning(paste0(opnm, " not set, value provided has different class (", class(args[[opnm]]), 
                    ") then current option value (", class(getDadaOpt(opnm)), ")"))
      } else {
        assign(opnm, args[[opnm]], envir=dada_opts)
      }
    } else {
      warning(opnm, " is not a valid DADA option.")
    }
  }
}

# Should add in more sanity checking here??
# matrix dimensions, general structure of score matrix
#  if(!(is.numeric(score) && dim(err) == c(4,4))) {
#    stop("dada: Invalid score matrix.")
#  }
#  
#  if(!(is.numeric(gap_penalty) && gap_penalty <=0)) {
#    stop("dada: Invalid gap penalty.")
#  }
#  if(gap_penalty > -1) warning("dada: Very small gap penalty.")

################################################################################
#' Get DADA options
#'
#' @param option (Optional). Character.
#'  The DADA option(s) to get.
#' 
#' @return Named list of option/value pairs.
#'  Returns NULL if an invalid option is requested.
#' 
#' @seealso 
#'  \code{\link{setDadaOpt}}
#'
#' @export
#'
#' @examples
#' getDadaOpt("BAND_SIZE")
#' getDadaOpt()
#' 
getDadaOpt <- function(option = NULL) {
  if(is.null(option)) option <- ls(dada_opts)
  
  if(!all(option %in% ls(dada_opts))) {
    warning("Tried to get an invalid DADA option: ", option[!(option %in% ls(dada_opts))])
    option <- option[option %in% ls(dada_opts)]
  }
  
  ropts <- lapply(option, function(x) get(x, envir=dada_opts))
  names(ropts) <- option
  if(length(ropts) == 1) ropts <- ropts[[1]]  # If just one option requested, return it alone
  return(ropts)
}
