dada_opts <- new.env()
assign("OMEGA_A", 1e-40, envir = dada_opts)
assign("USE_KMERS", TRUE, envir = dada_opts)
assign("KDIST_CUTOFF", 0.42, envir = dada_opts)
assign("MAX_CONSIST", 10, envir = dada_opts)
assign("SCORE_MATRIX", matrix(c(5L, -4L, -4L, -4L, -4L, 5L, -4L, -4L, -4L, -4L, 5L, -4L, -4L, -4L, -4L, 5L),
                              nrow=4, byrow=TRUE), envir = dada_opts)
assign("GAP_PENALTY", -8L, envir = dada_opts)
assign("BAND_SIZE", 16, envir = dada_opts)
assign("VECTORIZED_ALIGNMENT", TRUE, envir = dada_opts)
assign("MAX_CLUST", 0, envir=dada_opts)
assign("MIN_FOLD", 1, envir=dada_opts)
assign("MIN_HAMMING", 1, envir=dada_opts)
assign("USE_QUALS", TRUE, envir=dada_opts)
assign("VERBOSE", FALSE, envir=dada_opts)
assign("HOMOPOLYMER_GAP_PENALTY", NULL, envir = dada_opts)
# assign("FINAL_CONSENSUS", FALSE, envir=dada_opts) # NON-FUNCTIONAL AT THE MOMENT

#' High resolution sample inference from amplicon data.
#' 
#' The dada function takes as input dereplicated amplicon sequencing reads and returns the inferred composition
#'  of the sample (or samples). Put another way, dada removes all sequencing errors to reveal the members of the
#'  sequenced community.
#'  
#' If dada is run in selfConsist=TRUE mode, the algorithm will infer both the sample composition and
#'  the parameters of its error model from the data.
#'  
#' @param derep (Required). A \code{\link{derep-class}} object, the output of \code{\link{derepFastq}}.
#'  A list of such objects can be provided, in which case each will be denoised with a shared error model.
#'  
#' @param err (Required). 16xN numeric matrix. Each entry must be between 0 and 1.
#' 
#'  The matrix of estimated rates for each possible nucleotide transition (from sample nucleotide to read nucleotide).
#'  
#'  Rows correspond to the 16 possible transitions (t_ij) indexed such that 1:A->A, 2:A->C, ..., 16:T->T
#'    
#'  Columns correspond to quality scores. Typically there are 41 columns for the quality scores 0-40.
#'  However, if USE_QUALS=FALSE, the matrix must have only one column.
#'  
#'  If selfConsist = TRUE, \code{err} can be set to NULL and an initial error matrix will be estimated from the data
#'  by assuming that all reads are errors away from one true sequence.
#'    
#' @param errorEstimationFunction (Optional). Function. Default \code{\link{loessErrfun}}.
#' 
#'  If USE_QUALS = TRUE, \code{errorEstimationFunction(dada()$trans_out)} is computed after sample inference, 
#'    and the return value is used as the new estimate of the err matrix in $err_out.
#'    
#'  If USE_QUALS = FALSE, this argument is ignored, and transition rates are estimated by maximum likelihood (t_ij = n_ij/n_i).
#'  
#' @param selfConsist (Optional). \code{logical(1)}. Default FALSE.
#' 
#'  If selfConsist = TRUE, the algorithm will alternate between sample inference and error rate estimation 
#'    until convergence. Error rate estimation is performed by \code{errorEstimationFunction}.
#'    
#'  If selfConsist=FALSE the algorithm performs one round of sample inference based on the provided \code{err} matrix.
#'   
#' @param pool (Optional). \code{logical(1)}. Default is FALSE.
#' 
#'  If pool = TRUE, the algorithm will pool together all samples prior to sample inference.
#'  If pool = FALSE, sample inference is performed on each sample individually.
#'  
#'  This argument has no effect if only 1 sample is provided, and \code{pool} does not affect
#'   error rates, which are always estimated from pooled observations across samples.
#'   
#' @param nthreads (Optional). Default is 1.
#'  The number of threads to use in the core algorithm.   
#'   
#' @param ... (Optional). All dada_opts can be passed in as arguments to the dada() function.
#'  See \code{\link{setDadaOpt}} for a full list and description of these options. 
#'
#' @return A \code{\link{dada-class}} object or list of such objects if a list of dereps was provided. 
#'   
#' @details
#' 
#' Briefly, \code{dada} implements a statiscal test for the notion that a specific sequence was seen too many times
#'  to have been caused by amplicon errors from currently inferred sample sequences. Overly-abundant
#'  sequences are used as the seeds of new clusters of sequencing reads, and the final set of clusters
#'  is taken to represent the denoised composition of the sample. A more detailed explanation of the algorithm
#'  is found in two publications:
#' 
#' \itemize{ 
#'  \item{Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP (2015). DADA2: High resolution sample inference from amplicon data. bioRxiv, 024034.}
#'  \item{Rosen MJ, Callahan BJ, Fisher DS, Holmes SP (2012). Denoising PCR-amplified metagenome data. BMC bioinformatics, 13(1), 283.}
#' }
#'  
#' \code{dada} depends on a parametric error model of substitutions. Thus the quality of its sample inference is affected
#'  by the accuracy of the estimated error rates. \code{selfConsist} mode allows these error rates to be inferred 
#'  from the data.
#'  
#' All comparisons between sequences performed by \code{dada} depend on pairwise alignments. This step is the most 
#'  computationally intensive part of the algorithm, and two alignment heuristics have been implemented for speed:
#'  A kmer-distance screen and banded Needleman-Wunsch alignmemt. See \code{\link{setDadaOpt}}.
#'  
#' @seealso 
#'  \code{\link{derepFastq}}, \code{\link{setDadaOpt}}
#'
#' @export
#'
#' @examples
#' derep1 = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' derep2 = derepFastq(system.file("extdata", "sam2F.fastq.gz", package="dada2"))
#' dada(derep1, err=tperr1)
#' dada(list(sam1=derep1, sam2=derep2), err=tperr1, selfConsist=TRUE)
#' dada(derep1, err=inflateErr(tperr1,3), BAND_SIZE=32, OMEGA_A=1e-20)
#'
dada <- function(derep,
                 err,
                 errorEstimationFunction = loessErrfun,
                 selfConsist = FALSE, 
                 pool = FALSE,
                 nthreads = 1, ...) {
  
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
  if(!is.list.of(derep, "derep")) { stop("The derep argument must be a derep-class object or list of derep-class objects.") }
  if(opts$USE_QUALS && any(is.null(lapply(derep, function(x) x$quals)))) { stop("The input derep-class object(s) must include quals if USE_QUALS is TRUE.") }
  
  # Validate derep object(s)
  for(i in seq_along(derep)) {
    if(!(is.integer(derep[[i]]$uniques))) {
      stop("Invalid derep$uniques vector. Must be integer valued.")
    }
    if(!(all(C_isACGT(names(derep[[i]]$uniques))))) {
      stop("Invalid derep$uniques vector. Names must be sequences made up only of A/C/G/T.")
    }
    if(sum(tabulate(nchar(names(derep[[i]]$uniques)))>0) > 1) {
      stop("Invalid derep$uniques vector. All sequences must be the same length.")
    }
  }

  # Validate quals matrix(es)
  qmax <- 0
  if(opts$USE_QUALS) {
    for(i in seq_along(derep)) {
      if(nrow(derep[[i]]$quals) != length(derep[[i]]$uniques)) {
        stop("derep$quals matrices must have one row for each derep$unique sequence.")
      }
      if(any(sapply(names(derep[[i]]$uniques), nchar) > ncol(derep[[i]]$quals))) {
        stop("derep$quals matrices must have as many columns as the length of the derep$unique sequences.")
      }
      if(any(is.na(derep[[i]]$quals))) {
        stop("NAs in derep$quals matrix. Check that all input sequences were the same length.")
      }
      if(min(derep[[i]]$quals) < 0) {
        stop("Invalid derep$quals matrix. Quality values must be positive integers.")
      }
      qmax <- max(qmax, max(derep[[i]]$quals))
    }
  }
  qmax <- ceiling(qmax) # Only getting averages from derep$quals
  if(qmax > 45) {
    if(qmax > 62) {
      stop("drep$quals matrix has an invalid maximum Phred Quality Scores of ", qmax) 
    }
    warning("derep$quals matrix has Phred Quality Scores >45. For Illumina 1.8 or earlier, this is unexpected.")
  }
  
  # Pool the derep objects if so indicated
  if(length(derep) <= 1) { pool <- FALSE }
  if(pool) { # Make derep a length 1 list of pooled derep object
    derep.in <- derep
    derep <- list(combineDereps2(derep))
  }
  
  # Validate err matrix
  initializeErr <- FALSE
  if(class(err) == "dada") { err <- err$err_out }
  if(is.null(err) && selfConsist) {
    message("Initial error matrix unspecified. Error rates will be initialized to the maximum possible estimate from this data.")
    initializeErr <- TRUE
  } else {
    if(!is.numeric(err)) stop("Error matrix must be numeric.")
    if(!(nrow(err)==16)) stop("Error matrix must have 16 rows (A2A, A2C, ...).")
    if(!all(err>=0)) stop("All error matrix entries must be >= 0.")
    if(!all(err<=1)) stop("All error matrix entries must be <=1.")
    if(any(err==0)) warning("Zero in error matrix.")
    if(ncol(err) < qmax+1) { # qmax = 0 if USE_QUALS = FALSE
      message("The supplied error matrix does not extend to maximum observed Quality Scores in derep (", qmax, ").
  Extending error rates by repeating the last column of the Error Matrix (column ", ncol(err), ").
  In selfConsist mode this should converge to the proper error rates, otherwise this may not be what you want.")
      for (q in seq(ncol(err), qmax)) { 
        err <- cbind(err, err[1:16, q])
        colnames(err)[q+1] <- q
      }
    }
  }

  # Might want to check for summed transitions from NT < 1 also.
  
  # Validate errorEstimationFunction
  if(!opts$USE_QUALS) {
    if(!is.null(errorEstimationFunction)) message("The errorEstimationFunction argument is ignored when USE_QUALS is FALSE.")
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
  
  # Validate alignment parameters
  if(opts$GAP_PENALTY>0) opts$GAP_PENALTY = -opts$GAP_PENALTY
  if(is.null(opts$HOMOPOLYMER_GAP_PENALTY)) { # Don't use homopolymer gapping
    opts$HOMOPOLYMER_GAP_PENALTY <- 99
  } else { # Use homopolymer gapping
    opts$VECTORIZED_ALIGNMENT <- FALSE # No homopolymer gapping in vectorized aligner
    if(opts$HOMOPOLYMER_GAP_PENALTY > 0) opts$HOMOPOLYMER_GAP_PENALTY = -opts$HOMOPOLYMER_GAP_PENALTY
  }
  if(opts$VECTORIZED_ALIGNMENT) {
    if(length(unique(diag(opts$SCORE)))!=1 || 
           length(unique(opts$SCORE[upper.tri(opts$SCORE) | lower.tri(opts$SCORE)]))!=1) {
      message("The vectorized aligner requires that the score matrix reduces to match/mismatch. Turning off vectorization.")
      opts$VECTORIZED_ALIGNMENT=FALSE
    }
    if(opts$BAND_SIZE > 0 && opts$BAND_SIZE<8) {
      message("Warning: The vectorized aligner is slower for very small band sizes.")
    }
    if(opts$BAND_SIZE == 0) opts$VECTORIZED_ALIGNMENT=FALSE
  }
  
  # Check nthreads
  if(nthreads < 1) { nthreads <- 1 }
  if(nthreads > 255) { nthreads <- 255 }
  
  # Initialize
  cur <- NULL
  if(initializeErr) { nconsist <- 0 } else { nconsist <- 1 }
  errs <- list()
  # The main loop, run once, or repeat until err repeats if selfConsist=T

  repeat{
    clustering <- list()
    clusterquals <- list()
    birth_subs <- list()
    trans <- list()
    map <- list()
#    exp <- list()
    prev <- cur
    if(nconsist > 0) errs[[nconsist]] <- err

    for(i in seq_along(derep)) {
      if(!opts$USE_QUALS) { qi <- matrix(0, nrow=0, ncol=0) }
      else { qi <- unname(t(derep[[i]]$quals)) } # Need transpose so that sequences are columns

      if(nconsist == 1) {
        if(pool) {
          cat(length(derep.in), "samples were pooled:", sum(derep[[i]]$uniques), "reads in", length(derep[[i]]$uniques), "unique sequences.\n")
        } else {
          cat("Sample", i, "-", sum(derep[[i]]$uniques), "reads in", length(derep[[i]]$uniques), "unique sequences.\n")
        }
      } else if(i==1) {
        if(nconsist == 0) {
          cat("Initializing error rates to maximum possible estimate.\n")
        } else {
          cat("   selfConsist step", nconsist, "\n")
        }
      }
      res <- dada_uniques(names(derep[[i]]$uniques), unname(derep[[i]]$uniques), 
                          if(initializeErr) { matrix(1, nrow=16, ncol=max(41,qmax+1)) } else { err },
                          qi, 
                          opts[["SCORE_MATRIX"]], opts[["GAP_PENALTY"]],
                          opts[["USE_KMERS"]], opts[["KDIST_CUTOFF"]],
                          opts[["BAND_SIZE"]],
                          opts[["OMEGA_A"]], 
                          if(initializeErr) { 1 } else { opts[["MAX_CLUST"]] },
                          opts[["MIN_FOLD"]], opts[["MIN_HAMMING"]],
                          opts[["USE_QUALS"]],
                          FALSE,
#                          opts[["FINAL_CONSENSUS"]],
                          opts[["VECTORIZED_ALIGNMENT"]],
                          opts[["HOMOPOLYMER_GAP_PENALTY"]],
                          nthreads,
                          opts[["VERBOSE"]])
      
      # Augment the returns
      res$clustering$sequence <- as.character(res$clustering$sequence)
      # ... nothing here for now
      
      # List the returns
      clustering[[i]] <- res$clustering
      clusterquals[[i]] <- t(res$clusterquals) # make sequences rows and positions columns
      birth_subs[[i]] <- res$birth_subs
      trans[[i]] <- res$subqual
      map[[i]] <- res$map
#      exp[[i]] <- res$exp
      rownames(trans[[i]]) <- c("A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T", "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T")
      if(opts$USE_QUALS) colnames(trans[[i]]) <- seq(0, ncol(trans[[i]])-1)  # Assumes C sides is returning one col for each integer starting at 0
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
    if(initializeErr) {
      initializeErr <- FALSE
      err[c(1,6,11,16),] <- 1.0 # Set self-transitions (A2A, C2C, G2G, T2T) to max of 1
    }

    if(selfConsist) { # Validate err matrix
      if(!is.numeric(err)) stop("Error matrix returned by errorEstimationFunction not numeric.")
      if(!(nrow(err)==16)) stop("Error matrix returned by errorEstimationFunction does not have 16 rows.")
      if(!all(err>=0)) stop("Error matrix returned by errorEstimationFunction has entries <0.")
      if(!all(err<=1)) stop("Error matrix returned by errorEstimationFunction has entries >1.")
      if(any(err==0)) warning("Error matrix returned by errorEstimationFunction has 0s in some entries.")      
    }
    
    # Termination condition for selfConsist loop
    if((!selfConsist) || any(sapply(errs, identical, err)) || (nconsist >= opts$MAX_CONSIST)) {
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
  rval2 = replicate(length(derep), list(denoised=NULL, clustering=NULL, sequence=NULL, quality=NULL, birth_subs=NULL, trans=NULL, map=NULL,
                                        err_in=NULL, err_out=NULL, opts=NULL, call=NULL), simplify=FALSE)
  for(i in seq_along(derep)) {
    rval2[[i]]$denoised <- getUniques(clustering[[i]])
    rval2[[i]]$clustering <- clustering[[i]]
    rval2[[i]]$sequence <- names(rval2[[i]]$denoised)
    rval2[[i]]$quality <- clusterquals[[i]]
    rval2[[i]]$birth_subs <- birth_subs[[i]]
    rval2[[i]]$trans <- trans[[i]]
    rval2[[i]]$map <- map[[i]]
#    rval2[[i]]$exp <- exp[[i]]
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

  # If pool=TRUE, expand the rval and prune the individual return objects
  if(pool) {
    # Expand rval into a list of the proper length
    rval1 <- rval2[[1]]
    rval2 = replicate(length(derep.in), list(denoised=NULL, clustering=NULL, sequence=NULL, quality=NULL, birth_subs=NULL, trans=NULL, map=NULL,
                                          err_in=NULL, err_out=NULL, opts=NULL, call=NULL), simplify=FALSE)
    # Make map named by the pooled unique sequence
    map <- map[[1]]
    names(map) <- names(derep[[1]]$uniques)
    for(i in seq_along(derep.in)) {
      rval2[[i]] <- rval1
      # Identify which output clusters to keep
      keep <- unique(map[names(derep[[1]]$uniques) %in% names(derep.in[[i]]$uniques)])
      keep <- seq(length(rval2[[i]]$denoised)) %in% keep # -> logical
      newBi <- cumsum(keep) # maps pooled cluster index to individual index
      # Prune $denoised, $clustering, $sequence, $quality
      rval2[[i]]$denoised <- rval2[[i]]$denoised[keep]
      rval2[[i]]$clustering <- rval2[[i]]$clustering[keep,] # Leaves old (char of integer) rownames!
      rownames(rval2[[i]]$clustering) <- as.character(newBi[as.integer(rownames(rval2[[i]]$clustering))])
      rval2[[i]]$sequence <- rval2[[i]]$sequence[keep]
      rval2[[i]]$quality <- rval2[[i]]$quality[keep,,drop=FALSE] # Not the qualities for this sample alone!
      # Prune birth_subs and remap its $clust column
      rval2[[i]]$birth_subs <- rval2[[i]]$birth_subs[keep[rval2[[i]]$birth_subs$clust],,drop=FALSE]
      rval2[[i]]$birth_subs$clust <- newBi[rval2[[i]]$birth_subs$clust]      
      # Remap $map
      rval2[[i]]$map <- newBi[map[names(derep.in[[i]]$uniques)]]
      # Recalculate abundances (both $denoised and $clustering$abundance)
      rval2[[i]]$denoised[] <- tapply(derep.in[[i]]$uniques, rval2[[i]]$map, sum)
      rval2[[i]]$clustering$abundance <- rval2[[i]]$denoised
    }
    derep <- derep.in
    rm(derep.in)
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
#' @return NULL.  
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
#' HOMOPOLYMER_GAP_PENALTY: The cost of gaps in homopolymer regions (>=3 repeated bases). Default is NULL, which causes homopolymer
#'  gaps to be treated as normal gaps.  
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
