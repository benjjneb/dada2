dada_opts <- new.env()
assign("OMEGA_A", 1e-40, envir = dada_opts)
assign("OMEGA_P", 1e-4, envir = dada_opts)
assign("OMEGA_C", 1e-40, envir=dada_opts)
assign("USE_KMERS", TRUE, envir = dada_opts)
assign("KDIST_CUTOFF", 0.42, envir = dada_opts)
assign("MAX_CONSIST", 10, envir = dada_opts)
#assign("SCORE_MATRIX", matrix(c(5L, -4L, -4L, -4L, -4L, 5L, -4L, -4L, -4L, -4L, 5L, -4L, -4L, -4L, -4L, 5L),
#                              nrow=4, byrow=TRUE), envir = dada_opts)
assign("MATCH", 5L, envir = dada_opts)
assign("MISMATCH", -4L, envir = dada_opts)
assign("GAP_PENALTY", -8L, envir = dada_opts)
assign("BAND_SIZE", 16, envir = dada_opts)
assign("VECTORIZED_ALIGNMENT", TRUE, envir = dada_opts)
assign("MAX_CLUST", 0, envir=dada_opts)
assign("MIN_FOLD", 1, envir=dada_opts)
assign("MIN_HAMMING", 1, envir=dada_opts)
assign("MIN_ABUNDANCE", 1, envir=dada_opts)
assign("USE_QUALS", TRUE, envir=dada_opts)
assign("HOMOPOLYMER_GAP_PENALTY", NULL, envir = dada_opts)
assign("SSE", 2, envir = dada_opts)
assign("GAPLESS", TRUE, envir=dada_opts)
assign("GREEDY", TRUE, envir=dada_opts)
assign("PSEUDO_PREVALENCE", 2, envir=dada_opts)
assign("PSEUDO_ABUNDANCE", Inf, envir=dada_opts)
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
#' @param err (Required). 16xN numeric matrix, or an object coercible by \code{\link{getErrors}} 
#'  such as the output of the \code{\link{learnErrors}} function.
#'    
#'  The matrix of estimated rates for each possible nucleotide transition (from sample nucleotide to read nucleotide).
#'  Rows correspond to the 16 possible transitions (t_ij) indexed such that 1:A->A, 2:A->C, ..., 16:T->T
#'  Columns correspond to quality scores. Each entry must be between 0 and 1.
#'    
#'  Typically there are 41 columns for the quality scores 0-40.
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
#'  If pool = "pseudo", the algorithm will perform pseudo-pooling between individually processed samples.
#'  
#'  This argument has no effect if only 1 sample is provided, and \code{pool} does not affect
#'   error rates, which are always estimated from pooled observations across samples.
#'   
#' @param priors (Optional). \code{character}. Default is character(0), i.e. no prior sequences.
#'  
#' The priors argument provides a set of sequences for which there is prior information suggesting they may
#'  truly exist, i.e. are not errors. The abundance p-value of dereplicated sequences that exactly match one
#'  of the priors are calculated without conditioning on presence, allowing singletons to be detected,
#'  and are compared to a reduced threshold `OMEGA_P` when forming new partitions.
#'   
#' @param multithread (Optional). Default is FALSE.
#'  If TRUE, multithreading is enabled and the number of available threads is automatically determined.   
#'  If an integer is provided, the number of threads to use is set by passing the argument on to
#'  \code{\link{setThreadOptions}}.
#'   
#' @param verbose (Optional). Default TRUE. 
#'  Print verbose text output. More fine-grained control is available by providing an integer argument.
#' \itemize{ 
#'  \item{0: Silence. No text output (same as FALSE). }
#'  \item{1: Basic text output (same as TRUE). }
#'  \item{2: Detailed text output, mostly intended for debugging. }
#' }
#'  
#' @param ... (Optional). All dada_opts can be passed in as arguments to the dada() function.
#'  See \code{\link{setDadaOpt}} for a full list and description of these options. 
#'
#' @return A \code{\link{dada-class}} object or list of such objects if a list of dereps was provided. 
#'   
#' @details
#' 
#' Briefly, \code{dada} implements a statistical test for the notion that a specific sequence was seen too many times
#'  to have been caused by amplicon errors from currently inferred sample sequences. Overly-abundant
#'  sequences are used as the seeds of new partitions of sequencing reads, and the final set of partitions
#'  is taken to represent the denoised composition of the sample. A more detailed explanation of the algorithm
#'  is found in two publications:
#' 
#' \itemize{ 
#'  \item{Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP (2016). DADA2: High resolution sample inference from Illumina amplicon data. Nature Methods, 13(7), 581-3.}
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
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom RcppParallel setThreadOptions
#' @importFrom methods as
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
                 priors = character(0),
                 multithread = FALSE, 
                 verbose=TRUE, ...) {
  
  call <- sys.call(1)
  # Read in default opts and then replace with any that were passed in to the function
  opts <- getDadaOpt()
  args <- list(...)
  # Catch the deprecated SCORE_MATRIX option
  if("SCORE_MATRIX" %in% names(args)) {
    stop("DEFUNCT: The SCORE_MATRIX option has been replaced by the MATCH/MISMATCH arguments. Please update your code.")
  }
  # Catch the defunct VERBOSE option
  if("VERBOSE" %in% names(args)) {
    stop("DEFUNCT: The VERBOSE option has been replaced by the verbose argument. Please update your code.")
  }
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) {
      opts[[opnm]] <- args[[opnm]]
    } else {
      warning(opnm, " is not a valid DADA option.")
    }
  }
  
  # Parse verbose
  if(is.logical(verbose)) {
    if(verbose == FALSE) { verbose <- 0 }
    else { verbose <- 1 }
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
  }

  # Validate quals matrix(es)
  qmax <- 0
  for(i in seq_along(derep)) {
    if(nrow(derep[[i]]$quals) != length(derep[[i]]$uniques)) {
      stop("derep$quals matrices must have one row for each derep$unique sequence.")
    }
    if(any(sapply(names(derep[[i]]$uniques), nchar) > ncol(derep[[i]]$quals))) { ###ITS
      stop("derep$quals matrices must have as many columns as the length of the derep$unique sequences.")
    }
    if(any(sapply(seq(nrow(derep[[i]]$quals)), 
                  function(row) any(is.na(derep[[i]]$quals[row,1:nchar(names(derep[[i]]$uniques)[[row]])]))))) { ###ITS
      stop("NAs in derep$quals matrix. Check that all input sequences had valid associated qualities assigned.")
    }
    if(min(derep[[i]]$quals, na.rm=TRUE) < 0) {
      stop("Invalid derep$quals matrix. Quality values must be positive integers.")
    }
    qmax <- max(qmax, max(derep[[i]]$quals, na.rm=TRUE))
  }

  qmax <- ceiling(qmax) # Only getting averages from derep$quals
  if(qmax > 250) {
    stop("derep$quals matrix has an invalid maximum Phred Quality Scores of ", qmax) 
  }

  # Get prior sequences
  priors <- getSequences(priors)
  
  # Pool the derep objects if so indicated
  pseudo <- FALSE; pseudo_priors <- character(0)
  if(length(derep) <= 1) { pool <- FALSE }
  if(is.logical(pool)) {
    if(pool) { # Make derep a length 1 list of pooled derep object
      derep.in <- derep
      derep <- list(combineDereps2(derep))
    }
  } else if(is.character(pool) && pool == "pseudo") {
    pool <- FALSE
    pseudo <- TRUE
  } else { stop("Invalid pool argument.") }
  
  # Validate err matrix
  initializeErr <- FALSE
  if(selfConsist && (missing(err) || is.null(err))) {
    err <- NULL
    initializeErr <- TRUE
  } else {
    err <- getErrors(err, enforce=TRUE)
    if(ncol(err) < qmax+1 && verbose) { # qmax = 0 if USE_QUALS = FALSE
      message("The supplied error matrix does not extend to maximum observed Quality Scores in derep (", qmax, ").
  Extending error rates by repeating the last column of the Error Matrix (column ", ncol(err), ").
  In selfConsist mode this should converge to the proper error rates, otherwise this may not be what you want.")
      for (q in seq(ncol(err), qmax)) { 
        err <- cbind(err, err[1:16, q])
        colnames(err)[q+1] <- q
      }
    }
  }
  
  # Validate OMEGA parameters
  if(opts$OMEGA_A < 0 || opts$OMEGA_A >= 1) stop("OMEGA_A must be between zero and one.")
  if(opts$OMEGA_P < 0 || opts$OMEGA_P >= 1) stop("OMEGA_P must be between zero and one.")
  if(opts$OMEGA_P < opts$OMEGA_A && length(priors) > 0) warning("OMEGA_P should generally be larger than OMEGA_A.")
  if(opts$OMEGA_C > 1e-10 && selfConsist) warning("Strict error correction (OMEGA_C < 1e-10) is not recommended when learning error rates.")
  if(opts$OMEGA_C >= 1 && selfConsist) stop("Some error correction required when learning error rates.")
  
  # Validate errorEstimationFunction
  if(!opts$USE_QUALS) {
    if(!missing(errorEstimationFunction) && verbose) message("The errorEstimationFunction argument is ignored when USE_QUALS is FALSE.")
    errorEstimationFunction <- noqualErrfun  # NULL error function has different meaning depending on USE_QUALS
  } else {
    if(!is.function(errorEstimationFunction)) stop("Must provide a function for errorEstimationFunction.")
  }
  
  # Validate alignment parameters
  if(opts$GAP_PENALTY>0) opts$GAP_PENALTY = -opts$GAP_PENALTY
  if(is.null(opts$HOMOPOLYMER_GAP_PENALTY)) { # Set gap penalties equal
    opts$HOMOPOLYMER_GAP_PENALTY <- opts$GAP_PENALTY
  }
  if(opts$HOMOPOLYMER_GAP_PENALTY > 0) opts$HOMOPOLYMER_GAP_PENALTY = -opts$HOMOPOLYMER_GAP_PENALTY

  if(opts$HOMOPOLYMER_GAP_PENALTY != opts$GAP_PENALTY) { # Use homopolymer gapping
    opts$VECTORIZED_ALIGNMENT <- FALSE # No homopolymer gapping in vectorized aligner
  }
  if(opts$VECTORIZED_ALIGNMENT) {
    if(opts$BAND_SIZE > 0 && opts$BAND_SIZE<8) {
      if(verbose) message("The vectorized aligner is slower for very small band sizes.")
    }
    if(opts$BAND_SIZE == 0) opts$VECTORIZED_ALIGNMENT=FALSE
  }
  
  # Parse multithreading argument
  if(is.logical(multithread)) {
    if(multithread==TRUE) { RcppParallel::setThreadOptions(numThreads = "auto") }
  } else if(is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
    multithread <- TRUE
  } else {
    if(verbose) message("Invalid multithread parameter. Running as a single thread.")
    multithread <- FALSE
  }

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
    pval <- list()
    prev <- cur
    if(nconsist > 0) errs[[nconsist]] <- err

    for(i in seq_along(derep)) {
      qi <- unname(t(derep[[i]]$quals)) # Need transpose so that sequences are columns

      if(nconsist == 1 && verbose) {
        if(selfConsist) {
          if(i==1) cat("selfConsist step 1 ")
          cat(".")
        } else if(pool) {
          cat(length(derep.in), "samples were pooled:", sum(derep[[i]]$uniques), "reads in", 
              length(derep[[i]]$uniques), "unique sequences.\n")
        } else {
          cat("Sample", i, "-", sum(derep[[i]]$uniques), "reads in", 
              length(derep[[i]]$uniques), "unique sequences.\n")
        }
      } else if(i==1 && verbose) {
        if(nconsist == 0) {
          cat("Initializing error rates to maximum possible estimate.\n")
        } else {
          cat("\n   selfConsist step", nconsist)
        }
      }
      # Initialize error matrix if necessary
      if(initializeErr) {
        err <- matrix(1, nrow=16, ncol=max(41,qmax+1))
      }
      res <- dada_uniques(names(derep[[i]]$uniques), unname(derep[[i]]$uniques), names(derep[[i]]$uniques) %in% c(priors, pseudo_priors),
                          err,
                          qi, 
                          opts[["MATCH"]], opts[["MISMATCH"]], opts[["GAP_PENALTY"]],
                          opts[["USE_KMERS"]], opts[["KDIST_CUTOFF"]],
                          opts[["BAND_SIZE"]],
                          opts[["OMEGA_A"]], opts[["OMEGA_P"]], opts[["OMEGA_C"]],
                          if(initializeErr) { 1 } else { opts[["MAX_CLUST"]] },
                          opts[["MIN_FOLD"]], opts[["MIN_HAMMING"]], opts[["MIN_ABUNDANCE"]],
                          TRUE, #opts[["USE_QUALS"]],
                          FALSE,
                          opts[["VECTORIZED_ALIGNMENT"]],
                          opts[["HOMOPOLYMER_GAP_PENALTY"]],
                          multithread,
                          (verbose>=2),
                          opts[["SSE"]],
                          opts[["GAPLESS"]],
                          opts[["GREEDY"]])
      
      # Augment the returns
      res$clustering$sequence <- as.character(res$clustering$sequence)

      # List the returns
      clustering[[i]] <- res$clustering
      clusterquals[[i]] <- t(res$clusterquals) # make sequences rows and positions columns
      birth_subs[[i]] <- res$birth_subs
      trans[[i]] <- res$subqual
      map[[i]] <- res$map
      pval[[i]] <- res$pval
      rownames(trans[[i]]) <- c("A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T", "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T")
      colnames(trans[[i]]) <- seq(0, ncol(trans[[i]])-1)  # Assumes C sides is returning one col for each integer starting at 0
    }
    # Accumulate the sub matrix
    cur <- Reduce("+", trans) # The only thing that changes is err(trans), so this is sufficient
    
    # Estimate the new error model (if applicable)
    if(is.null(errorEstimationFunction)) {
      err <- NULL
    } else {
      err <- tryCatch(suppressWarnings(errorEstimationFunction(cur)),
              error = function(cond) {
                if((verbose && selfConsist) || verbose >= 2) {
                  message("Error rates could not be estimated (this is usually because of very few reads).")
                }
                return(NULL)
      })
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
      if(!pseudo || (pseudo && nconsist >= 2)) { # If pseudo, must go through first (full) loop to get pseudo priors
        break
      }
    }
    
    # Get pseudo priors
    if(pseudo && nconsist >= 1) { # Don't bother if nconsist=0, i.e. max error init
      st <- makeSequenceTable(clustering)
      pseudo_priors <- colnames(st)[colSums(st>0) >= opts$PSEUDO_PREVALENCE | colSums(st) >= opts$PSEUDO_ABUNDANCE]
      rm(st)
    }
    
    nconsist <- nconsist+1
  } # repeat

  if(selfConsist && verbose) {
    cat("\n")
    if(nconsist >= opts$MAX_CONSIST) {
      message("Self-consistency loop terminated before convergence.")
    } else {
      cat("Convergence after ", nconsist, " rounds.\n")
    }
  }
  
  # Construct return object
  # A single dada-class object if one derep object provided.
  # A list of dada-class objects if multiple derep objects provided.
  rval2 = replicate(length(derep), list(denoised=NULL, clustering=NULL, sequence=NULL, quality=NULL, birth_subs=NULL, trans=NULL, map=NULL,
                                        err_in=NULL, err_out=NULL, opts=NULL), simplify=FALSE)
  for(i in seq_along(derep)) {
    rval2[[i]]$denoised <- getUniques(clustering[[i]])
    rval2[[i]]$clustering <- clustering[[i]]
    rval2[[i]]$sequence <- names(rval2[[i]]$denoised)
    rval2[[i]]$quality <- clusterquals[[i]]
    rval2[[i]]$birth_subs <- birth_subs[[i]]
    rval2[[i]]$trans <- trans[[i]]
    rval2[[i]]$map <- map[[i]]
    rval2[[i]]$pval <- pval[[i]]
    # Return the error rate(s) used as well as the final estimated error matrix
    if(selfConsist) { # Did a self-consist loop
      rval2[[i]]$err_in <- errs
    } else {
      rval2[[i]]$err_in <- errs[[1]]
    }
    rval2[[i]]$err_out <- err
    
    # Store the options that were used in the return object
    rval2[[i]]$opts <- opts
  }

  # If pool=TRUE, expand the rval and prune the individual return objects
  if(pool) {
    # Expand rval into a list of the proper length
    rval1 <- rval2[[1]]
    rval2 = replicate(length(derep.in), list(denoised=NULL, clustering=NULL, sequence=NULL, quality=NULL, birth_subs=NULL, trans=NULL, map=NULL,
                                          err_in=NULL, err_out=NULL, opts=NULL), simplify=FALSE)
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
      # Would need to add $pval back in here
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
#' @details
#' 
#' **Sensitivity**
#' 
#' OMEGA_A: This parameter sets the threshold for when DADA2 calls unique sequences significantly overabundant, and therefore creates a
#'  new partition with that sequence as the center. Default is 1e-40, which is a conservative setting to avoid making false
#'  positive inferences, but which comes at the cost of reducing the ability to identify some rare variants.
#' 
#' OMEGA_P: The threshold for unique sequences with prior evidence of existence (see `priors` argument). Default is 1e-4.
#' 
#' OMEGA_C: The threshold at which unique sequences inferred to contain errors are corrected in the final output. 
#' The probability that each unique sequence
#'  is generated at its observed abundance from the center of its final partition is evaluated, and compared to OMEGA_C. If that
#'  probability is >= OMEGA_C, it is "corrected", i.e. replaced by the partition center sequence. The special value of 0 corresponds
#'  to correcting all input sequences, and any value > 1 corresponds to performing no correction on sequences found to contain
#'  errors. Default is 1e-40 (same as OMEGA_A).
#' 
#' **Alignment**
#' 
#' MATCH: The score of a match in the Needleman-Wunsch alignment. Default is 4.
#'  
#' MISMATCH: The score of a mismatch in the Needleman-Wunsch alignment. Default is -5.
#'  
#' GAP_PENALTY: The cost of gaps in the Needleman-Wunsch alignment. Default is -8.
#'  
#' HOMOPOLYMER_GAP_PENALTY: The cost of gaps in homopolymer regions (>=3 repeated bases). Default is NULL, which causes homopolymer
#'  gaps to be treated as normal gaps.  
#'  
#' BAND_SIZE: When set, banded Needleman-Wunsch alignments are performed. Banding restricts the net cumulative number of insertion
#'  of one sequence relative to the other. The default value of BAND_SIZE is 16. If DADA is applied to sequencing technologies with 
#'  high rates of indels, such as 454 sequencing, the BAND_SIZE parameter should be increased. Setting BAND_SIZE to a negative number
#'  turns off banding (i.e. full Needleman-Wunsch).
#' 
#' **Sequence Comparison Heuristics**
#' 
#' USE_KMERS: If TRUE, a 5-mer distance screen is performed prior to performing each pairwise alignment, and if the 5mer-distance
#'  is greater than KDIST_CUTOFF, no alignment is performed. Default is TRUE.
#' 
#' KDIST_CUTOFF: The default value of 0.42 was chosen to screen pairs of sequences that differ by >10\%, and was
#'  calibrated on Illumina sequenced 16S amplicon data. The assumption is that sequences that differ by such a large
#'  amount cannot be linked by amplicon errors (i.e. if you sequence one, you won't get a read of other) and so
#'  careful (and costly) alignment is unnecessary.
#' 
#' GAPLESS: If TRUE, the ordered kmer identity between pairs of sequences is compared to their unordered
#'  overlap. If equal, the optimal alignment is assumed to be gapless. Default is TRUE.
#'  Only relevant if USE_KMERS is TRUE.
#' 
#' GREEDY: The DADA2 algorithm is not greedy, but a very restricted form of greediness can be turned
#'  on via this option. If TRUE, unique sequences with reads less than those expected to be generated
#'  by resequencing just the central unique in their partition are "locked" to that partition.
#'  Modest (~30%) speedup, and almost no impact on output. Default is TRUE.
#' 
#' **New Partition Conditions**
#' 
#' MIN_FOLD: The minimum fold-overabundance for sequences to form new partitions. Default value is 1, which means this
#'  criteria is ignored.
#'  
#' MIN_HAMMING: The minimum hamming-separation for sequences to form new partitions. Default value is 1, which means this
#'  criteria is ignored.
#'
#' MIN_ABUNDANCE: The minimum abundance for unique sequences form new partitions. Default value is 1, which means this
#'  criteria is ignored.
#'
#' MAX_CLUST: The maximum number of partitions. Once this many partitions have been created, the algorithm terminates regardless
#'  of whether the statistical model suggests more real sequence variants exist. If set to 0 this argument is ignored. Default
#'  value is 0.
#'
#' **Self Consistency**
#' 
#' MAX_CONSIST: The maximum number of steps when selfConsist=TRUE. If convergence is not reached in MAX_CONSIST steps,
#'  the algorithm will terminate with a warning message. Default value is 10.
#'  
#' **Pseudo-pooling Behavior**
#'  
#' PSEUDO_PREVALENCE: When performing pseudo-pooling, all sequence variants found in at least two
#'  samples are used as priors for a subsequent round of sample inference. 
#'  Only relevant if `pool="pseudo"`. Default is 2.
#' 
#' PSEUDO_ABUNDANCE: When performing pseudo-pooling, all denoised sequence variants with total
#'  abundance (over all samples) greater than this are used as priors for a subsequent round 
#'  of sample inference.
#'  Only relevant if `pool="pseudo"`. Default is Inf (i.e. abundance ignored for this purpose).
#' 
#' **Error Model**
#' 
#' USE_QUALS: If TRUE, the dada(...) error model takes into account the consensus quality score of the dereplicated unique sequences.
#'  If FALSE, quality scores are ignored. Default is TRUE.
#' 
#' **Technical**
#' 
#' SSE: Controls the level of explicit SSE vectorization for kmer calculations. Default 2. Maintained for development reasons,
#'    should have no impact on output.   
#'      
#' \itemize{ 
#'  \item{0: No explicit vectorization (but modern compilers will auto-vectorize the code).}
#'  \item{1: Explicit SSE2. }
#'  \item{2: Explicit, packed SSE2 using 8-bit integers. Slightly faster than SSE=1. }
#' }
#' 
#' @seealso 
#'  \code{\link{getDadaOpt}}
#'
#' @export
#'
#' @examples
#' setDadaOpt(OMEGA_A = 1e-20)
#' setDadaOpt(MATCH=1, MISMATCH=-4, GAP_PENALTY=-6)
#' setDadaOpt(GREEDY=TRUE, GAPLESS=TRUE)
#'   
setDadaOpt <- function(...) {
  opts <- getDadaOpt()
  args <- list(...)
  if(length(args)==1 && class(args[[1]]) == "list") { # Arguments were passed in as a list, as returned by getDadaOpt
    args <- args[[1]]
  }
  for(opnm in names(args)) {
    if(opnm %in% names(opts)) {
      if( (class(getDadaOpt(opnm)) == class(args[[opnm]])) ||
         (opnm == "HOMOPOLYMER_GAP_PENALTY" && # Allow numeric or NULL for HOMOPOLYMER_GAP_PENALTY
          (is.numeric(args[["HOMOPOLYMER_GAP_PENALTY"]])) || is.null(args[["HOMOPOLYMER_GAP_PENALTY"]])) ) {
        assign(opnm, args[[opnm]], envir=dada_opts)
      } else {
        warning(paste0(opnm, " not set, value provided has different class (", class(args[[opnm]]), 
                       ") then current option value (", class(getDadaOpt(opnm)), ")"))
      }
    } else {
      warning(opnm, " is not a valid DADA option.")
    }
  }
}

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
