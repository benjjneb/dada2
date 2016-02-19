################################################################################
#' Merge paired forward and reverse reads after DADA denoising.
#' 
#' This function attempts each denoised pair of forward and reverse reads, rejecting any 
#' which do not sufficiently overlap or which contain too many (>0 by default) mismatches in the
#' overlap region. Note: This function assumes that the fastq files for the 
#' forward and reverse reads were in the same order. Use of the concatenate option will
#' result in concatenating forward and reverse reads without attempting a merge/alignment step.
#' 
#' @param dadaF (Required). A \code{\link{dada-class}} object.
#'  The output of dada() function on the forward reads.
#' 
#' @param derepF (Required). A \code{\link{derep-class}} object.
#'  The derep-class object returned by derepFastq() that was used as the input to the
#'  dada-class object passed to the dadaF argument.
#'  
#'  Alternatively the map itself (derepFastq()$map) can be provided in place of the
#'  derep-class object.
#'   
#' @param dadaR (Required). A \code{\link{dada-class}} object.
#'  The output of dada() function on the reverse reads.
#' 
#' @param derepR (Required). A \code{\link{derep-class}} object.
#'  See derepF description, but for the reverse reads.
#'
#' @param minOverlap (Optional). A \code{numeric(1)} of the minimum length of the overlap (in nucleotides)
#'  required for merging the forward and reverse reads. Default is 20.
#'
#' @param maxMismatch (Optional). A \code{numeric(1)} of the maximum mismatches allowed in the overlap region.
#'  Default is 0 (i.e. only exact matches in the overlap region are accepted).
#'  
#' @param returnRejects (Optional). A \code{logical(1)}. Default is False.
#'  If true, the pairs that that were rejected based on mismatches in the overlap
#'  region are retained in the return data.frame.
#'
#' @param propagateCol (Optional). \code{character}. Default is \code{character(0)}.
#'  The mergePairs return data.frame will include copies of columns with names specified
#'  in the dada-class$clustering data.frame.
#'
#' @param justConcatenate (Optional). \code{logical(1)}, Default FALSE.
#'  If TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged,
#'    with a NNNNNNNNNN (10 Ns) spacer inserted between them.
#' 
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#'
#' @return Dataframe. One row for each unique pairing of forward/reverse denoised sequences.
#' \itemize{
#'  \item{$abundance: }{Number of reads corresponding to this forward/reverse combination.}
#'  \item{$sequence: The merged sequence.}
#'  \item{$forward: The index of the forward denoised sequence.}
#'  \item{$reverse: The index of the reverse denoised sequence.}
#'  \item{$nmatch: Number of matching nts in the overlap region.}
#'  \item{$nmismatch: Number of mismatching nts in the overlap region.}
#'  \item{$nindel: Number of indels in the overlap region.}
#'  \item{$prefer: The sequence used for the overlap region. 1=forward; 2=reverse.}
#'  \item{$accept: TRUE if overlap between forward and reverse denoised sequences was at least minOverlap and had at most maxMismatch differences.
#'          FALSE otherwise.}
#'  \item{$...: Additional columns specified in the propagateCol argument.}
#' }
#' 
#' @seealso \code{\link{derepFastq}}, \code{\link{dada}}
#' @export
#' 
#' @examples
#' derepF = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#' derepR = derepFastq(system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#' dadaF <- dada(derepF, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' dadaR <- dada(derepR, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' mergePairs(dadaF, derepF, dadaR, derepR)
#' mergePairs(dadaF, derepF, dadaR, derepR, returnRejects=TRUE, propagateCol=c("n0", "birth_ham"))
#' mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE)
#' 
mergePairs <- function(dadaF, derepF, dadaR, derepR, minOverlap = 20, maxMismatch=0, returnRejects=FALSE, propagateCol=character(0), justConcatenate=FALSE, verbose=FALSE) {
  if(class(derepF) == "derep") mapF <- derepF$map
  else mapF <- derepF
  if(class(derepR) == "derep") mapR <- derepR$map
  else mapR <- derepR
  if(!(is.integer(mapF) && is.integer(mapR))) stop("Incorrect format of derep arguments.")
  rF <- dadaF$map[mapF]
  rR <- dadaR$map[mapR]
  if(any(is.na(rF)) || any(is.na(rR))) stop("Non-corresponding maps and dada-outputs.")
  
  pairdf <- data.frame(sequence = "", abundance=0, forward=rF, reverse=rR)
  ups <- unique(pairdf) # The unique forward/reverse pairs of denoised sequences
  Funqseq <- unname(as.character(dadaF$clustering$sequence[ups$forward]))
  Runqseq <- rc(unname(as.character(dadaR$clustering$sequence[ups$reverse])))
  
  if (justConcatenate == TRUE) {
    # Simply concatenate the sequences together
    ups$sequence <- mapply(function(x,y) paste0(x,"NNNNNNNNNN", y), Funqseq, Runqseq, SIMPLIFY=FALSE);  
    ups$nmatch <- 0
    ups$nmismatch <- 0
    ups$nindel <- 0
    ups$prefer <- NA
    ups$accept <- TRUE
  } else {
    # Align forward and reverse reads.
    # Use unbanded N-W align to compare forward/reverse
    # May want to adjust align params here, but for now just using dadaOpt
    alvecs <- mapply(function(x,y) nwalign(x,y,band=-1), Funqseq, Runqseq, SIMPLIFY=FALSE)
    outs <- t(sapply(alvecs, function(x) C_eval_pair(x[1], x[2])))
    ups$nmatch <- outs[,1]
    ups$nmismatch <- outs[,2]
    ups$nindel <- outs[,3]
    ups$prefer <- 1 + (dadaR$clustering$n0[ups$reverse] > dadaR$clustering$n0[ups$forward])
    ups$accept <- (ups$nmatch > minOverlap) & ((ups$nmismatch + ups$nindel) <= maxMismatch)
    # Make the sequence
    ups$sequence <- mapply(C_pair_consensus, sapply(alvecs,`[`,1), sapply(alvecs,`[`,2), ups$prefer);
    # Additional param to indicate whether 1:forward or 2:reverse takes precedence
    # Must also strip out any indels in the return
    # This function is only used here.
  }
  
  # Add abundance and sequence to the output data.frame
  tab <- table(pairdf$forward, pairdf$reverse)
  ups$abundance <- tab[cbind(ups$forward, ups$reverse)]
  ups$sequence[!ups$accept] <- ""
  # Add columns from forward/reverse clustering
  propagateCol <- propagateCol[propagateCol %in% colnames(dadaF$clustering)]
  for(col in propagateCol) {
    ups[,paste0("F.",col)] <- dadaF$clustering[ups$forward,col]
    ups[,paste0("R.",col)] <- dadaR$clustering[ups$reverse,col]
  }
  # Sort output by abundance and name
  ups <- ups[order(ups$abundance, decreasing=TRUE),]
  rownames(ups) <- paste0("s", ups$forward, "_", ups$reverse)
  if(verbose) {
    message(sum(ups$abundance[ups$accept]), " paired-reads (in ", sum(ups$accept), " unique pairings) successfully merged out of ", sum(ups$abundance), " (in ", nrow(ups), " pairings) input.")
  }
  if(!returnRejects) { ups <- ups[ups$accept,] }
  
  if(any(duplicated(ups$sequence))) {
    message("Duplicate sequences in merged output.")
  }
  
  ups
}

#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead id
#' @importFrom ShortRead yield
sameOrder <- function(fnF, fnR) {
  matched <- TRUE
  fF <- FastqStreamer(fnF)
  on.exit(close(fF))
  fR <- FastqStreamer(fnR)
  on.exit(close(fR), add=TRUE)
  
  while( length(suppressWarnings(fqF <- yield(fF))) && length(suppressWarnings(fqR <- yield(fR))) ) {
    idF <- trimTails(id(fqF), 1, " ")
    idR <- trimTails(id(fqR), 1, " ")
    matched <- matched && all(idF == idR)
  }
  return(matched)
}

isMatch <- function(al, minOverlap, verbose=FALSE) {
  out <- C_eval_pair(al[1], al[2]) # match, mismatch, indel
  if(verbose) { cat("Match/mismatch/indel:", out, "\n") }
  if(out[1] >= minOverlap && out[2] == 0 && out[3] == 0) {
    return(TRUE);
  } else {
    return(FALSE);
  }
}

################################################################################
#' Map denoised sequence to each read.
#' 
#' Takes the
#' \code{\link{dada}} result,
#' \code{\link{derepFastq}} result, and
#' \code{\link{ShortReadQ-class}} object
#' (an R representation of the original input sequence data),
#' and creates a \code{\link{data.table}}
#' in which each entry (row) contains
#' a unique read ID and the denoised sequence
#' to which it corresponds, as inferred by \code{\link{dada}}.
#' 
#' @param dadaRes (Required). A \code{\link{dada-class}} object.
#'  The output of \code{\link{dada}} function.
#' 
#' @param derep (Required).
#'  A \code{\link{derep-class}} object.
#'  An object returned by \code{\link{derepFastq}()}.
#'  Should be the same object that served as input
#'  to the \code{\link{dada}()} function.
#'  
#' @param sr (Required). The trimmed and filtered reads 
#'  that you used as input for \code{\link{derepFastq}},
#'  prior to running \code{\link{dada}()} on \code{derep}.
#'  More generally, this is an object that inherits from the 
#'  \code{\link{ShortRead-class}}. 
#'  In most cases this will be \code{\link{ShortReadQ-class}}.
#'  Objects from this class are the result of \code{\link{readFastq}}.
#'  Alternatively, this can be a character string
#'  that provides the path to your forward reads fastq file.
#'  
#' @param idRegExpr (Optional). 
#'  Exact same as for \code{\link{mergePairsByID}}.
#'
#' @param includeCol (Optional). 
#'  Exact same as for \code{\link{mergePairsByID}}.
#'  
#' @return A \code{\link{data.table}}
#' in which each entry (row) contains
#' a unique read ID and the denoised sequence
#' to which it corresponds, as inferred by \code{\link{dada}}.
#' 
#' @importFrom ShortRead id
#' @importFrom ShortRead readFastq
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom data.table uniqueN
#' 
#' @examples
#' exFileF = system.file("extdata", "sam1F.fastq.gz", package="dada2")
#' show(exFileF)
#' srF = ShortRead::readFastq(exFileF)
#' derepF = derepFastq(exFileF)
#' dadaF <- dada(derepF, err=tperr1,
#'               errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' dada2:::dada_to_seq_table(dadaF, derepF, srF)
#' 
dada_to_seq_table = function(dadaRes, derep, sr, 
                             idRegExpr = c("\\s.+$", ""), 
                             includeCol = character(0)){
  # derep <-> sr sanity check
  if( !inherits(derep, "derep") ){
    stop("`derep` argument must be of the derep class.
         See ?'derep-class' for more details.")
  }
  if( length(derep$map) != length(sr) ){
    stop("`derep` and `sr` arguments do not indicate the same number of reads.
         Check origin and resolve before trying again.")
  }
  
  map = derep$map
  if(!is.integer(map)){
    stop("Something wrong with `map`.")
  }
  # Define the read IDs
  readIDs = gsub(idRegExpr[1], idRegExpr[2], as.character(id(sr)))
  # Genotype index
  genotypeID <- dadaRes$map[map]
  if(any(is.na(genotypeID))){
    stop("Non-corresponding maps and dada-outputs.")
  }
  # A data.table containing these organized results, and more.
  dt = data.table(uniqueIndex = map, 
                  genotypeIndex = genotypeID,
                  # dereplicated aka unique sequence
                  derepSeq = names(derep$uniques)[map],
                  # genotype sequence
                  seq = names(dadaRes$denoised)[genotypeID],
                  # read id
                  id = readIDs)
  # Join with some information from `clustering`
  clusteringdt = data.table(dadaRes$clustering)
  setnames(clusteringdt, "sequence", "seq")
  # Keep only requested columns in `propogateCol`,
  # and the requisite `seq` and `n0` cols.
  clusteringdt <- clusteringdt[, c("seq", "n0", includeCol), with = FALSE]
  setkey(clusteringdt, seq)
  setkey(dt, seq)
  if(uniqueN(clusteringdt) != uniqueN(dt)){
    stop("Problem matching denoised genotypes to $clustering metadata")
  }
  dt2 = clusteringdt[dt]
  # Re-key back to `id` before handoff
  setkey(dt2, id)
  return(dt2)
}
################################################################################
#' Merge forward and reverse reads after DADA denoising,
#' even if reads were not originally ordered together.
#' 
#' This function attempts to merge each pair of denoised forward and reverse reads,
#' rejecting any which do not sufficiently overlap 
#' or which contain too many (>0 by default) mismatches in the overlap region. 
#' Note: This function does not assume that the fastq files 
#' for the forward and reverse reads were in the same order. 
#' If they are already in the same order, use \code{\link{mergePairs}}.
#' 
#' Not yet implemented: 
#' Use of the concatenate option 
#' will result in concatenating forward and reverse reads 
#' without attempting a merge/alignment step.
#' 
#' 
#' @param dadaF (Required). A \code{\link{dada-class}} object.
#'  The output of dada() function on the forward reads.
#' 
#' @param derepF (Required). A \code{\link{derep-class}} object.
#'  The derep-class object returned by derepFastq() that was used as the input to the
#'  dada-class object passed to the dadaF argument.
#'  
#'  @param srF (Required). The trimmed and filtered forward reads 
#'   that you used as input for \code{\link{derepFastq}}.
#'   More generally, this is an object that inherits from the 
#'   \code{\link{ShortRead-class}}. 
#'   In most cases this will be \code{\link{ShortReadQ-class}}.
#'   Objects from this class are the result of \code{\link{readFastq}}.
#'   Alternatively, this can be a character string
#'   that provides the path to your forward reads fastq file.
#'   
#' @param dadaR (Required). A \code{\link{dada-class}} object.
#'  The output of dada() function on the reverse reads.
#' 
#' @param derepR (Required). A \code{\link{derep-class}} object.
#'  See derepF description, but for the reverse reads.
#'  
#'  @param srR (Required). 
#'   See srF description, but in this case provide for the reverse reads.
#'
#' @param minOverlap (Optional). A \code{numeric(1)} of the minimum length of the overlap (in nucleotides)
#'  required for merging the forward and reverse reads. Default is 20.
#'
#' @param maxMismatch (Optional). A \code{numeric(1)} of the maximum mismatches allowed in the overlap region.
#'  Default is 0 (i.e. only exact matches in the overlap region are accepted).
#'  
#' @param returnRejects (Optional).
#'  A \code{\link{logical}(1)}. Default is \code{FALSE}.
#'  If \code{TRUE}, the pairs that that were rejected
#'  based on mismatches in the overlap
#'  region are retained in the return \code{\link{data.frame}}.
#'  
#' @param idRegExpr (Optional).
#'  A length 2 \code{\link{character}()} vector.
#'  This is passed along in order as the first two arguments
#'  to a \code{\link{gsub}} call that defines
#'  how each read \code{\link[ShortRead]{id}} is parsed.
#'  The default is \code{c("\\s.+$", "")},
#'  which is a \code{\link{gsub}} directive to keep 
#'  the \code{id} string from the beginning 
#'  up to but not including the first space.
#'  For some sequencing platforms and/or read ID schemes,
#'  an alternative parsing of the IDs may be appropriate.
#'
#' @param includeCol (Optional). \code{character}. 
#'  Default is \code{character(0)}.
#'  The returned \code{\link{data.table}} 
#'  will include columns with names specified
#'  by the \code{\link{dada-class}$clustering} data.frame.
#'
#' @param justConcatenate (Optional). 
#'  NOT CURRENTLY SUPPORTED.
#'  \code{logical(1)}, Default FALSE.
#'  If TRUE, the forward and reverse-complemented reverse read 
#'  are concatenated rather than merged,
#'  with a NNNNNNNNNN (10 Ns) spacer inserted between them.
#' 
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
#'
#' @return A \code{\link{data.table}}. 
#'  One row for each unique pairing of forward/reverse denoised sequences.
#' \itemize{
#'  \item{$abundance: }{Number of reads corresponding to this forward/reverse combination.}
#'  \item{$sequence: The merged sequence.}
#'  \item{$forward: The index of the forward denoised sequence.}
#'  \item{$reverse: The index of the reverse denoised sequence.}
#'  \item{$nmatch: Number of matching nts in the overlap region.}
#'  \item{$nmismatch: Number of mismatching nts in the overlap region.}
#'  \item{$nindel: Number of indels in the overlap region.}
#'  \item{$prefer: The sequence used for the overlap region. 1=forward; 2=reverse.}
#'  \item{$accept: TRUE if overlap between forward and reverse denoised sequences was at least minOverlap and had at most maxMismatch differences.
#'          FALSE otherwise.}
#'  \item{$...: Additional columns specified in the includeCol argument.}
#' }
#' 
#' @seealso \code{\link{derepFastq}}, \code{\link{dada}}
#' 
#' @importFrom ShortRead id
#' @importFrom ShortRead readFastq
#' @importFrom data.table data.table
#' @importFrom data.table uniqueN
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom data.table .N
#' 
#' @examples
#' # For the following example files, there are two ways to merge denoised directions.
#' # Because the read sequences are in order, `mergePairs()` works.
#' # `mergePairsByID` always works,
#' # because it uses the read IDs to match denoised pairs.
#' exFileF = system.file("extdata", "sam1F.fastq.gz", package="dada2")
#' exFileR = system.file("extdata", "sam1R.fastq.gz", package="dada2")
#' srF = ShortRead::readFastq(exFileF)
#' srR = ShortRead::readFastq(exFileR)
#' derepF = derepFastq(exFileF)
#' derepR = derepFastq(exFileR)
#' dadaF <- dada(derepF, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' dadaR <- dada(derepR, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
#' # Run and compare
#' ex1time = system.time({
#' ex1 <- mergePairs(dadaF, derepF, dadaR, derepR, verbose = TRUE)
#'     ex1 <- data.table::data.table(ex1)
#'  })
#' ex1time
#' # The new function, based on read IDs.
#' ex2time = system.time({
#'   ex2 = dada2:::mergePairsByID(dadaF = dadaF, derepF = derepF, srF = srF,
#'                        dadaR = dadaR, derepR = derepR, srR = srR, verbose = TRUE)
#' })
#' ex2time
#' # Compare results (should be identical)
#' ex2[(accept)]
#' data.table::setkey(ex2, sequence)
#' ex2[(accept), list(abundance = sum(abundance)), by = sequence]
#' # Same sequence set (exactly)
#' setequal(x = ex1$sequence,
#'          y = ex2[(accept)]$sequence)
#' # Test concatenation functionality
#' ex1cattime = system.time({
#' ex1cat <- mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate = TRUE, verbose = TRUE)
#' sapply(ex1cat, class)
#'   # need to convert to a character
#'   ex1cat$sequence <- unlist(ex1cat$sequence)
#'   ex1cat <- data.table::data.table(ex1cat)
#' })
#' ex1cattime
#' ex2cattime = system.time({
#'   ex2cat <- dada2:::mergePairsByID(dadaF = dadaF, derepF = derepF, srF = srF,
#'                            dadaR = dadaR, derepR = derepR, srR = srR,
#'                            justConcatenate = TRUE, verbose = TRUE)
#' })
#' ex2cattime
#' ex2cat[(accept)]
#' # Compare results (should be identical)
#' data.table::setkey(ex1cat, sequence)
#' ex1cat[(accept), list(abundance = sum(abundance)), by = sequence]
#' data.table::setkey(ex2cat, sequence)
#' ex2cat[(accept), list(abundance = sum(abundance)), by = sequence]
#' # Same sequence set (exactly)
#' setequal(x = ex1cat$sequence,
#'          y = ex2cat$sequence)
#' intersect(x = ex1cat$sequence,
#'           y = ex2cat$sequence)
#' ex1cat[, nchar(sequence)]
#' ex2cat[, nchar(sequence)]
mergePairsByID = function(dadaF, derepF, srF,
                          dadaR, derepR, srR, 
                          minOverlap = 20, 
                          maxMismatch = 0, 
                          returnRejects = FALSE,
                          idRegExpr = c("\\s.+$", ""),
                          includeCol = character(0), 
                          justConcatenate = FALSE,
                          verbose = FALSE) {
  # Interpret reads object | path
  # Forward Reads
  if(inherits(srF, "character")){
    if(verbose){
      message("`srF` interpreted as path to forward reads fastq file.
              Attempting to read...")
    }
    srF <- readFastq(srF)
    }
  if(!inherits(srF, "ShortRead")){
    stop("`srF` was not properly read, or is wrong object type.
         Please check documentation and try again.")
  }
  # Reverse Reads
  if(inherits(srR, "character")){
    if(verbose){
      message("`srR` interpreted as path to forward reads fastq file.
              Attempting to read...")
    }
    srR <- readFastq(srR)
    }
  if(!inherits(srR, "ShortRead")){
    stop("`srR` was not properly read, or is wrong object type.
         Please check documentation and try again.")
  }
  # Map read IDs to denoised sequence, `n0`, and optional columns
  dtF = dada_to_seq_table(dadaF, derepF, srF, 
                          idRegExpr = idRegExpr,
                          includeCol = includeCol)
  dtR = dada_to_seq_table(dadaR, derepR, srR,
                          idRegExpr = idRegExpr,
                          includeCol = includeCol)
  if(verbose){message(uniqueN(dtF), " unique forward read IDs.")}
  if(verbose){message(uniqueN(dtR), " unique reverse read IDs.")}
  # Rename the optional propagated columns before join
  if(length(includeCol) > 0){
    includeColF = paste0(includeCol, "_F")
    includeColR = paste0(includeCol, "_R")
    setnames(dtF, includeCol, includeColF)
    setnames(dtR, includeCol, includeColR)
    includeCol <- c(includeColF, includeColR)
  } else {
    includeColF = includeColR = includeCol
  }
  # Omit unneeded columns before join
  dtF <- dtF[, c("seq", "id", "n0", includeColF), with = FALSE]
  dtR <- dtR[, c("seq", "id", "n0", includeColR), with = FALSE]
  # Reverse complement the reverse reads
  dtR[, seq := rc(seq)]
  # Sanity check
  if(nrow(dtF[, .N, by = seq]) != length(dadaF$denoised)){
    stop("Forward read mapping to denoised sequence index failed.")
  }
  if(nrow(dtR[, .N, by = seq]) != length(dadaR$denoised)){
    stop("Reverse read mapping to denoised sequence index failed.")
  } 
  # Rename columns to indicate read direction
  setnames(dtF, "seq", "seqF")
  setnames(dtR, "seq", "seqR")
  setnames(dtF, "n0", "n0F")
  setnames(dtR, "n0", "n0R")  
  # Join each dada results mapping table to create ID data.table, `iddt`
  # `nomatch = 0` means ids that are not present in both will be dropped.
  setkey(dtF, id)
  setkey(dtR, id)
  iddt = dtF[dtR, nomatch = 0]
  # Define the unique version of `iddt` table.
  # Unique is defined by seqF and seqR, the pair of denoised sequences.
  # Each unique pair of denoised sequences must be evaluated independently
  # for overlap acceptance. However, we don't want to do this more than once per pair.
  # Will map back to iddt later to fill out the original table
  setkey(iddt, seqF, seqR)
  # Tally abundance of each pair
  # (still long form, nrow == nreads)
  iddt[, abundance := .N, by = list(seqF, seqR)]
  # Unique Pair iddt
  upiddt = unique(iddt)
  setkey(upiddt, seqF, seqR)
  if(verbose){
    message(nrow(iddt), " paired reads, corresponding to ",
            nrow(upiddt), " unique pairs that must be assessed for overlap merge")
  }
  # If `justConcatenate` is TRUE, don't need to align or evaluate.
  # Can have an early return
  if(justConcatenate){
    upiddt[, sequence := paste0(seqF, rep("N", times = 10), seqR), by = list(seqF, seqR)]
    upiddt[, accept := TRUE]
  } else {
    # Otherwise, alignment needed. More information considered and returned.
    # for functions that return multiple values...
    # DT[, c("new1","new2") := myfun(y,v)]
    # http://stackoverflow.com/questions/11308754/add-multiple-columns-to-r-data-table-in-one-function-call
    # (1) run nw unbanded alignment
    # `als` - alignment sequence, 1 forward, 2 reverse
    upiddt[, c("als1", "als2") := as.list(nwalign(seqF, seqR, band=-1)),
           by = list(seqF, seqR)]
    # (2) Evaluate alignment
    upiddt[, c("match", "mismatch", "indel") := as.list(C_eval_pair(als1, als2)),
           by = list(seqF, seqR)]
    upiddt[, prefer := 1 + (n0R > n0F)]
    upiddt[, allMismatch := mismatch + indel,
           by = list(seqF, seqR)]
    # (3) Define Acceptance
    upiddt[, accept := (match > minOverlap) & (allMismatch <= maxMismatch),
           by = list(seqF, seqR)]
    # (4) Make the consensus sequence, 
    #     but only for the pairs that passed (accepted merges)
    upiddt[(accept), sequence := C_pair_consensus(als1, als2, prefer),
           by = list(seqF, seqR)]
  }
  # Optionally add column details from iddt table (includeCol)
  if( !is.null(includeCol) ){
    iddtProp = iddt[, c("seqF", "seqR", includeCol), with = FALSE]
    setkey(iddtProp, seqF, seqR)
    iddtProp <- unique(iddtProp)
    setkey(upiddt, seqF, seqR)
    upiddt <- iddtProp[upiddt]
  }
  if(verbose) {
    # Report summary
    message(sum(upiddt[(accept)]$abundance),
            " paired-reads (in ", 
            nrow(upiddt[(accept)]), " unique pairings) successfully merged\n",
            "from ", nrow(iddt), " read pairs.")
  }
  # ID column doesn't make sense to include
  if("id" %in% colnames(upiddt)){upiddt[, id := NULL]}
  return(upiddt)
}
