#'
#' Classifies sequences against reference training dataset.
#' 
#' assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm described in
#' Wang et al. Applied and Environmental Microbiology 2007, with kmer size 8 and 100 bootstrap
#' replicates. Properly formatted reference files for several popular taxonomic databases
#' are available \url{http://benjjneb.github.io/dada2/training.html}
#' 
#' @param seqs (Required). A character vector of the sequences to be assigned, or an object 
#' coercible by \code{\link{getUniques}}.
#'   
#' @param refFasta (Required). The path to the reference fasta file, or an 
#' R connection Can be compressed.
#' This reference fasta file should be formatted so that the id lines correspond to the
#' taxonomy (or classification) of the associated sequence, and each taxonomic level is 
#' separated by a semicolon. Eg.
#' 
#'  >Kingom;Phylum;Class;Order;Family;Genus;   
#'  ACGAATGTGAAGTAA......   
#' 
#' @param minBoot (Optional). Default 50. 
#' The minimum bootstrap confidence for assigning a taxonomic level.
#'   
#' @param tryRC (Optional). Default FALSE. 
#' If TRUE, the reverse-complement of each sequences will be used for classification if it is a better match to the reference
#' sequences than the forward sequence.
#'   
#' @param outputBootstraps (Optional). Default FALSE.
#'  If TRUE, bootstrap values will be retained in an integer matrix. A named list containing the assigned taxonomies (named "taxa") 
#'  and the bootstrap values (named "boot") will be returned. Minimum bootstrap confidence filtering still takes place,
#'  to see full taxonomy set minBoot=0
#'   
#' @param taxLevels (Optional). Default is c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species").
#' The taxonomic levels being assigned. Truncates if deeper levels not present in
#' training fasta.
#'   
#' @param multithread (Optional). Default is FALSE.
#'  If TRUE, multithreading is enabled and the number of available threads is automatically determined.   
#'  If an integer is provided, the number of threads to use is set by passing the argument on to
#'  \code{\link{setThreadOptions}}.
#'   
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, print status to standard output.
#'   
#' @return A character matrix of assigned taxonomies exceeding the minBoot level of
#'   bootstrapping confidence. Rows correspond to the provided sequences, columns to the
#'   taxonomic levels. NA indicates that the sequence was not consistently classified at
#'   that level at the minBoot threshhold.
#'   
#'   If outputBootstraps is TRUE, a named list containing the assigned taxonomies (named "taxa") 
#'   and the bootstrap values (named "boot") will be returned.
#' 
#' @export
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead id
#' 
#' @examples
#' seqs <- getSequences(system.file("extdata", "example_seqs.fa", package="dada2"))
#' training_fasta <- system.file("extdata", "example_train_set.fa.gz", package="dada2")
#' taxa <- assignTaxonomy(seqs, training_fasta)
#' taxa80 <- assignTaxonomy(seqs, training_fasta, minBoot=80, multithread=2)
#' 
assignTaxonomy <- function(seqs, refFasta, minBoot=50, tryRC=FALSE, outputBootstraps=FALSE,
                           taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                           multithread=FALSE, verbose=FALSE) {
  MIN_REF_LEN <- 20 # Enforced minimum length of reference seqs. Must be bigger than the kmer-size used (8).
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  # Read in the reference fasta
  refsr <- readFasta(refFasta)
  lens <- width(sread(refsr))
  if(any(lens<MIN_REF_LEN)) {
    refsr <- refsr[lens>=MIN_REF_LEN]
    warning(paste0("Some reference sequences were too short (<", MIN_REF_LEN, "nts) and were excluded."))
  }
  refs <- as.character(sread(refsr))
  tax <- as.character(id(refsr))
  tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x)) # Remove leading/trailing whitespace
  # Sniff and parse UNITE fasta format
  UNITE <- FALSE
  if(all(grepl("FU\\|re[pf]s", tax[1:10]))) {
    UNITE <- TRUE
    cat("UNITE fungal taxonomic reference detected.\n")
    tax <- sapply(strsplit(tax, "\\|"), `[`, 5)
    tax <- gsub("[pcofg]__unidentified;", "_DADA2_UNSPECIFIED;", tax)
    tax <- gsub(";s__(\\w+)_", ";s__", tax)
    tax <- gsub(";s__sp$", ";_DADA2_UNSPECIFIED", tax)
  }
  # Crude format check
  if(!grepl(";", tax[[1]])) {
    if(length(unlist(strsplit(tax[[1]], "\\s")))==3) {
      stop("Incorrect reference file format for assignTaxonomy (this looks like a file formatted for assignSpecies).")
    } else {
      stop("Incorrect reference file format for assignTaxonomy.")
    }
  }
  # Parse the taxonomies from the id string
  tax.depth <- sapply(strsplit(tax, ";"), length)
  td <- max(tax.depth)
  for(i in seq(length(tax))) {
    if(tax.depth[[i]] < td) {
      for(j in seq(td - tax.depth[[i]])) {
        tax[[i]] <- paste0(tax[[i]], "_DADA2_UNSPECIFIED;")
      }
    }
  }
  # Create the integer maps from reference to type ("genus") and for each tax level
  genus.unq <- unique(tax)
  ref.to.genus <- match(tax, genus.unq)
  tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol=td, byrow=TRUE)
  tax.df <- as.data.frame(tax.mat)
  for(i in seq(ncol(tax.df))) {
    tax.df[,i] <- factor(tax.df[,i])
    tax.df[,i] <- as.integer(tax.df[,i])
  }
  tax.mat.int <- as.matrix(tax.df)
  # Assign
  # Parse multithreading argument
  if(is.logical(multithread)) {
    if(multithread==TRUE) { RcppParallel::setThreadOptions(numThreads = "auto") }
  } else if(is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
    multithread <- TRUE
  } else {
    warning("Invalid multithread parameter. Running as a single thread.")
    multithread <- FALSE
  }
  if(multithread) {
    assignment <- C_assign_taxonomy2(seqs, rc(seqs), refs, ref.to.genus, tax.mat.int, tryRC, verbose)
  } else {
    assignment <- C_assign_taxonomy(seqs, rc(seqs), refs, ref.to.genus, tax.mat.int, tryRC, verbose)
  }
  # Parse results and return tax consistent with minBoot
  bestHit <- genus.unq[assignment$tax]
  boots <- assignment$boot
  taxes <- strsplit(bestHit, ";")
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,]>=minBoot])
  # Convert to character matrix
  tax.out <- matrix(NA_character_, nrow=length(seqs), ncol=td)
  for(i in seq(length(seqs))) {
    if(length(taxes[[i]]) > 0) {
      tax.out[i,1:length(taxes[[i]])] <- taxes[[i]]
    }
  }
  rownames(tax.out) <- seqs
  colnames(tax.out) <- taxLevels[1:ncol(tax.out)]
  tax.out[tax.out=="_DADA2_UNSPECIFIED"] <- NA_character_
  if(outputBootstraps){
      # Convert boots to integer matrix
      boots.out <- matrix(boots, nrow=length(seqs), ncol=td)
      rownames(boots.out) <- seqs
      colnames(boots.out) <- taxLevels[1:ncol(boots.out)]
      list(tax=tax.out, boot=boots.out)
  } else {
    tax.out
  }
}

# Helper function for assignSpecies
mapHits <- function(x, refs, keep, sep="/") {
  hits <- refs[x]
  hits[grepl("Escherichia", hits, fixed=TRUE) | grepl("Shigella", hits, fixed=TRUE)] <- "Escherichia/Shigella"
  if(length(unique(hits))<=keep) {
    rval <- do.call(paste, c(as.list(sort(unique(hits))), sep=sep))
  } else { rval <- NA_character_ }
  if(length(rval)==0) rval <- NA_character_
  rval
}

# Match curated genus names to binomial genus names
# Handles Clostridium groups and split genera names
matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
  if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
  if((gen.tax==gen.binom) || 
     grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) || 
     grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#'
#' Taxonomic assignment to the species level by exact matching.
#' 
#' \code{assignSpecies} uses exact matching against a reference fasta to identify the 
#' genus-species binomial classification of the input sequences.
#' 
#' @param seqs (Required). A character vector of the sequences to be assigned, or an object 
#' coercible by \code{\link{getUniques}}.
#'   
#' @param refFasta (Required). The path to the reference fasta file, or an 
#' R connection. Can be compressed.
#' This reference fasta file should be formatted so that the id lines correspond to the
#' genus-species of the associated sequence:
#'   
#'  >SeqID genus species  
#'  ACGAATGTGAAGTAA......
#' 
#' @param allowMultiple (Optional). Default FALSE.
#' Defines the behavior when multiple exact matches against different species are returned.
#' By default only unambiguous identifications are return. If TRUE, a concatenated string
#' of all exactly matched species is returned. If an integer is provided, multiple
#' identifications up to that many are returned as a concatenated string.
#'   
#' @param tryRC (Optional). Default FALSE. 
#' If TRUE, the reverse-complement of each sequences will also be tested for exact matching 
#' to the reference sequences.
#'   
#' @param n (Optional). Default \code{2000}.
#' The number of sequences to perform assignment on at one time. 
#' This controls the peak memory requirement so that large numbers of sequences are supported. 
#'
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, print status to standard output.
#' 
#' @return A two-column character matrix. Rows correspond to the provided sequences,
#'   columns to the genus and species taxonomic levels. NA indicates that the sequence
#'   was not classified at that level. 
#' 
#' @export
#' 
#' @importFrom Biostrings vcountPDict
#' @importFrom Biostrings PDict
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead reverseComplement
#' @importFrom ShortRead id
#' @importFrom methods as
#' 
#' @examples
#' seqs <- getSequences(system.file("extdata", "example_seqs.fa", package="dada2"))
#' species_fasta <- system.file("extdata", "example_species_assignment.fa.gz", package="dada2")
#' spec <- assignSpecies(seqs, species_fasta)
#' 
assignSpecies <- function(seqs, refFasta, allowMultiple=FALSE, tryRC=FALSE, n=2000, verbose=FALSE) {
  # Define number of multiple species to return
  if(is.logical(allowMultiple)) {
    if(allowMultiple) keep <- Inf
    else keep <- 1
  } else {
    keep <- as.integer(allowMultiple)
  }
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  # Read in the reference fasta
  refsr <- readFasta(refFasta)
  ids <- as(id(refsr), "character")
  # Crude format check
  if(!length(unlist(strsplit(ids[[1]], "\\s"))) >= 3) {
    if(length(unlist(gregexpr(";", ids[[1]]))) >= 3) {
      stop("Incorrect reference file format for assignSpecies (this looks like a file formatted for assignTaxonomy).")
    } else {
      stop("Incorrect reference file format for assignSpecies.")
    }
  }
  genus <- sapply(strsplit(ids, "\\s"), `[`, 2)
  species <- gsub("^\\S*\\s\\S*\\s", "", ids)
  # Identify the exact hits
  hits <- vector("list", length(seqs))
  lens <- nchar(seqs)
  for(len in unique(lens)) { # Requires all same length sequences
    i.len <- which(lens==len); n.len <- length(i.len)
    j.lo<-1; j.hi<-min(n,n.len)
    while(j.lo <= n.len) {
      i.loop <- i.len[j.lo:j.hi]
      seqdict <- PDict(seqs[i.loop])
      vhit <- (vcountPDict(seqdict, sread(refsr))>0)
      if(tryRC) vhit <- vhit | (vcountPDict(seqdict, reverseComplement(sread(refsr)))>0)
      hits[i.loop] <- lapply(seq(nrow(vhit)), function(x) vhit[x,])
      j.lo <- j.lo + n; j.hi <- min(j.hi+n, n.len)
      rm(seqdict)
      gc()
    }
  }
  # Get genus species return strings
  rval <- cbind(unlist(sapply(hits, mapHits, refs=genus, keep=1)),
                unlist(sapply(hits, mapHits, refs=species, keep=keep)))
  colnames(rval) <- c("Genus", "Species")
  rownames(rval) <- seqs
  gc()
  if(verbose) cat(sum(!is.na(rval[,"Species"])), "out of", length(seqs), "were assigned to the species level.\n")
  rval
}

#'
#' Add species-level annotation to a taxonomic table.
#' 
#' \code{addSpecies} wraps the \code{\link{assignSpecies}} function to assign genus-species 
#' binomials to the input sequences by exact matching against a reference fasta. Those binomials
#' are then merged with the input taxonomic table with species annotations appended as an 
#' additional column to the input table.
#' Only species identifications where the genera in the input table and the binomial 
#' classification are consistent are included in the return table.
#' 
#' @param taxtab (Required). A taxonomic table, the output of \code{\link{assignTaxonomy}}.
#'   
#' @param refFasta (Required). The path to the reference fasta file, or an 
#' R connection. Can be compressed.
#' This reference fasta file should be formatted so that the id lines correspond to the
#' genus-species binomial of the associated sequence:
#'   
#'  >SeqID genus species  
#'  ACGAATGTGAAGTAA......
#' 
#' @param allowMultiple (Optional). Default FALSE.
#' Defines the behavior when multiple exact matches against different species are returned.
#' By default only unambiguous identifications are return. If TRUE, a concatenated string
#' of all exactly matched species is returned. If an integer is provided, multiple
#' identifications up to that many are returned as a concatenated string.
#'   
#' @param tryRC (Optional). Default FALSE. 
#' If TRUE, the reverse-complement of each sequences will be used for classification if it is a better match to the reference
#' sequences than the forward sequence.
#'   
#' @param n (Optional). Default \code{1e5}.
#' The number of records (reads) to read in and filter at any one time. 
#' This controls the peak memory requirement so that very large fastq files are supported. 
#' See \code{\link{FastqStreamer}} for details.
#'
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, print status to standard output.
#' 
#' @return A character matrix one column larger than input. Rows correspond to
#'   sequences, and columns to the taxonomic levels. NA indicates that the sequence
#'   was not classified at that level. 
#' 
#' @seealso 
#'  \code{\link{assignTaxonomy}}, \code{\link{assignSpecies}}
#'  
#' @export
#' 
#' @examples
#' 
#' seqs <- getSequences(system.file("extdata", "example_seqs.fa", package="dada2"))
#' training_fasta <- system.file("extdata", "example_train_set.fa.gz", package="dada2")
#' taxa <- assignTaxonomy(seqs, training_fasta)
#' species_fasta <- system.file("extdata", "example_species_assignment.fa.gz", package="dada2")
#' taxa.spec <- addSpecies(taxa, species_fasta)
#' taxa.spec.multi <- addSpecies(taxa, species_fasta, allowMultiple=TRUE)
#' 
addSpecies <- function(taxtab, refFasta, allowMultiple=FALSE, tryRC=FALSE, n=2000, verbose=FALSE) {
  seqs <- rownames(taxtab)
  binom <- assignSpecies(seqs, refFasta=refFasta, allowMultiple=allowMultiple, tryRC=tryRC, n=n, verbose=verbose)
  # Merge tables
  if("Genus" %in% colnames(taxtab)) gcol <- which(colnames(taxtab) == "Genus")
  else gcol <- ncol(taxtab)
  # Match genera
  gen.match <- mapply(matchGenera, taxtab[,gcol], binom[,1])
  taxtab <- cbind(taxtab, binom[,2])
  colnames(taxtab)[ncol(taxtab)] <- "Species"
  taxtab[!gen.match,"Species"] <- NA_character_
  if(verbose) cat("Of which", sum(!is.na(taxtab[,"Species"])),"had genera consistent with the input table.")
  taxtab
}

#' This function creates the dada2 assignTaxonomy training fasta for the RDP trainset .fa file
#' 
#' ## RDP Trainset 16
#' path <- "~/Desktop/RDP/RDPClassifier_16S_trainsetNo16_rawtrainingdata"
#' dada2:::makeTaxonomyFasta_RDP(file.path(path, "trainset16_022016.fa"), 
#'     file.path(path, "trainset16_db_taxid.txt"), 
#'     "~/tax/rdp_train_set_16.fa.gz")
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead writeFasta
#' @importFrom ShortRead sread
#' @importFrom Biostrings BStringSet
#' @importFrom utils read.table
#' @keywords internal
makeTaxonomyFasta_RDP <- function(fin, fdb, fout, compress=TRUE) {
  # Read in the fasta and pull out the taxonomy entries
  sr <- readFasta(fin)
  id <- as.character(gsub("\"", "", id(sr)))
  tax <- sapply(strsplit(id, "\\t"), `[`, 2)
  tax <- gsub("^Root;", "", tax)
  tax <- strsplit(tax, ";")
  # Get the names of the standard 6 taxonomic levels
  db <- read.table(file=fdb, sep="*", stringsAsFactors=FALSE)
  colnames(db) <- c("Index", "Name", "L", "R", "Level")
  keep <- db$Name[db$Level %in% c("domain", "phylum", "class", "order", "family", "genus")]
  # Cut down to just the 6 level taxonomy
  tax <- sapply(tax, function(x) x[x %in% keep])
  tax <- lapply(tax, paste, collapse = ";")
  tax <- unlist(tax)
  # Final formatting
  tax <- paste0(tax, ";") # Ending semicolon
  tax <- gsub("[^;]*_incertae_sedis;$", "", tax) # Uncertain lowest-level assignment is better to leave blank
  tax <- gsub(" ", "_", tax)
  # Write to disk
  writeFasta(ShortRead(sread(sr), BStringSet(tax)), fout,
             width=20000L, compress=compress)
}

#' This function creates the dada2 assignSpecies fasta file for the RDP
#' from the RDP's _Bacteria_unaligned.fa file.
#' 
#' ## RDP Trainset 16/Release 11.5
#' dada2:::makeSpeciesFasta_RDP("~/Desktop/RDP/current_Bacteria_unaligned.fa", "~/tax/rdp_species_assignment_16.fa.gz")
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead writeFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead narrow
#' @importFrom IRanges narrow
#' @importFrom Biostrings BStringSet
#' @keywords internal
makeSpeciesFasta_RDP <- function(fin, fout, compress=TRUE) {
  # Read in and remove records not assigned to species
  sr <- readFasta(fin)
  is.uncult <- grepl("[Uu]ncultured", id(sr))
  sr <- sr[!is.uncult]
  is.unclass <- grepl("[Uu]nclassified", id(sr))
  sr <- sr[!is.unclass]
  is.outgroup <- (grepl("Outgroup", id(sr)))
  sr <- sr[!is.outgroup]
  is.unident <- grepl("[Uu]nidentified", id(sr))
  sr <- sr[!is.unident]
  
  # Pull out the genus species binomial string
  binom <- sapply(strsplit(as.character(id(sr)), ";"), `[`, 1)
  binom <- sapply(strsplit(binom, "\\t"), `[`, 1)
  binom <- gsub(" \\(T\\)", "", binom)
  binom <- gsub("\\[", "", binom)
  binom <- gsub("\\]", "", binom)
  
  # Match genera between binomial and the curated taxonomy
  bar <- strsplit(as.character(id(sr)), ";")
  barlens <- sapply(bar, length)
  geni <- mapply(function(x,y) x[[y]], bar, barlens-1)
  
  # Get rid of SXXX id
  binom <- gsub("^S[0123456789]{9} ", "", binom)
  binom <- gsub("\'" , "", binom)
  # Drop Candidatus strings
  binom <- gsub("Candidatus ", "", binom)
  geni <- gsub("Candidatus ", "", geni)
  
  # Subset down to those binomials which match the curated genus
  binom.geni <- sapply(strsplit(binom, "\\s"), `[`, 1)
  gen.match <- mapply(matchGenera, geni, binom.geni)
  sr <- sr[gen.match]
  binom <- binom[gen.match]
  geni <- geni[gen.match]
  
  # Make matrix of genus/species
  binom[sapply(strsplit(binom, "\\s"), length)==1] <- paste(binom[sapply(strsplit(binom, "\\s"), length)==1], "sp.")
  binom2 <- cbind(sapply(strsplit(binom, "\\s"), `[`, 1),
                  sapply(strsplit(binom, "\\s"), `[`, 2))
  # Keep only those with a named species
  has.spec <- !grepl("sp\\.", binom2[,2])
  sum(has.spec)
  binom2 <- binom2[has.spec,]
  sr <- sr[has.spec]
  binom <- binom[has.spec]
  geni <- geni[has.spec]
  cat(length(binom), "sequences with genus/species binomial annotation output.\n")
  
  # Write to disk
  ids <- as.character(narrow(id(sr),1,10))
  writeFasta(ShortRead(sread(sr), BStringSet(paste(ids, binom))), fout,
             width=20000L, compress=compress)
}

#' This function creates the dada2 assignTaxonomy training fasta for the Silva .align file
#' generated by the Mothur project.
#' 
#' ## Silva release v128
#' path <- "~/Desktop/Silva/Silva.nr_v128"
#' dada2:::makeTaxonomyFasta_Silva(file.path(path, "silva.nr_v128.align"), 
#'     file.path(path, "silva.nr_v128.tax"), 
#'     "~/tax/silva_nr_v128_train_set.fa.gz")
#' 
#' ## Silva release v132
#' path <- "~/Desktop/Silva/Silva.nr_v132"
#' dada2:::makeTaxonomyFasta_Silva(file.path(path, "silva.nr_v132.align"), 
#'     file.path(path, "silva.nr_v132.tax"), 
#'     "~/tax/silva_nr_v132_train_set.fa.gz")
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead writeFasta
#' @importFrom ShortRead sread
#' @importFrom Biostrings BStringSet
#' @importFrom utils read.table
#' @keywords internal
makeTaxonomyFasta_Silva <- function(fin, ftax, fout, compress=TRUE) {
  # Read in the fasta and pull out the taxonomy entries
  sr <- readFasta(fin) # ~10GB to read in
  ids <- sapply(strsplit(as.character(id(sr)), "\\t"), `[`, 1)
  seqs <- gsub("[.-]", "", sread(sr)) # Takes a while
  rm(sr);gc()
  # Read in the taxnoomy file
  taxdf <- read.table(ftax, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  colnames(taxdf) <- c("id", "tax")
  taxdf$tax <- gsub("^\\s+|\\s+$", "", taxdf$tax)
  if(!identical(taxdf$id, ids)) stop("Input align and taxonomy files don't match.")
  # Final formatting
  tax <- taxdf$tax
  tax <- gsub("Escherichia-Shigella", "Escherichia/Shigella", tax)
  # Remove faux-assignments added by new Mothur processing script
  tax <- gsub("[^;]*_ge;$", "", tax)
  tax <- gsub("[^;]*_fa;$", "", tax)
  tax <- gsub("[^;]*_or;$", "", tax)
  tax <- gsub("[^;]*_cl;$", "", tax)
  tax <- gsub("[^;]*_ph;$", "", tax)
  tax <- gsub(";uncultured;$", ";", tax)
  # Write to disk
  writeFasta(ShortRead(DNAStringSet(seqs), BStringSet(tax)), fout,
             width=20000L, compress=compress)
}

#' This function creates the dada2 assignSpecies fasta file for Silva
#' from the SILVA_[VERSION]_SSURef_tax_silva.fasta file
#' 
#' ## Silva release v128
#' dada2:::makeSpeciesFasta_Silva("~/Desktop/Silva/SILVA_128_SSURef_tax_silva.fasta.gz", 
#'     "~/tax/silva_species_assignment_v128.fa.gz")
#' 
#' ## Silva release v132
#' dada2:::makeSpeciesFasta_Silva("~/Desktop/Silva/SILVA_132_SSURef_tax_silva.fasta.gz", 
#'     "~/tax/silva_species_assignment_v132.fa.gz")
#' 
#' Output: 313502 sequences with genus/species binomial annotation output.
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead writeFasta
#' @importFrom ShortRead sread
#' @importFrom Biostrings BStringSet
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings readRNAStringSet
#' @keywords internal
makeSpeciesFasta_Silva <- function(fin, fout, compress=TRUE) {
  # Read in and remove records not assigned to species and non-bacteria
  xset <- DNAStringSet(readRNAStringSet(fin, format="fasta"))
  is.bact <- grepl("Bacteria;", names(xset), fixed=TRUE)
  xset <- xset[is.bact]
  is.uncult <- grepl("[Uu]ncultured", names(xset))
  xset <- xset[!is.uncult]
  is.unident <- grepl("[Uu]nidentified", names(xset))
  xset <- xset[!is.unident]
  is.complete <- sapply(strsplit(as.character(names(xset)), ";"), length)==7
  xset <- xset[is.complete]
  
  # Pull out binomial strings
  binom <- strsplit(as.character(names(xset)), ";")
  genus <- sapply(binom, `[`, 6)
  binom <- sapply(binom, `[`, 7)
  
  genus <- gsub("Candidatus ", "", genus)
  binom <- gsub("Candidatus ", "", binom)
  genus <- gsub("\\[", "", genus)
  genus <- gsub("\\]", "", genus)
  binom <- gsub("\\[", "", binom)
  binom <- gsub("\\]", "", binom)
  
  # Subset down to those binomials which match the curated genus
  genus.binom <- sapply(strsplit(binom, "\\s"), `[`, 1)
  gen.match <- mapply(matchGenera, genus, genus.binom, split.glyph="-")
  # Note that raw Silva files use Escherichia-Shigella, but this is changed to Escherichia/Shigella in dada2 version
  xset <- xset[gen.match]
  binom <- binom[gen.match]
  genus <- genus[gen.match]

  # Make matrix of genus/species
  binom[sapply(strsplit(binom, "\\s"), length)==1] <- paste(binom[sapply(strsplit(binom, "\\s"), length)==1], "sp.")
  binom2 <- cbind(sapply(strsplit(binom, "\\s"), `[`, 1),
                  sapply(strsplit(binom, "\\s"), `[`, 2))
  # Keep only those with a named species
  has.spec <- !grepl("sp\\.", binom2[,2]) & !(binom2[,2]=="endosymbiont")
  binom2 <- binom2[has.spec,]
  xset <- xset[has.spec]
  binom <- binom[has.spec]
  genus <- genus[has.spec]
  cat(length(binom), "sequences with genus/species binomial annotation output.\n")
  
  # Write to disk
  ids <- sapply(strsplit(as.character(names(xset)), "\\s"), `[`, 1)
  writeFasta(ShortRead(unname(xset), BStringSet(paste(ids, binom))), fout,
             width=20000L, compress=compress)
}
