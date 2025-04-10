#'
#' Classifies sequences against reference training dataset.
#' 
#' assignTaxonomy implements the Naive Bayesian Classifier algorithm described in
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
                           multithread=FALSE, verbose=FALSE, nboot=100) {
  MIN_REF_LEN <- 20 # Enforced minimum length of reference seqs. Must be bigger than the kmer-size used (8).
  MIN_TAX_LEN <- 50 # Minimum length of input sequences to get a taxonomic assignment
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  if(min(nchar(seqs)) < MIN_TAX_LEN) {
    warning("Some sequences were shorter than ", MIN_TAX_LEN, " nts and will not receive a taxonomic classification.")
  }
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
  ### Assign
  # Parse multithreading argument
  if(is.logical(multithread)) {
    if(multithread==TRUE) { RcppParallel::setThreadOptions(numThreads = "auto") }
    else { RcppParallel::setThreadOptions(numThreads = 1) }
  } else if(is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
  } else {
    warning("Invalid multithread parameter. Running as a single thread.")
    RcppParallel::setThreadOptions(numThreads = 1)
  }
  # Run C assignemnt code
  assignment <- C_assign_taxonomy2(seqs, rc(seqs), refs, ref.to.genus, tax.mat.int, tryRC, verbose, nboot)
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
  if(nchar(gen.tax) == 0 || nchar(gen.binom) == 0) { return(FALSE) }
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
#' coercible by \code{\link{getUniques}}. Sequences must be A/C/G/T only.
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
  if(!all(C_isACGT(seqs))) stop("Non-ACGT characters present in the query sequences.")
  if(!length(unlist(strsplit(ids[[1]], "\\s"))) >= 3) {
    if(length(unlist(gregexpr(";", ids[[1]]))) >= 3) {
      stop("Incorrect reference file format for assignSpecies (this looks like a file formatted for assignTaxonomy).")
    } else {
      stop("Incorrect reference file format for assignSpecies.")
    }
  }
  genus <- sapply(strsplit(ids, "\\s"), `[`, 2)
  species <- sapply(strsplit(ids, "\\s"), `[`, 3)
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

#' This function creates the dada2 assignTaxonomy training fasta from the speciesrank RDP trainset .fa file
#' The RDP trainset data was downloaded from: https://sourceforge.net/projects/rdp-classifier/
#' 
#' ## RDP Trainset 19
#' path <- "~/tax/rdp/v19"
#' dada2:::makeTaxonomyFasta_RDP(file.path(path, "trainset19_072023_speciesrank.fa"), 
#'                               file.path(path, "trainset19_db_taxid.txt"), 
#'                               "~/Desktop/rdp_19_toGenus_trainset.fa.gz", include.species=FALSE,
#'                               compress=TRUE)
#' dada2:::tax.check("~/Desktop/rdp_19_toGenus_trainset.fa.gz")
#' 
#' dada2:::makeTaxonomyFasta_RDP(file.path(path, "trainset19_072023_speciesrank.fa"), 
#'                               file.path(path, "trainset19_db_taxid.txt"), 
#'                               "~/Desktop/rdp_19_toSpecies_trainset.fa.gz", include.species=TRUE,
#'                               compress=TRUE)
#' dada2:::tax.check("~/Desktop/rdp_19_toSpecies_trainset.fa.gz")
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead writeFasta
#' @importFrom ShortRead sread
#' @importFrom Biostrings BStringSet
#' @importFrom utils read.table
#' @keywords internal
makeTaxonomyFasta_RDP <- function(fin, fdb, fout, include.species=FALSE, compress=TRUE) {
  # Read in the fasta and pull out the taxonomy entries
  sr <- readFasta(fin)
  id <- as.character(id(sr)) # ID has 3 fields separated by tabs: Accession -- Species binomial+ -- Taxonomy to genus level
  tax <- sapply(strsplit(id, "\\t"), `[`, 3)
  tax <- gsub("[a-z]{5,8}__", "", tax)
  tax <- gsub("; ", ";", tax)
  tax <- strsplit(tax, ";")
  # Get the names of the standard 6 taxonomic levels
  db <- read.table(file=fdb, sep="*", stringsAsFactors=FALSE)
  colnames(db) <- c("Index", "Name", "L", "R", "Level")
  keep <- db$Name[db$Level %in% c("domain", "phylum", "class", "order", "family", "genus")]
  # Cut down to just the 6 level taxonomy
  tax <- sapply(tax, function(x) x[x %in% keep])
  if(max(sapply(tax, length)) > 6) stop("Taxonomy with >6 levels detected.")
  # Add handling of species binomial
  if(include.species) {
    binom <- sapply(strsplit(id, "\\t"), `[`, 2) # Format: Genus species additional info for strain and accession (space separators)
    gen.binom <- sapply(strsplit(binom, "\\s"), `[`, 1)
    spc.binom <- sapply(strsplit(binom, "\\s"), `[`, 2)
    gen <- sapply(tax, `[`, 6) # Will be NA if absent
    gen.match <- mapply(matchGenera, gen, gen.binom)
    spc.binom[!gen.match] <- NA
    # If not NA and tax has 6 levels, append the species
    nspc <- 0
    for(i in seq_along(tax)) {
      if(!is.na(spc.binom[[i]]) && length(tax[[i]] == 6)) {
        tax[[i]] <- c(tax[[i]], spc.binom[[i]])
        nspc <- nspc+1
      }
    }
  }
  tax <- lapply(tax, paste, collapse = ";")
  tax <- unlist(tax)
  # Final formatting
  tax <- paste0(tax, ";") # Ending semicolon
  tax <- gsub("[^;]*_incertae_sedis;$", "", tax) # Uncertain lowest-level assignment is better to leave blank
  tax <- gsub(" ", "_", tax)
  ## Add some verbose output describing what happened.
  cat(length(tax), "reference sequences were output.\n")
  if(include.species) cat(nspc, "had valid species names.\n")
  # Write to disk
  writeFasta(ShortRead(sread(sr), BStringSet(tax)), fout,
             width=20000L, compress=compress)
}

#' This function creates the dada2 assignSpecies fasta file for the RDP
#' from the RDP's _Bacteria_unaligned.fa file.
#' 
#' THE WAY RDP RELEASES SPECIES LEVEL INFORMATION APPEARS TO HAVE CHANGED IN RELEASE 19
#' AS A RESULT, THIS OUTPUT IS NOT CURRENTLY BEING MAINTAINED
#' 
#' ## RDP Trainset 18/Release 11.5
#' ## The RDP documentation does not make clear whether the updates to the taxonomy from training set release 18 were
#' ## propagated to the current Bacterial alignment.
#' dada2:::makeSpeciesFasta_RDP("~/Desktop/RDP/current_Bacteria_unaligned.fa", "~/tax/rdp_species_assignment_18.fa.gz")
#' dada2:::tax.check("~/tax/rdp_species_assignment_18.fa.gz", mode="species")
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

#' This function creates the dada2 assignTaxonomy training fasta for the official Silva NR99
#' release files. If `include.species`=TRUE, a 7th taxonomic level (species) will be added based on the
#' Genus species binomial in the Silva taxonomy string (if consistent with the genus assignment).
#' 
#' ## Silva release v138.2
#' path <- "~/tax/Silva/v138_2"
#' dada2:::makeTaxonomyFasta_SilvaNR(file.path(path, "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"), 
#'     file.path(path, "tax_slv_ssu_138.2.txt"), 
#'     "~/Desktop/silva_nr99_v138.2_toGenus_trainset.fa.gz")
#' dada2:::tax.check("~/Desktop/silva_nr99_v138.2_toGenus_trainset.fa.gz")
#'     
#' dada2:::makeTaxonomyFasta_SilvaNR(file.path(path, "SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"), 
#'     file.path(path, "tax_slv_ssu_138.2.txt"), 
#'     include.species=TRUE, "~/silva_nr99_v138.2_toSpecies_trainset.fa.gz")
#' dada2:::tax.check("~/Desktop/silva_nr99_v138.2_toSpecies_trainset.fa.gz")
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead writeFasta
#' @importFrom ShortRead sread
#' @importFrom Biostrings BStringSet
#' @importFrom utils read.table
#' @keywords internal
makeTaxonomyFasta_SilvaNR <- function(fin, ftax, fout, include.species=FALSE, compress=TRUE) {
  xset <- DNAStringSet(readRNAStringSet(fin, format="fasta"))
  # taxl: The taxonmic strings or (l)ines associated with each entry. Named by the sequence ID/accession.
  taxl <- names(xset)
  names(taxl) <- sapply(strsplit(names(xset), "\\s"), `[`, 1)
  if(any(duplicated(names(taxl)))) stop("Duplicated sequence IDs detected.")
  names(xset) <- names(taxl)
  taxl <- gsub("^[A-Za-z0-9.]+\\s", "", taxl)
  # taxa: A list of the ordered taxonomic levels corresponding to each reference sequence. Named by the sequence ID/accession.
  taxa <- strsplit(taxl, ";")
  # Read in the defined Silva taxonomic levels, e.g. Bacteria;Desulfobacterota;Desulfobulbia;Desulfobulbales;Desulfurivibrionaceae;
  silva.taxa <- read.table(ftax, sep="\t", col.names=c("Taxon", "V2", "Level", "V4", "V5"), stringsAsFactors=FALSE)
  silva.taxa <- silva.taxa[,c("Taxon", "Level")]
  # Subset down to Bacteria and Archaea
  kingdom <- sapply(strsplit(taxl, ";"), `[`, 1)
  taxl.ba <- taxl[kingdom %in% c("Bacteria", "Archaea")]
  taxa.ba <- taxa[names(taxl.ba)]
  # Create 6-column matrix with Silva taxonomic assignment for each sequence at each level from Kingdom to Genus (NA if no assignment)
  taxa.ba.mat <- matrix(sapply(taxa.ba, function(flds) {
    c(flds[1], flds[2], flds[3], flds[4], flds[5], flds[6])
  }), ncol=6, byrow=TRUE)
  rownames(taxa.ba.mat) <- names(taxl.ba)
  # Create 6-column matrix with full Silva taxonomic string at each level for each sequence, from Kingdom to Genus
  # Strings will include NA levels if no assignment at that level, e.g. Bacteria;Firmicutes;NA;NA
  taxa.ba.mat.string <- matrix("UNDEF", nrow=nrow(taxa.ba.mat), ncol=ncol(taxa.ba.mat))
  rownames(taxa.ba.mat.string) <- names(taxl.ba)
  taxa.ba.mat.string[,1] <- paste0(taxa.ba.mat[,1],";")
  for(col in seq(2,6)) {
    taxa.ba.mat.string[,col] <- paste0(taxa.ba.mat.string[,col-1], taxa.ba.mat[,col],";")
  }
  if(any(taxa.ba.mat.string == "UNDEF")) stop("Taxon string matrix was not fully initialized.")

  ##### SILVA 138_2 CLEANUP NO LONGER NEEDS THIS BLOCK, DEPRECATING IN COMMENTS FOR NOW, SHOULD REMOVE LATER
  # Define the set of valid taxonomic assignment by their appearance in the list of valid Silva taxonomic levels
  taxa.ba.mat.is_valid <- matrix(taxa.ba.mat.string %in% silva.taxa$Taxon, ncol=6)
  # Update taxa.ba.mat matrix by replacing invalid entries with NAs
  taxa.ba.mat[!taxa.ba.mat.is_valid] <- NA
  # Also replace "uncultured" taxonomic ranks with NAs (note, uncultured only shows up as the terminal "assigned" rank)
  taxa.ba.mat[taxa.ba.mat %in% c("Uncultured", "uncultured")] <- NA
  #####
  
  ##### PROCESSING FOR UPDATES USE OF INCERTAE SEDIS IN SILVA 138_2
  # Change terminal "Incertae Sedis" assignments to NA
  # Note: Non-terminal "Incertae Sedis" assignments will remain
  #  That is, when "Incertae Sedis" appears at an intermediate rank, but lower ranks are assigned a valid name
  # Example: AF282253.1.1503  "Bacteria" "Bacillota" "Incertae Sedis" "Thermolithobacterales" "Thermolithobacteraceae" "Thermolithobacter"
  taxa.ba.mat.is_incertae <- matrix(taxa.ba.mat %in% "Incertae Sedis", ncol=ncol(taxa.ba.mat))
  # Define make_na as TRUE when rank is "Incertae Sedis" and all lower ranks are also "Incertae Sedis"
  taxa.ba.mat.make_na <- taxa.ba.mat.is_incertae
  for(col in seq(6,2)) {
    taxa.ba.mat.make_na[,col-1] <- taxa.ba.mat.make_na[,col-1] & taxa.ba.mat.make_na[,col]
  }
  taxa.ba.mat[taxa.ba.mat.make_na] <- NA
  #####
  
  ######### ADD SPECIES PART HERE ##############
  if(include.species) {
    # Add the 7th column, which will be the species column
    taxa.ba.mat <- cbind(taxa.ba.mat, 
                         matrix(sapply(taxa.ba, `[`, 7), ncol=1, byrow=TRUE))
    # Get validated genus from the matrix
    genus <- taxa.ba.mat[,6]
    genus <- gsub("Candidatus ", "", genus)
    genus <- gsub("\\[", "", genus)
    genus <- gsub("\\]", "", genus)
    # Get the "binomial" string from the 7th field in the Silva taxonomic annotation
    # The "binomial" field is not curated like the other Silva taxonomic levels, and can have varying info
    # We assume that the first two words are the Genus species binomial, when there is a valid one in the field
    # NOTE: the binomial is actually not always in the 7th field, so this isn't strictly correct.
    # the binomial is in the "last" field, which may be <7 when not all the levels down to genus are assigned.
    # But we are throwing away everything that doesn't match the genus anyway, so that case
    # doesn't need to be handled correctly here.
    binom <- taxa.ba.mat[,7]
    binom <- gsub("Candidatus ", "", binom)
    binom <- gsub("\\[", "", binom)
    binom <- gsub("\\]", "", binom)
    # Pull out the first two fields, and turn binom into a two column matrix (Genus, species)
    binom <- cbind(sapply(strsplit(binom, "\\s"), `[`, 1),
                   sapply(strsplit(binom, "\\s"), `[`, 2))
    # Identify binomials that match the curated genus
    gen.match <- mapply(dada2:::matchGenera, genus, binom[,1], split.glyph="-")
    # Identify some other types of invalid species names
    is.NA <- apply(binom, 1, function(x) any(is.na(x)))
    is.sp <- grepl("sp\\.", binom[,2]) # "sp." is not a valid species name, just a generic
    is.endo <- binom[,1] %in% "endosymbiont" | binom[,2] %in% "endosymbiont"
    is.uncult <- grepl("[Uu]ncultured", binom[,1]) | grepl("[Uu]ncultured", binom[,2])
    is.unident <- grepl("[Uu]nidentified", binom[,1]) | grepl("[Uu]nidentified", binom[,2])
    # Define the "valid" species, and set invalid species to NA in the taxonomic matrix
    valid.spec <- gen.match & !is.NA & !is.sp & !is.endo & !is.uncult & !is.unident
    binom[!valid.spec,2] <- NA
    taxa.ba.mat[,7] <- binom[,2]
  }
  # Organize a small number of Eukaryota sequences for outgroup purposes, keeping only the Eukaryota Kingdom taxonomic assignment
  set.seed(500); N_EUK <- 500
  euk.keep <- sample(names(taxl)[kingdom %in% "Eukaryota"], N_EUK)
  taxa.euk.mat <- matrix("", nrow=N_EUK, ncol=ncol(taxa.ba.mat))
  rownames(taxa.euk.mat) <- euk.keep
  taxa.euk.mat[,1] <- "Eukaryota"
  taxa.euk.mat[,2:ncol(taxa.euk.mat)] <- NA
  
  # Now need to make the final training fasta in DADA2 format.
  taxa.mat.final <- rbind(taxa.ba.mat, taxa.euk.mat)
  taxa.string.final <- apply(taxa.mat.final, 1, function(x) {
    tst <- paste(x, collapse=";")
    tst <- paste0(tst, ";")
    tst <- gsub("NA;", "", tst)
    tst
  })
  
  if(any(is.na(names(taxa.string.final)))) stop("NA names in the final set of taxon strings.")
  if(!all(names(taxa.string.final) %in% names(xset))) stop("Some names of the final set of taxon strings don't match sequence names.")
  xset.out <- xset[names(taxa.string.final)]
  
  ## Add some verbose output describing what happened.
  cat(length(xset.out), "reference sequences were output.\n")
  print(table(taxa.mat.final[,1], useNA="ifany"))
  if(include.species) cat(sum(!is.na(taxa.mat.final[,7])), "entries include species names.\n")
  
  writeFasta(ShortRead(unname(xset.out), BStringSet(taxa.string.final)), fout,
             width=20000L, compress=compress)
}

#' This function creates the dada2 assignSpecies fasta file for Silva
#' from the SILVA_[VERSION]_SSURef_tax_silva.fasta file (NOT the NR99 file).
#' 
#' ## Silva release v138.2
#' dada2:::makeSpeciesFasta_Silva("~/tax/silva/v138_2/SILVA_138.2_SSURef_tax_silva.fasta.gz", 
#'     "~/Desktop/silva_v138.2_assignSpecies.fa.gz")
#' 
#' Output: 352047 sequences with genus/species binomial annotation output.
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
  tax <- strsplit(as.character(names(xset)), ";") ###!
  genus <- sapply(tax, `[`, 6) ###!
  binom <- sapply(tax, `[`, 7) ###!
  
  # Remove parens/brackets, which do not seem to be evenly used in the formal taxonomy and the binomial
  # ...and that mess up the `matchGenera` regex
  # Also fix the white-space issue with "Candidatus XXX" genus names
  genus <- gsub("Candidatus ", "Candidatus_", genus)
  binom <- gsub("Candidatus ", "Candidatus_", binom)
  genus <- gsub("\\[", "", genus)
  genus <- gsub("\\]", "", genus)
  binom <- gsub("\\[", "", binom)
  binom <- gsub("\\]", "", binom)
  genus <- gsub("\\(", "", genus)
  genus <- gsub("\\)", "", genus)
  binom <- gsub("\\(", "", binom)
  binom <- gsub("\\)", "", binom)
  
  # Subset down to those binomials which match the curated genus
  genus.binom <- sapply(strsplit(binom, "\\s"), `[`, 1)
  gen.match <- mapply(matchGenera, genus, genus.binom, split.glyph="-")
  # Note on split.glyph: raw Silva files use Escherichia-Shigella, but this is changed to Escherichia/Shigella in dada2 version
  xset <- xset[gen.match]
  binom <- binom[gen.match]
  genus <- genus[gen.match]
  
  # Make matrix of genus/species
  only.genus.in.binom <- sapply(strsplit(binom, "\\s"), length)==1
  binom[only.genus.in.binom] <- paste(binom[only.genus.in.binom], "sp.")
  binom2 <- cbind(sapply(strsplit(binom, "\\s"), `[`, 1),
                  sapply(strsplit(binom, "\\s"), `[`, 2))
  # Keep only those with a named species
  has.spec <- !grepl("sp\\.$", binom2[,2]) & !(binom2[,2]=="endosymbiont")
  binom2 <- binom2[has.spec,]
  xset <- xset[has.spec]
  binom <- binom[has.spec]
  genus <- genus[has.spec]
  cat(length(binom), "sequences with genus/species binomial annotation output.\n")
  
  # Write to disk
  ids <- sapply(strsplit(as.character(names(xset)), "\\s"), `[`, 1)
  writeFasta(ShortRead(unname(xset), BStringSet(paste(ids, binom2[,1], binom2[,2]))), fout,
             width=20000L, compress=compress)
}

#' This function creates the dada2 assignTaxonomy training fasta from the GreenGenes2
#' release files. If `include.species`=TRUE, the 7th taxonomic level (species) will be output.
#' Otherwise only the first 6 taxonomic levels (down to genus) will output.
#' 
#' ## Greengenes2 release 2024_09
#' path <- "~/tax/GG2/2024_09"
#' setwd(path)
#' # download.file("http://ftp.microbio.me/greengenes_release/current/2024.09.backbone.full-length.fna.qza", 
#'               "2024.09.backbone.full-length.fna.qza")
#' download.file("http://ftp.microbio.me/greengenes_release/current/2024.09.backbone.tax.qza",
#'               "2024.09.backbone.tax.qza")
#' unzip("2024.09.backbone.full-length.fna.qza")
#' unzip("2024.09.backbone.tax.qza")
#' fn <- "5b42d9b6-2f24-4f01-b989-9b4dafca7d5e/data/dna-sequences.fasta"
#' txfn <- "b7c3e691-ea51-4547-94dd-f79f49e41a36/data/taxonomy.tsv"
#' 
#' fn.out <- "~/Desktop/gg2_2024_09_toGenus_trainset.fa.gz"
#' dada2:::makeTaxonomyFasta_GG2(fn, txfn, fn.out, include.species=FALSE, compress=TRUE)
#' dada2:::tax.check(fn.out)
#' 
#' fn.out.spc <- "~/Desktop/gg2_2024_09_toSpecies_trainset.fa.gz"
#' dada2:::makeTaxonomyFasta_GG2(fn, txfn, fn.out.spc, include.species=TRUE, compress=TRUE)
#' dada2:::tax.check(fn.out.spc)
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead writeFasta
#' @importFrom utils read.csv
#' @keywords internal
makeTaxonomyFasta_GG2 <- function(fn, txfn, fout, include.species=FALSE, output.binomials = FALSE, compress=TRUE) {
  # tx is 2 column, id and taxonomy, taxonomy is 7 level, domain -- species
  sq <- getSequences(fn)
  tdf <- read.csv(txfn, sep="\t", header=TRUE)
  tax <- tdf[,2]
  names(tax) <- tdf[,1]
  sq <- sq[names(tax)]
  if(!identical(names(sq), names(tax))) stop("Mismatch between reference sequences and taxonomy file.")  ###!
  
  # Parse the taxonomies from the id string
  taxes <- strsplit(tax, "; ")
  tax.depth <- sapply(taxes, length)
  table(tax.depth) 
  # All taxonomies are 7-level
  # Note: A significant number of unassigned taxonomic levels here, which are encoded as e.g. "g__"
  # Note: Species has genus name duplicated, i.e. species level has genus [SPACE] species binomial, rather than just species.
  # Note: Enforcing consistency between the binomial genus in order to keep the species assignment.
  genus <- gsub("^g__", "", sapply(taxes, `[`, 6))
  genus <- gsub("Escherichia", "Escherichia_Shigella", genus)
  binom <- gsub("^s__", "", sapply(taxes, `[`, 7))
  binom.depth <- sapply(strsplit(binom, " "), length)
###  table(binom.depth) # All either 0 (s__) or 2
  has.binom <- binom.depth==2
  genus.binom <- sapply(strsplit(binom, "\\s"), `[`, 1)
  spec.binom <- sapply(strsplit(binom, "\\s"), `[`, 2)
  gen.match <- mapply(matchGenera, genus, genus.binom, split.glyph="_")
  # Replace the species field with just the species name, but only for concordant binomials
  for(i in seq_along(taxes)) {
    if(has.binom[i]) {
      if(gen.match[i]) { # Define the formatted species field
        if(output.binomials) { # Keep the GG2 species binomials in the species rank (but replace spaces).
          taxes[[i]][[7]] <- taxes[[i]][[7]] <- gsub(" ", "_", taxes[[i]][[7]])
        } else { # DEFAULT, just the species name
          taxes[[i]][[7]] <- paste0("s__", spec.binom[[i]]) 
        }
      } else { # !gen.match: scrub the species assignment
        taxes[[i]][[7]] <- "s__"
      }
    }
  }
  if(include.species) {
    cat(sum(has.binom), "out of", length(taxes), "sequences had a binomial species name assigned.\n", 
        sum(has.binom & !gen.match), "species assignments were removed as discordant with the genus assignment.")
  }
  
  # Identify unassigned taxonomic levels
  tax_pre <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  is.unassigned <- sapply(taxes, function(tx) {
    tx == tax_pre
  }) |> t()
  # Note: A relatively small number of entries have lower taxonomic levels assigned, even though higher taxonomic levels aren't assigned
  # e.g. `MJ030-2-barcode58-umi83452bins-ubs-6`, which is assigned s__Spirochaeta aurantia, but "o__;f__;g__" for order/family/genus designation
  # These assignments will be dropped in accordance with `assignTaxonmy` expectations
  tax.depth <- apply(is.unassigned, 1, function(isu) { min(which(isu)-1L, 7L) })
  tax.ids <- sapply(seq_along(taxes), function(i) {
    td <- tax.depth[[i]]
    if(!include.species) { td <- min(td, 6L) }
    id.str <- paste(taxes[[i]][1:td], collapse=";")
    id.str <- paste0(id.str, ";") # Add terminal semicolon
    id.str
  })
  names(tax.ids) <- names(taxes)
  # Write out the training fasta file
  sq.out <- sq
  names(sq.out) <- tax.ids
  writeFasta(sq.out, fout, width=20000L, compress=compress)
  ## Add some verbose output describing what happened.
  cat(length(sq.out), "reference sequences were output.\n")
}

## This uses the "ten_16s.100.fa" originally from Robert Edgar's taxonomy testing page: https://drive5.com/taxxi/doc/fasta_index.html
## This file is relicensed here under the DADA2 LGPL2 license on permission from Robert Edgar.
## Test file only contains taxonomy assigned to genus level (level=6), no species information
tax.check <- function(fn.tax, fn.test=system.file("extdata", "ten_16s.100.fa.gz", package="dada2"), nseq=100, level=6, mode="taxonomy") {
  sq.test <- sample(getSequences(fn.test), nseq)
  if(mode == "taxonomy") {
    tax <- assignTaxonomy(sq.test, fn.tax, multi=TRUE)
    rval <- cbind(unname(tax[,level]), sapply(strsplit(names(sq.test), ":"), `[`, level+1))
  } else if (mode=="species") {
    sq.acgt <- sq.test[dada2:::C_isACGT(sq.test)]
    spc <- assignSpecies(sq.acgt, fn.tax)
    rval <- cbind(unname(spc[,level-5]), sapply(strsplit(names(sq.acgt), ":"), `[`, level+1))
  } else { stop("Valid modes are taxonomy or species.") }
  colnames(rval) <- c("assigned", "reference")
  rval
}

