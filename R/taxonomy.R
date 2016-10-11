#'
#' Classifies sequences against reference training dataset.
#' 
#' assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm described in
#' Wang et al. Applied and Environmental Microbiology 2007, with kmer size 8 and 100 bootstrap
#' replicates.
#' 
#' @param seqs (Required). A character vector of the sequences to be assigned, or an object 
#' coercible by \code{\link{getUniques}}.
#'   
#' @param refFasta (Required). The path to the reference fasta file, or an 
#' R connection Can be compresssed.
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
#' @param taxLevels (Optional). Default is c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species").
#' The taxonomic levels being assigned. Truncates if deeper levels not present in
#' training fasta.
#'   
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, print status to standard output.
#' 
#' @return A character matrix of assigned taxonomies exceeding the minBoot level of
#'   bootstrapping confidence. Rows correspond to the provided sequences, columns to the
#'   taxonomic levels. NA indicates that the sequence was not consistently classified at
#'   that level at the minBoot threshhold. 
#' 
#' @export
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead id
#' 
#' @examples
#' \dontrun{
#'  taxa <- assignTaxonomy(dadaF, "gg_13_8_train_set_97.fa.gz")
#'  taxa <- assignTaxonomy(dadaF, "rdp_train_set_14.fa.gz", minBoot=80)
#' }
#' 
assignTaxonomy <- function(seqs, refFasta, minBoot=50,
                           taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                           verbose=FALSE, outputBootstraps=FALSE) {
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  # Read in the reference fasta
  refsr <- readFasta(refFasta)
  refs <- as.character(sread(refsr))
  tax <- as.character(id(refsr))
  tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x)) # Remove leading/trailing whitespace
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
  genus.unq
  ref.to.genus <- match(tax, genus.unq)
  tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol=td, byrow=TRUE)
  tax.df <- as.data.frame(tax.mat)
  for(i in seq(ncol(tax.df))) {
    tax.df[,i] <- factor(tax.df[,i])
    tax.df[,i] <- as.integer(tax.df[,i])
  }
  tax.mat.int <- as.matrix(tax.df)
  # Assign  
  assignment <- C_assign_taxonomy(seqs, refs, ref.to.genus, tax.mat.int, verbose)
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
      list(tax.out, boots.out)
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
#' R connection. Can be compresssed.
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
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, print status to standard output.
#' 
#' @return A two-column character matrix. Rows correspond to the provided sequences,
#'   columns to the genus and species taxonomic levels. NA indicates that the sequence
#'   was not classified at that level. 
#' 
#' @export
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead id
#' 
#' @examples
#' \dontrun{
#'  taxa <- assignSpecies(dadaF, "rdp_species.fa.gz")
#' }
#' 
assignSpecies <- function(seqs, refFasta, allowMultiple=FALSE, verbose=FALSE) {
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  # Read in the reference fasta
  refsr <- readFasta(refFasta)
  refs <- as(sread(refsr), "character")
  ids <- as(id(refsr), "character")
  genus <- sapply(strsplit(ids, "\\s"), `[`, 2)
  species <- sapply(strsplit(ids, "\\s"), `[`, 3)
  # Identify the exact hits
  hits <- lapply(seqs, function(x) grepl(x, refs, fixed=TRUE))
  # Define number of multiple species to return
  if(is.logical(allowMultiple)) {
    if(allowMultiple) keep <- Inf
    else keep <- 1
  } else {
    keep <- as.integer(allowMultiple)
  }
  # Get genus species return strings
  rval <- cbind(unlist(sapply(hits, mapHits, refs=genus, keep=1)),
                unlist(sapply(hits, mapHits, refs=species, keep=keep)))
  colnames(rval) <- c("Genus", "Species")
  rownames(rval) <- seqs
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
#' R connection. Can be compresssed.
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
#' \dontrun{
#'  taxtab <- assignTaxonomy(dadaF, "rdp_train_set_14.fa.gz")
#'  taxtab <- addSpecies(taxtab, "rdp_species_assignment_14.fa.gz")
#' }
#' 
addSpecies <- function(taxtab, refFasta, allowMultiple=FALSE, verbose=FALSE) {
  seqs <- rownames(taxtab)
  binom <- assignSpecies(seqs, refFasta=refFasta, allowMultiple=allowMultiple, verbose=verbose)
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

#' This function creates the dada2 assignSpecies fasta file for the RDP
#' from the rdp_Bacteria_unaligned.fa file.
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

#' This function creates the dada2 assignSpecies fasta file for Silva
#' from the SILVA_123.1_SSURef_3_3_2016_fix.fasta file.
#' @keywords internal
makeSpeciesFasta_Silva <- function(fin, fout, compress=TRUE) {
  # Read in and remove records not assigned to species and non-bacteria
  sr <- readFasta(fin)
  is.bact <- grepl("Bacteria;", id(sr), fixed=TRUE)
  sr <- sr[is.bact]
  is.uncult <- grepl("[Uu]ncultured", id(sr))
  sr <- sr[!is.uncult]
  is.unident <- grepl("[Uu]nidentified", id(sr))
  sr <- sr[!is.unident]
  is.complete <- sapply(strsplit(as.character(id(sr)), ";"), length)==7
  sr <- sr[is.complete]
  
  # Pull out binomial strings
  binom <- strsplit(as.character(id(sr)), ";")
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
  # Note that raw Silva files use Escherichia-Shigella, but this is changes to E/S in dada2 version
  sr <- sr[gen.match]
  binom <- binom[gen.match]
  genus <- genus[gen.match]

  # Make matrix of genus/species
  binom[sapply(strsplit(binom, "\\s"), length)==1] <- paste(binom[sapply(strsplit(binom, "\\s"), length)==1], "sp.")
  binom2 <- cbind(sapply(strsplit(binom, "\\s"), `[`, 1),
                  sapply(strsplit(binom, "\\s"), `[`, 2))
  # Keep only those with a named species
  has.spec <- !grepl("sp\\.", binom2[,2]) & !(binom2[,2]=="endosymbiont")
  binom2 <- binom2[has.spec,]
  sr <- sr[has.spec]
  binom <- binom[has.spec]
  genus <- genus[has.spec]
  cat(length(binom), "sequences with genus/species binomial annotation output.\n")
  
  # Write to disk
  ids <- sapply(strsplit(as.character(id(sr)), "\\s"), `[`, 1)
  writeFasta(ShortRead(sread(sr), BStringSet(paste(ids, binom))), fout,
             width=20000L, compress=compress)
}
