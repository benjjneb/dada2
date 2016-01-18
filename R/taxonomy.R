#' assignTaxonomy returns the taxonomic assignment of a set of sequences according to the provided reference file.
#' 
#' assignTaxonomy implements the RDP classifier algorithm in Wang 2007 with kmer size 8.
#' 
#' @param seqs (Required). A character vector of the sequences to be assigned.
#'   
#' @param refFasta (Required). A character string naming the path to the reference fasta file, or an 
#' R connection. This reference fasta file should be formatted with the id line corresponding to the
#' taxonomy (or classification) of the reference sequence, with each level separated by a semicolon. Eg.
#' 
#'  >Kingom;Phylum;Class;Order;Family;Genus;
#'  ACGAATGTGAAGTAA......
#' 
#' @param minBoot (Optional). Default 80. The minimum bootstrap confidence for assigning a taxonomic level.
#'   
#' @return \code{character}. An character vector of assigned taxonomies exceeding the minBoot level of
#'   bootstrapping confidence. Levels are separated by semicolons.
#' 
#' @export
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead id
#' 
assignTaxonomy <- function(seqs, refFasta, minBoot=80) {
  refsr <- readFasta(refFasta)
  refs <- as.character(sread(refsr))
  tax <- as.character(id(refsr))
  
  tax.depth <- sapply(strsplit(tax, ";"), length)
  td <- tax.depth[[1]]
  if(!all(sapply(strsplit(tax, ";"), length) == td)) {
    stop("References must all be assigned to the same taxonomic depth.")
  }

  tax.mat <- matrix(unlist(strsplit(tax, ";")), ncol=td, byrow=TRUE)
  tax.df <- as.data.frame(tax.mat)
  for(i in seq(ncol(tax.df))) {
    tax.df[,i] <- factor(tax.df[,i])
    tax.df[,i] <- as.integer(tax.df[,i])
  }
  tax.mat.int <- as.matrix(tax.df)
  tax.mat.int.aug <- cbind(tax.mat.int, as.integer(factor(tax)))
  ind_to_tax <- levels(factor(tax))
  ntype <- length(unique(tax))
  
  assignment <- C_assign_taxonomy(seqs, refs, tax.mat.int.aug, ntype)
  
  bestHit <- ind_to_tax[assignment$tax]
  boots <- assignment$boot
  taxes <- strsplit(bestHit, ";")
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,]>minBoot])
  taxes <- unlist(lapply(taxes, paste, collapse=";"))
  taxes <- paste0(taxes, ";") # Add suffix ;
  
  taxes
}