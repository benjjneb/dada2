#'
#' Classifies sequences against reference training set.
#' 
#' assignTaxonomy implements the RDP classifier algorithm in Wang 2007 with kmer size 8 and
#'  100 bootstrap replicates.
#' 
#' @param seqs (Required). A character vector of the sequences to be assigned, or an object coercible by
#' \code{\link{getUniques}}.
#'   
#' @param refFasta (Required). A character string naming the path to the reference fasta file, or an 
#' R connection. This reference fasta file should be formatted with the id line corresponding to the
#' taxonomy (or classification) of the reference sequence, with each level separated by a semicolon. Eg.
#' 
#'  >Kingom;Phylum;Class;Order;Family;Genus;
#'  ACGAATGTGAAGTAA......
#' 
#' @param minBoot (Optional). Default 50. The minimum bootstrap confidence for assigning a taxonomic level.
#'   
#' @param verbose (Optional). \code{logical(1)} indicating verbose text output. Default FALSE.
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
#'  taxa <- assignTaxonomy(dadaF, "gg_13_8_97_dada.fa")
#' }
#' 
assignTaxonomy <- function(seqs, refFasta, minBoot=50, verbose=FALSE) {
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
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,]>minBoot])
  # Convert to character matrix
  tax.out <- matrix(NA_character_, nrow=length(seqs), ncol=td)
  for(i in seq(length(seqs))) {
    tax.out[i,1:length(taxes[[i]])] <- taxes[[i]]
  }
  rownames(tax.out) <- seqs
  tax.out[tax.out=="_DADA2_UNSPECIFIED"] <- NA_character_
  tax.out
}