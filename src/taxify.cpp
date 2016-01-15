#include "dada.h"
#include <Rcpp.h>
using namespace Rcpp;

unsigned int tax_kmer(const char *seq, unsigned int k) {
  unsigned int j, nti;
  unsigned int kmer=0;
  
  for(j=0; j<k; j++) {
    if(seq[j] == 'A') {
      nti = 0;
    } else if (seq[j] == 'C') {
      nti = 1;
    } else if (seq[j] == 'G') {
      nti = 2;
    } else if (seq[j] == 'T') {
      nti = 3;
    } else {
      Rcpp::stop("Unexpected nucleotide.");
    }
    kmer = 4*kmer + nti;
  }
  return(kmer);

}

void tax_kvec(const char *seq, unsigned int k, unsigned char *kvec) {  // Assumes a clean seq (just A/C/G/T)
  unsigned int i;
  unsigned int len = strlen(seq);
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  for(i=0;i<n_kmers;i++) { kvec[i] = 0; }

  for(i=0; i<len-k; i++) {
    kmer = tax_kmer(&seq[i], k);
    
    // Make sure kmer index is valid.
    if(kmer >= n_kmers) {
      Rcpp::stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer] = 1;
    }
  }
}

void tax_karray(const char *seq, unsigned int k, unsigned int *karray) {
  unsigned int i, kmer;
  unsigned int len = strlen(seq);
  
  for(i=0;i<len-k;i++) {
    kmer = tax_kmer(&seq[i], k);
    karray[i] = kmer;
  }
}

//------------------------------------------------------------------
// Assigns taxonomy to sequence based on provided ref seqs and corresponding taxonomies.
//
// [[Rcpp::export]]
Rcpp::CharacterVector C_taxify(std::vector<std::string> seqs, std::vector<std::string> refs, std::vector<int32_t> ref_to_genus, std::vector<std::string> taxs) {
  unsigned int i, j, pos, kmer, g;
  unsigned int k=8;
  unsigned int n_kmers = (1 << (2*k));
  unsigned int nseq = seqs.size();
  if(nseq == 0) Rcpp::stop("No seqs provided to classify.");
  unsigned int nref = refs.size();
  if(nref != ref_to_genus.size()) Rcpp::stop("Length mismatch between number of references and map to genus.");
  unsigned int ngenus = taxs.size();
  for(i=0;i<nref;i++) {
    if(ref_to_genus[i] <= 0 || ref_to_genus[i] > ngenus) Rcpp::stop("Invalid ref_to_genus map.");
    ref_to_genus[i]--; // R 1-index -> C 0-index
  } 
  // Make reference kvec array by sequence
  unsigned char *ref_kmers = (unsigned char *) malloc((nref * n_kmers) * sizeof(unsigned char));
  if(ref_kmers == NULL) Rcpp::stop("Memory allocation failed.");
  
  for(i=0;i<nref;i++) {
    tax_kvec(refs[i].c_str(), k, &ref_kmers[i*n_kmers]);
  }
  
  // Count seqs in each genus (M_g)
  uint16_t *genus_num = (uint16_t *) calloc(ngenus, sizeof(uint16_t));
  if(genus_num == NULL) Rcpp::stop("Memory allocation failed.");  
  for(i=0;i<nref;i++) {
    genus_num[ref_to_genus[i]]++;
  }
  
  // Make reference kvec array by genus
  uint16_t *genus_kmers = (uint16_t *) calloc((ngenus * n_kmers), sizeof(uint16_t));
  if(genus_kmers == NULL) Rcpp::stop("Memory allocation failed.");
  
  uint16_t *genus_kv;
  unsigned char *ref_kv;
  for(i=0;i<nref;i++) {
    g = ref_to_genus[i];
    ref_kv = &ref_kmers[i*n_kmers];
    genus_kv = &genus_kmers[g*n_kmers];
    for(kmer=0;kmer<n_kmers;kmer++) {
      if(ref_kv[kmer]) { genus_kv[kmer]++; }
    }
  }
  
  // Calculate word priors (Pi)
  unsigned int tot;
  double *kmer_prior = (double *) malloc(n_kmers * sizeof(double));
  if(kmer_prior == NULL) Rcpp::stop("Memory allocation failed.");
  for(kmer=0;kmer<n_kmers;kmer++) {
    tot=0;
    for(i=0;i<nref;i++) {
      tot += ref_kmers[i*n_kmers + kmer];
    }
    kmer_prior[kmer] = (tot + 0.5)/(1.0 + nref);
  }

  unsigned int arraylen = 0;
  unsigned int seqlen;
  for(i=0;i<nseq;i++) {
    seqlen = seqs[i].size();
    if(seqlen < 50) Rcpp::stop("Sequences must be at least 50 nts to classify.");
    if((seqlen-k) > arraylen) { arraylen = seqlen-k; }
  }

  // Allocate kmer array to be used by the seqs
  unsigned int *karray = (unsigned int *) malloc(arraylen * sizeof(unsigned int));
  if(karray == NULL) Rcpp::stop("Memory allocation failed.");

  Rcpp::CharacterVector rval;
  int max_g;
  double max_p;
  double p;
  for(j=0;j<nseq;j++) {
    seqlen = seqs[j].size();
    tax_karray(seqs[j].c_str(), k, karray);
    // Find best hit
    max_g = -1;
    max_p = 0.0;
    for(g=0;g<ngenus;g++) {
      p = 1.0;
      for(pos=0;pos<(seqlen-k);pos++) {
        kmer = karray[pos];
        p = p * ((genus_kmers[g*n_kmers + kmer] + kmer_prior[kmer])/(1.0 + genus_num[g]));
      }
      if(p>max_p) {
        max_p = p;
        max_g = g;
      }
    }
    if(max_g < 0) Rcpp::stop("Assignment failed.");
    rval.push_back(taxs[max_g]);
//    Rprintf("max_p = %.2e, max_g = %i.\n", max_p, max_g); // -> 1-index
  }
  
  free(ref_kmers);
  free(genus_num);
  free(genus_kmers);
  free(kmer_prior);
  free(karray);
  
//  return(Rcpp::CharacterVector::create(taxs[max_i]));
  return(rval);
}
