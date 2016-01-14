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
  size_t n_kmers = (2 << (2*k));  // 4^k kmers
  for(i=0;i<n_kmers;i++) { kvec[i] = 0; }

  for(i=0; i<len-k; i++) {
    kmer = tax_kmer(&seq[i], k);
    
    // Make sure kmer index is valid.
    if(kmer >= n_kmers) {
      Rcpp::stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
}

void tax_karray(const char *seq, unsigned int k, unsigned int *karray) {
  unsigned int i;
  unsigned int len = strlen(seq);
  
  for(i=0;i<len-k;i++) {
    karray[i] = tax_kmer(&seq[i], k);
  }
}

//------------------------------------------------------------------
// Assigns taxonomy to sequence based on provided ref seqs and corresponding taxonomies.
//
// [[Rcpp::export]]
Rcpp::CharacterVector C_taxify(std::string seq, std::vector<std::string> refs, std::vector<std::string> taxs) {
  unsigned int i, j;
  unsigned int k=8;
  unsigned int n_kmers = (2 << (2*k));
  unsigned int nref = refs.size();
  if(nref != taxs.size()) Rcpp::stop("Length mismatch between number of reference sequences and taxonomies.");
  unsigned int seqlen = seq.size();
  if(seqlen < 50) Rcpp::stop("Sequences must be at least 50 nts to classify.");
  unsigned int arraylen = seqlen-k;
  
  // Make reference kvec array
  unsigned char *ref_kmers = (unsigned char *) malloc((nref * n_kmers) * sizeof(unsigned char));
  if(ref_kmers == NULL) Rcpp::stop("Memory allocation failed.");
  
  for(i=0;i<nref;i++) {
    tax_kvec(refs[i].c_str(), k, &ref_kmers[i*n_kmers]);
  }

  // Make kmer list for seq
  unsigned int *karray = (unsigned int *) calloc(arraylen, sizeof(unsigned int));
  if(karray == NULL) Rcpp::stop("Memory allocation failed.");
  tax_karray(seq.c_str(), k, karray);
  
  int max_i = -1;
  unsigned int max_p = 0;
  unsigned int p;
  for(i=0;i<nref;i++) {
    p = 0;
    for(j=0;j<arraylen;j++) {
      p += ref_kmers[i*n_kmers + karray[j]];
    }
    if(p>max_p) {
      max_p = p;
      max_i = i;
    }
  }
  
  if(max_i < 0) Rcpp::stop("Assignment failed.");
  Rprintf("max_p = %i, max_i = %i.\n", max_p, max_i);
  
//  return(Rcpp::CharacterVector::create(taxs[max_i]));
  Rcpp::CharacterVector rval;
  rval.push_back(taxs[max_i]);
  return(rval);
}
