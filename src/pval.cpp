/*
Pval.cpp contains the functions related to calculating the abundance pval in DADA2.
*/

#include <Rcpp.h>
#include "dada.h"

// [[Rcpp::interfaces(r, cpp)]]

// Calculate abundance pval for given reads and expected number of reads
double calc_pA(int reads, double E_reads) {
  double norm, pval=1.;
  
  // Calculate norm (since conditioning on sequence being present).
  norm = (1.0 - exp(-E_reads));
  if(norm < TAIL_APPROX_CUTOFF) {
    norm = E_reads - 0.5*E_reads*E_reads; 
    // Assumption: TAIL_APPROX_CUTOFF is small enough to terminate taylor expansion at 2nd order
  }
  
  // Calculate pval from poisson cdf.
  Rcpp::IntegerVector n_repeats(1);
  n_repeats(0) = reads-1;
  Rcpp::NumericVector res = Rcpp::ppois(n_repeats, E_reads, false);  // lower.tail = false
  pval = Rcpp::as<double>(res);
  
  pval = pval/norm;
  return pval;
}

// Find abundance pval from a Raw in a Bi
double get_pA(Raw *raw, Bi *bi) {
  unsigned int hamming;
  double lambda, E_reads, pval = 1.;
  
  unsigned int ci = bi->comp_index[raw->index];
  lambda = bi->comp[ci].lambda;
  hamming = bi->comp[ci].hamming;
  
  if(raw->reads == 1) {   // Singleton. No abundance pval.
    pval=1.;
  } 
  else if(hamming == 0) { // Cluster center (or no mismatch to center)
    pval=1.;
  }
  else if(lambda == 0) { // Zero expected reads of this raw
    pval = 0.;
  }
  else { // Calculate abundance pval.
    // E_reads is the expected number of reads for this raw
    E_reads = lambda * bi->reads;
    pval = calc_pA(raw->reads, E_reads);
  }
  return pval;
}

// This calculates lambda from a lookup table index by transition (row) and rounded quality (col)
double compute_lambda(Raw *raw, Sub *sub, Rcpp::NumericMatrix errMat, bool use_quals) {
  // use_quals does nothing in this function, just here for backwards compatability for now
  int s, pos0, pos1, nti0, nti1, len1, ncol;
  double lambda;
  float prefactor, fqmin;
  int tvec[SEQLEN];
  int qind[SEQLEN];
  
  if(!sub) { // NULL Sub, outside Kmer threshold
    return 0.0;
  }
  
  // Make vector that indexes as integers the transitions at each position in seq1
  // Index is 0: exclude, 1: A->A, 2: A->C, ..., 5: C->A, ...
  len1 = raw->length;
  ncol = errMat.ncol();
  prefactor = ((float) (ncol-1))/((float) QMAX-QMIN);
  fqmin = (float) QMIN;
  for(pos1=0;pos1<len1;pos1++) {
    nti1 = ((int) raw->seq[pos1]) - 1;
    if(nti1 == 0 || nti1 == 1 || nti1 == 2 || nti1 == 3) {
      tvec[pos1] = nti1*4 + nti1;
    } else {
      Rcpp::stop("Error: Can't handle non ACGT sequences in CL3.");
    }
    if(raw->qual) {
      // Turn quality into the index in the array
      qind[pos1] = round(prefactor * (raw->qual[pos1] - fqmin));
    } else {
      qind[pos1] = 0;
    }
    
    if( qind[pos1] > (ncol-1) ) {
      Rcpp::stop("Error: rounded quality exceeded err lookup table.");
    }
  }

  // Now fix the ones where subs occurred
  for(s=0;s<sub->nsubs;s++) {
    pos0 = sub->pos[s];
    if(pos0 < 0 || pos0 >= len1) { Rcpp::stop("CL: Bad pos0."); }
    pos1 = sub->map[sub->pos[s]];
    if(pos1 < 0 || pos1 >= len1) { Rcpp::stop("CL: Bad pos1."); }
    
    nti0 = ((int) sub->nt0[s]) - 1;
    nti1 = ((int) sub->nt1[s]) - 1;
    tvec[pos1] = nti0*4 + nti1;
  }

  // And calculate lambda
  lambda = 1.0;
  for(pos1=0;pos1<len1;pos1++) {
    lambda = lambda * errMat(tvec[pos1], qind[pos1]);
  }

  if(lambda < 0 || lambda > 1) { Rcpp::stop("Bad lambda."); }

  return lambda;
}

