/*
Pval.cpp contains the functions related to calculating the abundance pval in DADA2.
*/

#include <Rcpp.h>
#include "dada.h"

// [[Rcpp::interfaces(r, cpp)]]

/* b_p_update:
 Calculates the abundance p-value for each raw in the clustering.
 Depends on the lambda between the raw and its cluster, and the reads of each.
 */
void b_p_update(B *b, bool greedy) {
  unsigned int i, r;
  Raw *raw;
  Bi *bi;
  double E_reads_center;
  for(i=0;i<b->nclust;i++) {
    bi = b->bi[i];
    if(bi->update_e) {
      for(r=0;r<bi->nraw;r++) {
        raw = bi->raw[r];
        raw->p = get_pA(raw, bi);
      } // for(r=0;r<b->bi[i]->nraw;r++)
      bi->update_e = false;
    } // if(bi->update_e)
    
    if(greedy && bi->check_locks) { // Lock raw if its reads are less than the expected reads from just the center 
      for(r=0;r<bi->nraw;r++) {
        raw = bi->raw[r];
        E_reads_center = b->bi[i]->center->reads * raw->comp.lambda;
        if(E_reads_center > raw->reads) { raw->lock = true; }
///!          if(raw->lock) { Rprintf("Raw %i (%i reads): E_Center = %.2f, P = %.2f\n", raw->index, raw->reads, E_reads_center, raw->p); }
        if(raw == bi->center) { raw->lock = true; }
      }
      bi->check_locks = false; // Locking can only happen first time around
    } // if(greedy && bi->check_locks)
  } // for(i=0;i<b->nclust;i++)
}

// Calculate abundance pval for given reads and expected number of reads
// Pval is conditional on sequnce being present, unless prior evidence is true
double calc_pA(int reads, double E_reads, bool prior) {
  double norm, pval=1.;
  
  // Calculate pval from poisson cdf.
  Rcpp::IntegerVector n_repeats(1);
  n_repeats(0) = reads-1; // -1 since strict > being calculated, and want to include the observed count
  Rcpp::NumericVector res = Rcpp::ppois(n_repeats, E_reads, false);  // lower.tail = false: P(X > x)
  pval = Rcpp::as<double>(res);

  if(!prior) {
    // Calculate norm (since conditioning on sequence being present).
    norm = (1.0 - exp(-E_reads));
    if(norm < TAIL_APPROX_CUTOFF) {
      norm = E_reads - 0.5*E_reads*E_reads; 
      // Assumption: TAIL_APPROX_CUTOFF is small enough to terminate taylor expansion at 2nd order
    }
    pval = pval/norm;
  }
  
  return pval;
}

// Find abundance pval from a Raw in a Bi
double get_pA(Raw *raw, Bi *bi) {
  unsigned int hamming;
  double lambda, E_reads, pval = 1.;
  
  lambda = raw->comp.lambda;
  hamming = raw->comp.hamming;

  if(raw->reads == 1 && !raw->prior) {   // Singleton. No abundance pval.
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
    pval = calc_pA(raw->reads, E_reads, raw->prior);
  }
  return pval;
}

// This calculates lambda from a lookup table index by transition (row) and rounded quality (col)
double compute_lambda(Raw *raw, Sub *sub, Rcpp::NumericMatrix errMat, bool use_quals, unsigned int ncol) {
  int s, pos0, pos1, nti0, nti1, len1;
  double lambda;
  int tvec[SEQLEN];
  unsigned int qind[SEQLEN];
  
  if(!sub) { // NULL Sub, outside Kmer threshold
    return 0.0;
  }
  
  // Make vector that indexes as integers the transitions at each position in seq1
  // Index is 0: A->A, 1: A->C, ..., 4: C->A, ...
  len1 = raw->length;
  for(pos1=0;pos1<len1;pos1++) {
    nti1 = ((int) raw->seq[pos1]) - 1;
    if(nti1 == 0 || nti1 == 1 || nti1 == 2 || nti1 == 3) {
      tvec[pos1] = nti1*4 + nti1;
    } else {
      Rcpp::stop("Non-ACGT sequences in compute_lambda.");
    }
    if(use_quals) {
      // Turn quality into the index in the array
      qind[pos1] = raw->qual[pos1]; // unsigned int
    } else {
      qind[pos1] = 0; // unsigned int
    }
  }
  
  // Now fix the ones where subs occurred
  for(s=0;s<sub->nsubs;s++) {
    pos0 = sub->pos[s];
    if(pos0 < 0 || pos0 >= sub->len0) { Rcpp::stop("CL: Bad pos0: %i (len0=%i).", pos0, sub->len0); }
    pos1 = sub->map[sub->pos[s]];
    if(pos1 < 0 || pos1 >= len1) { Rcpp::stop("CL: Bad pos1: %i (len1=%i).", pos1, len1); }
    
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

// This calculates lambda from a lookup table index by transition (row) and rounded quality (col)
double compute_lambda_ts(Raw *raw, Sub *sub, unsigned int ncol, double *err_mat, bool use_quals) {
  int s, pos0, pos1, nti0, nti1, len1;
  double lambda;
  unsigned int tvec[SEQLEN];
  unsigned int qind[SEQLEN];
  
  if(!sub) { // NULL Sub, outside Kmer threshold
    return 0.0;
  }
  
  // Make vector that indexes as integers the transitions at each position in seq1
  // Index is 0: A->A, 1: A->C, ..., 4: C->A, ...
  len1 = raw->length;
  for(pos1=0;pos1<len1;pos1++) {
    nti1 = ((int) raw->seq[pos1]) - 1;
    if(nti1 == 0 || nti1 == 1 || nti1 == 2 || nti1 == 3) {
      tvec[pos1] = nti1*4 + nti1;
    } else {
      Rcpp::stop("Non-ACGT sequences in compute_lambda.");
    }
    if(use_quals) {
      // Turn quality into the index in the array
      qind[pos1] = raw->qual[pos1]; // qind = unsigned int
    } else {
      qind[pos1] = 0;
    }
    
    if( qind[pos1] > (ncol-1) ) {
      Rcpp::stop("Rounded quality exceeded range of err lookup table.");
    }
  }
  
  // Now fix the ones where subs occurred
  for(s=0;s<sub->nsubs;s++) {
    pos0 = sub->pos[s];
    if(pos0 < 0 || pos0 >= sub->len0) { Rcpp::stop("CL: Bad pos0: %i (len0=%i).", pos0, sub->len0); }
    pos1 = sub->map[sub->pos[s]];
    if(pos1 < 0 || pos1 >= len1) { Rcpp::stop("CL: Bad pos1: %i (len1=%i).", pos1, len1); }
    
    nti0 = ((int) sub->nt0[s]) - 1;
    nti1 = ((int) sub->nt1[s]) - 1;
    tvec[pos1] = nti0*4 + nti1;
  }
  
  // And calculate lambda
  lambda = 1.0;
  for(pos1=0;pos1<len1;pos1++) {
    lambda = lambda * err_mat[tvec[pos1]*ncol+qind[pos1]];
  }
  
  if(lambda < 0 || lambda > 1) { Rcpp::stop("Bad lambda."); }
  
  return lambda;
}

/*
 *  Code below is modified from the R source code...
 *  https://github.com/wch/r-source/blob/af7f52f70101960861e5d995d3a4bec010bc89e6/src/nmath/ppois.c
 *  https://github.com/wch/r-source/blob/af7f52f70101960861e5d995d3a4bec010bc89e6/src/nmath/pgamma.c
 * 
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
 *  Copyright (C) 2005-10 The R Foundation
 *  Copyright (C) 2006-2015 The R Core Team
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  DESCRIPTION
 *
 *    The distribution function of the Poisson distribution.
 */

/*
#include "nmath.h"
#include "dpq.h"

// Rcpp::ppois(n_repeats, E_reads, false);
// Can assume...
// x is an integer and >= 1
// lambda is > 0
// lower_tail is false
// log_p is false
double ppois(double x, double lambda, int lower_tail, int log_p)
{
  if(lambda <= 0.) Rcpp::stop("Lambda must be > 0.");
  if (x < 0) Rcpp::stop("x must be >= 0.");
  x = floor(x + 1e-7); // Why?
  
  return pgamma(lambda, x + 1, 1., !lower_tail, log_p);
}

// x > 0, alph is an integer and >=2, scale=1., lower_tail=true, log_p=false
double pgamma(double x, double alph, double scale, int lower_tail, int log_p)
{
  if(alph <= 0. || scale <= 0.) Rcpp::stop("alph > 0 and scale > 0 are required.")
  x /= scale;
  return pgamma_raw (x, alph, lower_tail, log_p);
}

// x>0 (the original lambda), alph is an integer and >=2, lower_tail=true, log_p=false
double pgamma_raw (double x, double alph, int lower_tail, int log_p)
{
  // Here, assume that  (x,alph) are not NA  &  alph > 0 .
  
  double res;
  
///?  R_P_bounds_01(x, 0., ML_POSINF);
  
  if (x < 1) {
    res = pgamma_smallx (x, alph, lower_tail, log_p);
  } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
    // incl. large alph compared to x 
    double sum = pd_upper_series (x, alph, log_p);// = x/alph + o(x/alph) 
    double d = dpois_wrap (alph, x, log_p);
    
    if (!lower_tail)
      res = log_p
      ? R_Log1_Exp (d + sum)
        : 1 - d * sum;
    else
      res = log_p ? sum + d : sum * d;
  } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
    // incl. large x compared to alph 
    double sum;
    double d = dpois_wrap (alph, x, log_p);

    if (alph < 1) {
      if (x * DBL_EPSILON > 1 - alph)
        sum = R_D__1;
      else {
        double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
        // = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) 
        sum = log_p ? log (f) : f;
      }
    } else {
      sum = pd_lower_series (x, alph - 1);// = (alph-1)/x + o((alph-1)/x) 
        sum = log_p ? log1p (sum) : 1 + sum;
    }
    
    if (!lower_tail)
      res = log_p ? sum + d : sum * d;
    else
      res = log_p
      ? R_Log1_Exp (d + sum)
        : 1 - d * sum;
  } else { // x >= 1 and x fairly near alph. 
    res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
  }
  
  //
  // We lose a fair amount of accuracy to underflow in the cases
  // where the final result is very close to DBL_MIN.	 In those
  // cases, simply redo via log space.
  //
  if (!log_p && res < DBL_MIN / DBL_EPSILON) {
    // with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292
    return exp (pgamma_raw (x, alph, lower_tail, 1));
  } else
    return res;
}

/// dpois_wrap (x__1, lambda) := dpois(x__1 - 1, lambda);  where
 // dpois(k, L) := exp(-L) L^k / gamma(k+1)  {the usual Poisson probabilities}
 //
 // and  dpois*(.., give_log = TRUE) :=  log( dpois*(..) )
 //
static double
  dpois_wrap (double x_plus_1, double lambda, int give_log)
  {
    if (x_plus_1 > 1)
      return dpois_raw (x_plus_1 - 1, lambda, give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
      return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
      double d = dpois_raw (x_plus_1, lambda, give_log);

      return give_log
        ? d + log (x_plus_1 / lambda)
          : d * (x_plus_1 / lambda);
    }
  }
*/
