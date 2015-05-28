/*
Pval.cpp contains the functions related to calculating the abundance and singleton pvals in DADA.
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

// Find abundance pval from a Fam in a Bi
double get_pA(Fam *fam, Bi *bi) {
  double E_reads, pval = 1.;
  
  if(fam->reads == 1) {   // Singleton. No abundance pval.
    pval=1.;
  } 
  else if(!(fam->sub)) { // Outside kmer threshhold
    pval=0.;
  } 
  else if(fam->sub->nsubs == 0) { // Cluster center
    pval=1.;
  }
  else if(fam->lambda == 0) { // Zero expected reads of this fam
    pval = 0.;
  }
  else { // Calculate abundance pval.
    // E_reads is the expected number of reads for this fam
    E_reads = fam->lambda * bi->reads;
    pval = calc_pA(fam->reads, E_reads);
  }
  return pval;
}

/*
 compute_lambda: DEPRECATED
 
 parameters:
 char *seq, the reference sequence away from which a lambda will be computed
 Sub *sub, a substitution struct
 double self, the rate of self-production of the consensus sequence
 double t[4][4], a 4x4 matrix of content-independent error probabilities
 
 returns:
 the double lambda.
 */
double compute_lambda(Sub *sub, double self, double t[4][4], bool use_quals) {
  int s, nti0, nti1;
  double lambda, trans, qual;
  
  if(!sub) { // NULL Sub, outside Kmer threshold
    return 0.0;
  }

  lambda = self;
  for(s=0; s < sub->nsubs; s++) {
    nti0 = (int)sub->nt0[s] - 1;
    nti1 = (int)sub->nt1[s] - 1;
    trans = t[nti0][nti1] / t[nti0][nti0];
    
    if(use_quals) {
      if(!sub->q1) {
        printf("Warning: Missing quality information when computing lambda.\n");
      }
      else {     // Incorporate the qualities. HACKY FOR NOW
        qual = sub->q1[s];
        if(qual<36) {
          trans = trans * pow(10.0, (36.0-qual)/6.5);
        }
        if(trans > 0.33) { trans = 0.33; } // capping at 0.33 because?
      }
    }
    
    lambda = lambda * trans;
  }

  if(lambda < 0 || lambda > 1) { printf("Error: Over- or underflow OF lambda: %.4e\n", lambda); }

  return lambda;
}

// This calculates lambda from a lookup table index by transition (row) and rounded quality (col)
double compute_lambda3(Raw *raw, Sub *sub, Rcpp::NumericMatrix errMat, bool use_quals) {
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
      Rcpp::stop("Error: Can't handle non ACGT sequences in CL3.\n");
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
    if(pos0 < 0 || pos0 >= len1) {
      Rcpp::stop("CL3: Bad pos0 (%i)", pos0);
    }
    pos1 = sub->map[sub->pos[s]];
    if(pos1 < 0 || pos1 >= len1) {
      Rcpp::stop("CL3: Bad pos1 (%i)", pos1);
    }
    
    nti0 = ((int) sub->nt0[s]) - 1;
    nti1 = ((int) sub->nt1[s]) - 1;
    tvec[pos1] = nti0*4 + nti1;
  }

  // And calculate lambda
  lambda = 1.0;
  for(pos1=0;pos1<len1;pos1++) {
    lambda = lambda * errMat(tvec[pos1], qind[pos1]);
  }

  if(lambda < 0 || lambda > 1) { Rcpp::stop("Bad lambda"); }

  return lambda;
}

// Find singleton pval for a Fam in a Bi, also have to pass in B for lambda/cdf lookup
double get_pS(Fam *fam, Bi *bi, B *b) {
  size_t ifirst, imid, ilast;
  double pval;

  // Calculate singleton pval (pS) from cluster lookup
  if(fam->lambda >= b->lams[0]) {  // fam->lambda bigger than all lambdas in lookup
    pval = 1.0;
  }
  else if(fam->lambda <= b->lams[b->nlam-1]) { // fam->lam smaller than all lams in lookup
    pval = (1.0 - b->cdf[b->nlam-1]);
  }
  else { // Find lam in the lookup and assign pS
    ifirst = 0;
    ilast = b->nlam-1;
    while((ilast-ifirst) > 1) {
      imid = (ifirst+ilast)/2;
      if(b->lams[imid] > fam->lambda) {
        ifirst = imid;
      } else {
        ilast = imid;
      }
    }
    pval = (1.0 - b->cdf[ifirst]);
  }
  return pval;
}

void getCDF(std::vector<double>& ps, std::vector<double>& cdf, double err[4][4], int nnt[4], int maxD) {
  int i, j, k, d;
  size_t index;
  
  // Copy err "matrix" into relative errs vector
  k=0;
  double errs[NERRS];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      if(i!=j) { errs[k++] = err[i][j]; }
    }
  }

  std::vector<Prob> probs;

  // Get prob of error for each nt, and the self-trans prob
  double pa = errs[0]+errs[1]+errs[2];
  double pc = errs[3]+errs[4]+errs[5];
  double pg = errs[6]+errs[7]+errs[8];
  double pt = errs[9]+errs[10]+errs[11];
  double self = pow((1.-pa), nnt[0]) * pow((1.-pc), nnt[1]) * pow((1.-pg), nnt[2]) * pow((1.-pt), nnt[3]);

  // make errors relative to the non-err prob
  for(i=0;i<3;i++) { errs[i] = errs[i]/(1.0 - pa); }
  for(i=3;i<6;i++) { errs[i] = errs[i]/(1.0 - pc); }
  for(i=6;i<9;i++) { errs[i] = errs[i]/(1.0 - pg); }
  for(i=9;i<12;i++) { errs[i] = errs[i]/(1.0 - pt); }
  
  // Declare variables needed to iterate though all d-aways
  int nerr[NERRS];
  double nopen[4];
  int first, store;
  double p, n;

  for(d=0;d<=maxD;d++) {
    // init partition
    for(i=0;i<NERRS;i++) { nerr[i] = 0; }
    nerr[NERRS-1] = d;
    first = NERRS-1;
    if(VERBOSE) {
      Rcpp::Rcout << "---- D = " << d << " ----\n";
    }

    while(1) {
      // Calc and store p/n/pval contribution
      p = self;
      n = 1;
      for(i=0;i<4;i++) { nopen[i] = nnt[i]; }
      
      for(i=0;i<NERRS;i++) {
        if(nerr[i]>0) {
          for(j=0;j<nerr[i];j++) {
            p = p*errs[i];
            n = n*nopen[i/3]/(j+1.0);
            nopen[i/3]--; // one fewer of that base available for future errors
            // going below zero no prob, since the times zero makes everything zero
          }
        }
      }
      probs.push_back(Prob(p,n));

      if(nerr[0] >= d) { break; } // Should break when all d in first
      
      // Advance to next partition
      if(first > 0) {
        nerr[first]--;
        first--;
        nerr[first]++;
      } else {
        store = nerr[0];
        nerr[0] = 0;
        for(first=1;first<NERRS;first++) {
          if(nerr[first]>0) { break; }
          if(first==(NERRS-1)) {  // Error checking
            printf("Error: Partition iteration failed in getCDF!\n");
            return;
          }
        }
        nerr[first]--;
        first--;
        nerr[first] += (store+1);
      }
    } // while(1)

  } // for(d=0;d<maxD;d++)

  // Sort probs in descending order and create pval vector
  std::sort(probs.rbegin(), probs.rend());   // note: reverse iterators -> decreasing sort
  // Note that the 15-17 digit significand precision of doubles sets a limit on the precision of the 1-cdf tail.

  ps.resize(0);
  cdf.resize(0);
  double cum = 0;
  for(index=0;index < probs.size();index++) {
    cum += (probs[index].first * probs[index].second);
    ps.push_back(probs[index].first);
    cdf.push_back(cum);
  }

  // Find largest p beyond maxD (= max_err**(maxD+1))
  double max_err = 0.0;
  for(i=0;i<12;i++) { if(errs[i] > max_err) max_err = errs[i]; }
  double min_p = pow(max_err, maxD+1);
  
  for(index=0;index < ps.size();index++) {
    if(ps[index] <= min_p) { break; }
  }
  ps.resize(index);
  cdf.resize(index);
}

// I AM BROKEN RIGHT NOW BECAUSE OF CHANGE IN ERR MATRIX STUFF!!!!!!!!!
void b_make_pS_lookup(B *b) {
///!  static double err[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
  static int ave_nnt[4] = {0, 0, 0, 0};
  static double *lams = NULL;   // the lookup table
  static double *cdf = NULL;    // the lookup table
  static int nlam = 0;          // the lookup table
  static double omegaS = 0.0;
  int i, j, nti;
  int new_ave_nnt[4];
  size_t index;
  
  // Error and exit if requested OmegaS would exceed or near double precision
  if(b->reads * DBL_PRECISION > b->omegaS) {
    printf("Error: Doubles not precise enough to meet requested OmegaS.\n");
    printf("       Re-run DADA with singletons turned off or a less stringent OmegaS.\n");
    
    b_free(b);
    Rcpp::stop("Cannot meet requested OmegaS\n");
  }    

  // Calculate average sequence nnt
  // RIGHT NOW THIS IS THE AVERAGE OVER RAWS NOT OVER READS: RIGHT CHOICE?????
  double tot_nnt[] = {0.0,0.0,0.0,0.0};
  for (index = 0; index < b->nraw; index++) {
    for(i=0;i<b->raw[index]->length;i++) {
      nti = ((int) b->raw[index]->seq[i]) - 1;
      if(nti == 0 || nti == 1 || nti ==2 || nti == 3) {
        tot_nnt[nti]++;
      }
    }
  }
  int nnt=0;
  int del_ave_nnt = 0;
  for(nti=0;nti<4;nti++) {
    new_ave_nnt[nti] = (int) (0.4999 + tot_nnt[nti]/b->nraw);
    del_ave_nnt += abs(new_ave_nnt[nti] - ave_nnt[nti]);
    nnt += new_ave_nnt[nti];
  }
  
  int err_diffs = 0;
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
///!      if(err[i][j] != b->err[i][j]) { err_diffs++; }
    }
  }
  
  // Check if lookup parameters are same as last time, and if so use the stored lookup
  if(lams && cdf) {  // Already a lookup stored (not NULL)
    if(err_diffs == 0 && b->omegaS == omegaS && del_ave_nnt < nnt/10) { // same enough
      b->lams = lams;
      b->cdf = cdf;
      b->nlam = nlam;
      return;
    } else {
      free(lams);
      free(cdf);
    }
  }
  
  // Making a new lookup, so store the new params
  omegaS = b->omegaS;
  for(i=0;i<4;i++) {
    ave_nnt[i] = new_ave_nnt[i];
    for(j=0;j<4;j++) {
///!      err[i][j] = b->err[i][j];
    }
  }

  // Iterate over maxDs until going far enough to call significant singletons
  // Approximating Bonferonni correction by nraw (in place of total fams)
  // i.e. most significant possible pS* = (1.0 - temp_cdf.back()) * b->nraw
  // DADA_ML MATCH: maxD = 10
  int maxD=8;
  std::vector<double> temp_lambdas;
  std::vector<double> temp_cdf;
  do {
    maxD+=2;
///!    getCDF(temp_lambdas, temp_cdf, b->err, ave_nnt, maxD);
  } while((1.0 - temp_cdf.back()) * b->reads > b->omegaS && maxD < MAXMAXD);
  
  // Warn if couldnt make lookup big enough to get OmegaS
  if((1.0 - temp_cdf.back()) * b->reads > b->omegaS) {
    printf("Warning: Cannot calculate singleton pvals small enough to meet requested OmegaS.\n");
    printf("         Running DADA with singletons turned off.\n");
    b->lams = NULL;
    b->cdf = NULL;
    b->use_singletons = false;
    return;
  }
  
  // Copy into C style arrays
  // Kind of silly, at some point might be worthwhile doing the full C++ conversion
//  if(VERBOSE) { printf("b_new: The least most significant possible pval = %.4e, pS* ~ %.4e (maxD=%i, ave_nnt=%i,%i,%i,%i)\n", 1.0-(temp_cdf.back()), b->reads*(1.0-(temp_cdf.back())), maxD, ave_nnt[0], ave_nnt[1], ave_nnt[2], ave_nnt[3]); }
  lams = (double *) malloc(temp_lambdas.size() * sizeof(double)); //E
  if (lams == NULL)  Rcpp::stop("Memory allocation failed!\n");
  cdf = (double *) malloc(temp_cdf.size() * sizeof(double)); //E
  if (cdf == NULL)  Rcpp::stop("Memory allocation failed!\n");
  
  nlam = temp_lambdas.size();
  for(index=0;index<nlam;index++) {
    lams[index] = temp_lambdas[index];
    cdf[index] = temp_cdf[index];
  }
  
  // Assign to the B object
  b->lams = lams;
  b->cdf = cdf;
  b->nlam = nlam;
}

// [[Rcpp::export]]
Rcpp::DataFrame getSingletonCDF(Rcpp::NumericMatrix err, std::vector<int> nnt, int maxD) {
  int i, j;
  // Copy err into a C style array
  if(err.nrow() != 4 || err.ncol() != 4) {
    Rcpp::Rcout << "C: Error matrix malformed:" << err.nrow() << ", " << err.ncol() << "\n";
    return R_NilValue;
  }

  double c_err[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_err[i][j] = err(i,j);
    }
  }
  
  int c_nnt[4];
  for(i=0;i<4;i++) { c_nnt[i] = nnt[i]; }

  std::vector<double> ps;
  std::vector<double> cdf;
  getCDF(ps, cdf, c_err, c_nnt, maxD);
  
  return Rcpp::DataFrame::create(Rcpp::_["p"]=ps, Rcpp::_["cdf"]=cdf);
}

