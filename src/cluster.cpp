#include <Rcpp.h>
#include <RcppParallel.h>
#include "dada.h"
// [[Rcpp::interfaces(cpp)]]

/********* ALGORITHM LOGIC *********/

/*
 compare:
Performs alignments and computes lambda for all raws to the specified Bi
Stores only those that can possibly be recruited to this Bi
*/
void b_compare(B *b, unsigned int i, Rcpp::NumericMatrix errMat, int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, bool vectorized_alignment, int SSE, bool gapless, bool verbose) {
  unsigned int index, cind;
  double lambda;
  Raw *raw;
//  Raw *center = b->bi[i]->center;
//  uint8_t *center_kmer8 = center->kmer8; // Store original kmer8 pointer
  Sub *sub;
  Comparison comp;
  /* Testing code for keeping comparative kmer8 close to others. Isn't improving perf thus far.
  // Caching current kmer8 w/in the k8 memory block every CACHE_SIZE strides
  size_t cached;
  size_t n_kmer = 1 << (2*KMER_SIZE);
  uint8_t *kcache = (uint8_t *) malloc(n_kmer * sizeof(uint8_t)); //E
  if (kcache == NULL)  Rcpp::stop("Memory allocation failed.");
  cached=CACHE_STRIDE-1;
  // Cache
  if(cached < b->nraw) {
    memcpy(kcache, b->raw[cached]->kmer8, n_kmer);
    memcpy(b->raw[cached]->kmer8, center->kmer8, n_kmer);
    center->kmer8 = b->raw[cached]->kmer8;
  } */

  // align all raws to this sequence and compute corresponding lambda
  if(verbose) { Rprintf("C%iLU:", i); }
  for(index=0, cind=0; index<b->nraw; index++) {
    raw = b->raw[index];
    
/*    if(index == cached) { // Return real data
      memcpy(raw->kmer8, kcache, n_kmer);
      if((cached+=CACHE_STRIDE) < b->nraw) { // Move cache forward
        memcpy(kcache, raw->kmer8, n_kmer);
        memcpy(raw->kmer8, center_kmer8, n_kmer);
        center->kmer8 = b->raw[cached]->kmer8;
      } else {
        center->kmer8 = center_kmer8;
      }
} */

    // get sub object
    sub = sub_new(b->bi[i]->center, raw, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless);
    b->nalign++;
    if(!sub) { b->nshroud++; }
    
    // Calculate lambda for that sub
    lambda = compute_lambda(raw, sub, errMat, b->use_quals, errMat.ncol());

    // Store lambda and set self
    if(index == b->bi[i]->center->index) { b->bi[i]->self = lambda; }

    // Store comparison if potentially useful
    if(lambda * b->reads > raw->E_minmax) { // This cluster could attract this raw
      if(lambda * b->bi[i]->center->reads > raw->E_minmax) { // Better E_minmax, set
        raw->E_minmax = lambda * b->bi[i]->center->reads;
      }
      comp.i = i;
      comp.index = index;
      comp.lambda = lambda;
      comp.hamming = sub->nsubs;
      b->bi[i]->comp.push_back(comp);
      if(i==0 || raw == b->bi[i]->center) { // Update on init (i=0) or if the center (as b_bud doesn't update raw->comp)
        raw->comp = comp;
      } /// Handle init better?
    }
    sub_free(sub);
  }
}

struct CompareParallel : public RcppParallel::Worker
{
  // source data
  B *b;
  unsigned int i;
  
  // destination comparison array
  Comparison *output;
  
  // parameters
  double *err_mat;
  unsigned int ncol;
  int match, mismatch, gap_pen, homo_gap_pen;
  bool use_kmers;
  double kdist_cutoff;
  int band_size;
  bool vectorized_alignment;
  int SSE;
  bool gapless;
  
  // initialize with source and destination
  CompareParallel(B *b, unsigned int i, double *err_mat, unsigned int ncol, Comparison *output,
                  int match, int mismatch, int gap_pen, int homo_gap_pen,
                  bool use_kmers, double kdist_cutoff, int band_size, bool vectorized_alignment, int SSE, bool gapless) 
    : b(b), i(i), err_mat(err_mat), ncol(ncol), output(output), match(match), mismatch(mismatch), gap_pen(gap_pen), homo_gap_pen(homo_gap_pen), 
      use_kmers(use_kmers), kdist_cutoff(kdist_cutoff), band_size(band_size), vectorized_alignment(vectorized_alignment), SSE(SSE), gapless(gapless) {}
  
  // Perform sequence comparison
  void operator()(std::size_t begin, std::size_t end) {
    Raw *raw;
    Sub *sub;
    
    for(std::size_t index=begin;index<end;index++) {
      raw = b->raw[index];
      sub = sub_new(b->bi[i]->center, raw, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless);

      // Make comparison object
      output[index].i = i;
      output[index].index = index;
      output[index].lambda = compute_lambda_ts(raw, sub, ncol, err_mat, b->use_quals);
      if(sub) {
        output[index].hamming = sub->nsubs;
      } else {
        output[index].hamming = -1;
      }

      // Free sub
      sub_free(sub);
    }
  }
};


void b_compare_parallel(B *b, unsigned int i, Rcpp::NumericMatrix errMat, int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, bool vectorized_alignment, int SSE, bool gapless, bool verbose) {
  unsigned int index, cind, row, col, ncol;
  double lambda;
  Raw *raw;
  Comparison comp;
  
  // Make thread-safe C-array for the error rate matrix
  double *err_mat = (double *) malloc(sizeof(double) * errMat.ncol() * errMat.nrow());
  if(err_mat==NULL) Rcpp::stop("Memory allocation failed.");
  ncol = errMat.ncol();
  if(errMat.nrow() != 16) { Rcpp::stop("Error matrix doesn't have 16 rows."); }
  for(row=0;row<errMat.nrow();row++) {
    for(col=0;col<errMat.ncol();col++) {
      err_mat[row*ncol + col] = errMat(row, col);
    }
  }
  
  // Parallelize for loop to perform all comparisons
  Comparison *comps = (Comparison *) malloc(sizeof(Comparison) * b->nraw);
  if(comps==NULL) Rcpp::stop("Memory allocation failed.");
  CompareParallel compareParallel(b, i, err_mat, ncol, comps, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless);
  RcppParallel::parallelFor(0, b->nraw, compareParallel, GRAIN_SIZE);
  
  // Selectively store
  for(index=0, cind=0; index<b->nraw; index++) {
    b->nalign++; ///t
    raw = b->raw[index];
    comp = comps[index];
    lambda = comp.lambda;
    if(lambda<0 || lambda>1) Rcpp::stop("Lambda out-of-range error.");

    // Store self-lambda
    if(index == b->bi[i]->center->index) { 
      b->bi[i]->self = lambda; 
    }
    
    // Store comparison if potentially useful
    if(lambda * b->reads > raw->E_minmax) { // This cluster could attract this raw
      if(lambda * b->bi[i]->center->reads > raw->E_minmax) { // Better E_minmax, set
        raw->E_minmax = lambda * b->bi[i]->center->reads;
      }
      b->bi[i]->comp.push_back(comp);
      if(i==0 || raw == b->bi[i]->center) { // Update on init (i=0) or if the center (as b_bud doesn't update raw->comp)
        raw->comp = comp;
      } // Handle init better?
    }
  }
  free(err_mat);
  free(comps);
}

/* b_shuffle2:
 move each sequence to the bi that produces the highest expected
 number of that sequence. The center of a Bi cannot leave.
*/
bool b_shuffle2(B *b) {
  unsigned int i, cind, index;
  double e;
  bool shuffled = false;
  Comparison *comp;
  Raw *raw;
  
  double *emax = (double *) malloc(b->nraw * sizeof(double)); //E
  Comparison **compmax = (Comparison **) malloc(b->nraw * sizeof(Comparison *)); //E
  if(emax==NULL || compmax==NULL) Rcpp::stop("Memory allocation failed.");

  // Initialize emax/imax off of cluster 0
  // Comparisons to all raws exist in cluster 0, in index order
  for(index=0;index<b->nraw;index++) {
    compmax[index] = &b->bi[0]->comp[index];
    emax[index] = compmax[index]->lambda * b->bi[0]->reads;
  }
  
  // Iterate over comparisons, find comparison with best E for each raw
  for(i=1;i<b->nclust;i++) {
    for(cind=0;cind<b->bi[i]->comp.size();cind++) {
      comp = &b->bi[i]->comp[cind];
      index = comp->index;
      e = comp->lambda * b->bi[i]->reads;
      if(e > emax[index]) { // better E
        compmax[index] = comp;
        emax[index] = e;
      }
    }
  }
  
  // Iterate over raws, if best i different than current, move
  for(i=0;i<b->nclust;i++) {
    // IMPORTANT TO ITERATE BACKWARDS DUE TO BI_POP_RAW!!!!!!
    for(int r=b->bi[i]->nraw-1; r>=0; r--) {
      raw = b->bi[i]->raw[r];
      // If a better cluster was found, move the raw to the new bi
      if(compmax[raw->index]->i != i) {
        if(raw->index == b->bi[i]->center->index) {  // Check if center
          if(VERBOSE) { Rprintf("Warning: Shuffle blocked the center of a Bi from leaving."); }
          continue;
        }
        // Moving raw
        bi_pop_raw(b->bi[i], r);
        bi_add_raw(b->bi[compmax[raw->index]->i], raw);
        // Assign raw the Comparison of its new cluster
        raw->comp = *compmax[raw->index];
        shuffled = true;  
      }  
    } // for(r=0;r<b->bi[i]->nraw;r++)
  }

  free(compmax);
  free(emax);
  
  return shuffled;
}

/* b_p_update:
 Calculates the abundance p-value for each raw in the clustering.
 Depends on the lambda between the raw and its cluster, and the reads of each.
*/
void b_p_update(B *b) {
  unsigned int i, r;
  Raw *raw;
  for(i=0;i<b->nclust;i++) {
    if(b->bi[i]->update_e) {
      for(r=0;r<b->bi[i]->nraw;r++) {
        raw = b->bi[i]->raw[r];
        raw->p = get_pA(raw, b->bi[i]);
      } // for(r=0;r<b->bi[i]->nraw;r++)
      b->bi[i]->update_e = false;
    }
  } // for(i=0;i<b->nclust;i++)
}

/* b_bud:
 Finds the minimum p-value. If significant, creates a new cluster and moves the
 raws from the raw with the minimum p-value to the new cluster.
 Returns index of new cluster, or 0 if no new cluster added.
*/

int b_bud(B *b, double min_fold, int min_hamming, int min_abund, bool verbose) {
  int i, r, hamming;
  int mini = -1, minr = -1, mini_prior = -1, minr_prior = -1; // Negative initialization
  double pA, pP;
  double lambda, expected;
  Raw *raw, *minraw, *minraw_prior;
  minraw = b->bi[0]->center; // Assumes that complete alignment/pval calcs were performed in the init cluster
  minraw_prior = b->bi[0]->center; // Assumes that complete alignment/pval calcs were performed in the init cluster
  
  // Find i, r indices and value of minimum pval.
  for(i=0;i<b->nclust;i++) {
    for(r=1; r<b->bi[i]->nraw; r++) { // r=0 is the center
      raw = b->bi[i]->raw[r];
///      if(b->bi[i]->center->index == raw->index) { continue; } // Don't bud centers
      if(raw->reads < min_abund) { continue; }
      hamming = raw->comp.hamming;
      lambda = raw->comp.lambda;

      // Calculate the fold over-abundance and the hamming distance to this raw
      if(hamming >= min_hamming) { // Only those passing the hamming/fold screens can be budded
        if(min_fold <= 1 || ((double) raw->reads) >= min_fold * lambda * b->bi[i]->reads) {  
          if((raw->p < minraw->p) ||
             ((raw->p == minraw->p && raw->reads > minraw->reads))) { // Most significant
            mini = i; minr = r;
            minraw = raw;
          }
          if(raw->prior && ((raw->p < minraw_prior->p) ||
             (raw->p == minraw_prior->p && raw->reads > minraw_prior->reads))) { // Most significant
            mini_prior = i; minr_prior = r;
            minraw_prior = raw;
          }
        }
      }
    }
  }

  // Bonferoni correct the abundance pval by the number of raws and compare to OmegaA
  // (quite conservative, although probably unimportant given the abundance model issues)
  pA = minraw->p * b->nraw;
  pP = minraw_prior->p;
  if(pA < b->omegaA && mini >= 0) {  // A significant abundance pval
    expected = minraw->comp.lambda * b->bi[mini]->reads;
    raw = bi_pop_raw(b->bi[mini], minr);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "A");
    b->bi[i]->birth_pval = pA;
    b->bi[i]->birth_fold = raw->reads/expected;
    b->bi[i]->birth_e = expected;
    b->bi[i]->birth_comp = minraw->comp;
    
    // Add raw to new cluster.
    bi_add_raw(b->bi[i], raw);
    bi_assign_center(b->bi[i]);
    return i;
  } else if (pP < b->omegaP && mini_prior >= 0) {  // A significant prior-abundance pval
    expected = minraw_prior->comp.lambda * b->bi[mini_prior]->reads;
    raw = bi_pop_raw(b->bi[mini_prior], minr_prior);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "P");
    b->bi[i]->birth_pval = pP;
    b->bi[i]->birth_fold = raw->reads/expected;
    b->bi[i]->birth_e = expected;
    b->bi[i]->birth_comp = minraw_prior->comp;
    
    // Add raw to new cluster.
    bi_add_raw(b->bi[i], raw);
    bi_assign_center(b->bi[i]);
    return i;
  }

  return 0;
}

/********* CONTAINER HOUSEKEEPING *********/

// Iterate over the raws in a bi and update reads/nraw
void bi_census(Bi *bi) {
  unsigned int r, reads=0, nraw=0;
  for(r=0;r<bi->nraw;r++) {
    reads += bi->raw[r]->reads;
    nraw++;
  }
  if(reads != bi->reads) {
    bi->update_e = true;
  }
  bi->reads = reads;
  bi->nraw = nraw;
}

// Takes a Bi object, and calculates and assigns its center Raw.
// Currently this is done trivially by choosing the most abundant raw.
// This function also currently assigns the cluster sequence (equal to center->seq).
void bi_assign_center(Bi *bi) {
  unsigned int r, max_reads;
  
  // Assign the raw with the most reads as the center
  bi->center = NULL;
  for(r=0,max_reads=0;r<bi->nraw;r++) {
    if(bi->raw[r]->reads > max_reads) { // Most abundant
      bi->center = bi->raw[r];
      max_reads = bi->center->reads;
    }
  }
  // Assign center sequence to bi->seq
  if(bi->center) { strcpy(bi->seq, bi->center->seq); }
}

