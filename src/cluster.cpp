#include <Rcpp.h>
#include <RcppParallel.h>
#include "dada.h"
// [[Rcpp::interfaces(cpp)]]

/*
 methods for "B" objects.
 The "B" object is the partition of a DNA data set into
 clusters that is updated until convergence during DADA's
 inner two loops.
 */

#define RAWBUF 50
#define CLUSTBUF 50

/* private function declarations */
Bi *bi_new(unsigned int totraw);
void bi_free(Bi *bi);

unsigned int b_add_bi(B *b, Bi *bi);
Raw *bi_pop_raw(Bi *bi, unsigned int r);
unsigned int bi_add_raw(Bi *bi, Raw *raw);

void bi_census(Bi *bi);
void bi_assign_center(Bi *bi);
void bi_make_consensus(Bi *bi, bool use_quals);


/********* CONSTRUCTORS AND DESTRUCTORS *********/

// The constructor for the Raw object.
Raw *raw_new(char *seq, double *qual, unsigned int reads) {
  // Allocate
  Raw *raw = (Raw *) malloc(sizeof(Raw)); //E
  if (raw == NULL)  Rcpp::stop("Memory allocation failed.");
  raw->seq = (char *) malloc(strlen(seq)+1); //E
  if (raw->seq == NULL)  Rcpp::stop("Memory allocation failed.");
  // Assign sequence and associated properties
  strcpy(raw->seq, seq);
  raw->length = strlen(seq);
///  raw->kmer = get_kmer(seq, KMER_SIZE);
///  raw->kmer8 = get_kmer8(seq, KMER_SIZE);
  raw->kord = get_kmer_order(seq, KMER_SIZE);
  raw->reads = reads;
  // Allocate and copy quals (quals downgraded to floats here for memory savings)
  if(qual) { 
    raw->qual = (uint8_t *) malloc(raw->length * sizeof(uint8_t)); //E
    if (raw->qual == NULL)  Rcpp::stop("Memory allocation failed.");
    for(size_t i=0;i<raw->length;i++) { raw->qual[i] = (uint8_t) round(qual[i]); }
  } else {
    raw->qual = NULL;
  }
  raw->E_minmax = -999.0;
  return raw;
}

// The destructor for the Raw object.
void raw_free(Raw *raw) {
  free(raw->seq);
  if(raw->qual) { free(raw->qual); }
///  free(raw->kmer);
///  free(raw->kmer8);
  free(raw->kord);
  free(raw);  
}

// The constructor for the Bi object.
Bi *bi_new(unsigned int totraw) {
  Bi *bi = new Bi;
  if (bi == NULL)  Rcpp::stop("Memory allocation failed!\n");
  bi->raw = (Raw **) malloc(RAWBUF * sizeof(Raw *)); //E
  if (bi->raw == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->maxraw = RAWBUF;
  bi->totraw = totraw;
  bi->center = NULL;
  strcpy(bi->seq, "");
  bi->update_e = true;
  bi->shuffle = true;
  
  bi->reads = 0;
  bi->nraw = 0;
  return bi;
}

// The destructor for the Bi object
void bi_free(Bi *bi) {
  free(bi->raw);
  delete bi;
}

// The constructor for the B object. Takes in array of Raws.
B *b_new(Raw **raws, unsigned int nraw, int score[4][4], int gap_pen, int homo_gap_pen, double omegaA, int band_size, bool vectorized_alignment, bool use_quals) {
  unsigned int i, j, index;

  // Allocate memory
  B *b = (B *) malloc(sizeof(B)); //E
  if (b == NULL)  Rcpp::stop("Memory allocation failed.");
  b->bi = (Bi **) malloc(CLUSTBUF * sizeof(Bi *)); //E
  if (b->bi == NULL)  Rcpp::stop("Memory allocation failed.");
  b->maxclust = CLUSTBUF;
  
  // Initialize basic values
  b->nclust = 0;
  b->reads = 0;
  b->nraw = nraw;
  b->gap_pen = gap_pen;
  b->homo_gap_pen = homo_gap_pen;
  b->omegaA = omegaA;
  b->band_size = band_size;
  b->vectorized_alignment = vectorized_alignment;
  b->use_quals = use_quals;
  
  // Copy the score matrix
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      b->score[i][j] = score[i][j];
    }
  }

  // Initialize array of raws
  // Remember that these were allocated outside the construction of B
  b->raw = raws;
  for (index = 0; index < b->nraw; index++) {
    b->raw[index]->index = index;
    b->reads += b->raw[index]->reads;
  }

  // Initialize with one cluster containing all the raws.
  b_init(b);
  return b;
}

// Initialize B object by placing all Raws into one Bi
// Called on B construction
void b_init(B *b) {
  unsigned int i, index;
  
  // Destruct existing clusters
  for(i=0; i<b->nclust; i++) {
    bi_free(b->bi[i]);
  }
  b->nclust=0;
  
  // Add the one cluster and initialize its "birth" information
  b_add_bi(b, bi_new(b->nraw));
  strcpy(b->bi[0]->birth_type, "I");
  b->bi[0]->birth_pval = 0.0;
  b->bi[0]->birth_fold = 1.0;
  b->bi[0]->birth_e = b->reads;
  b->nalign = 0;
  b->nshroud = 0;

  // Add all raws to that cluster
  for (index=0; index<b->nraw; index++) {
    bi_add_raw(b->bi[0], b->raw[index]);
  }

  bi_census(b->bi[0]);
  bi_assign_center(b->bi[0]); // Makes cluster center sequence
}

// The destructor for the B object.
void b_free(B *b) {
  for(int i=0;i<b->nclust;i++) { bi_free(b->bi[i]); }
  free(b->bi);
  free(b);
}


/********* CONTAINER OPERATIONS *********/

// Add a Raw to a Bi object. Update reads/nraw. Return index to raw in bi.
unsigned int bi_add_raw(Bi *bi, Raw *raw) {
  // Allocate more space if needed
  if(bi->nraw >= bi->maxraw) {    // Extend Raw* buffer
    bi->raw = (Raw **) realloc(bi->raw, (bi->maxraw+RAWBUF) * sizeof(Raw *)); //E
    if (bi->raw == NULL)  Rcpp::stop("Memory allocation failed.");
    bi->maxraw+=RAWBUF;
  }
  // Add raw and update reads/nraw
  bi->raw[bi->nraw] = raw;
  bi->reads += raw->reads;
  bi->update_e = true;
  return(bi->nraw++);
}

// Adds a Bi to a B object. Update nclust. Returns index to bi in b.
// This should only ever used to add new, empty Bi objects.
unsigned int b_add_bi(B *b, Bi *bi) {
  // Allocate more space if needed
  if(b->nclust >= b->maxclust) {    // Extend Bi* buffer
    b->bi = (Bi **) realloc(b->bi, (b->maxclust+CLUSTBUF) * sizeof(Bi *)); //E
    if (b->bi == NULL)  Rcpp::stop("Memory allocation failed.");
    b->maxclust+=CLUSTBUF;
  }
  // Add bi and update nclust
  b->bi[b->nclust] = bi;
  bi->i = b->nclust;
  return(b->nclust++);
}

// Removes a Raw from a Bi object. Updates reads/nraw. Returns pointer to that raw.
// Sets update_e flag. Frees the bi/raw sub object.
// THIS REORDERS THE BIs LIST OF RAWS. CAN BREAK INDICES.
// This is called in b_shuffle
Raw *bi_pop_raw(Bi *bi, unsigned int r) {
  Raw *pop;
  if(r<bi->nraw) {
    pop = bi->raw[r];
    bi->raw[r] = bi->raw[bi->nraw-1]; // POPPED RAW REPLACED BY LAST RAW
    bi->raw[bi->nraw-1] = NULL;
    bi->nraw--;
    bi->reads -= pop->reads;
    bi->update_e = true;
  } else {
    Rcpp::stop("Container Error (Bi): Tried to pop out-of-range raw.");
    pop = NULL;
  }
  return pop;
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

/********* ALGORITHM LOGIC *********/

/*
 compare:
Performs alignments and computes lambda for all raws to the specified Bi
Stores only those that can possibly be recruited to this Bi
*/
void b_compare(B *b, unsigned int i, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat, bool verbose, int SSE) {
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
    sub = sub_new(b->bi[i]->center, raw, b->score, b->gap_pen, b->homo_gap_pen, use_kmers, kdist_cutoff, b->band_size, b->vectorized_alignment, SSE);
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
  bool use_kmers;
  double kdist_cutoff;
  int SSE;
  unsigned int ncol;
  double *err_mat;
  
  // initialize with source and destination
  CompareParallel(B *b, unsigned int i, Comparison *output, bool use_kmers, double kdist_cutoff, int SSE,
                  unsigned int ncol, double *err_mat) 
    : b(b), i(i), output(output), use_kmers(use_kmers), kdist_cutoff(kdist_cutoff), SSE(SSE), ncol(ncol), err_mat(err_mat) {}
  
  // Perform sequence comparison
  void operator()(std::size_t begin, std::size_t end) {
    Raw *raw;
    Sub *sub;
    
    for(std::size_t index=begin;index<end;index++) {
      raw = b->raw[index];
      sub = sub_new(b->bi[i]->center, raw, b->score, b->gap_pen, b->homo_gap_pen, use_kmers, kdist_cutoff, b->band_size, b->vectorized_alignment, SSE);
      
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


void b_compare_parallel(B *b, unsigned int i, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat, bool verbose, int SSE) {
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
  CompareParallel compareParallel(b, i, comps, use_kmers, kdist_cutoff, SSE, ncol, err_mat);
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
  int i, r;
  int mini, minr, minreads, hamming;
  double minp = 1.0;
  double pA=1.0;
  double lambda, mine;
  Comparison mincomp;
  Raw *raw;

  // Find i, r indices and value of minimum pval.
  mini=-999; minr=-999; minreads=0; minp=1.0;
  for(i=0;i<b->nclust;i++) {
    for(r=0; r<b->bi[i]->nraw; r++) {
      raw = b->bi[i]->raw[r];
      if(b->bi[i]->center->index == raw->index) { continue; } // Don't bud centers
      if(raw->reads < min_abund) { continue; }
      hamming = raw->comp.hamming;
      lambda = raw->comp.lambda;

      // Calculate the fold over-abundance and the hamming distance to this raw
      if(hamming >= min_hamming) { // Only those passing the hamming/fold screens can be budded
        if(min_fold <= 1 || ((double) raw->reads) >= min_fold * lambda * b->bi[i]->reads) {  
          if((raw->p < minp) ||
            ((raw->p == minp && raw->reads > minreads))) { // Most significant
            mini = i; minr = r;
            minp = raw->p;
            mine = lambda * b->bi[i]->reads;
            mincomp = raw->comp;
            minreads = raw->reads;
          }
        }
      }
    }
  }

  // Bonferoni correct the abundance pval by the number of raws and compare to OmegaA
  // (quite conservative, although probably unimportant given the abundance model issues)
  pA = minp*b->nraw;
  if(pA < b->omegaA && mini >= 0 && minr >= 0) {  // A significant abundance pval
    raw = bi_pop_raw(b->bi[mini], minr);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "A");
    b->bi[i]->birth_pval = pA;
    b->bi[i]->birth_fold = raw->reads/mine;
    b->bi[i]->birth_e = mine;
    b->bi[i]->birth_comp = mincomp;
    
    // Add raw to new cluster.
    bi_add_raw(b->bi[i], raw);
/*
    if(verbose) { 
      Rprintf("\nNew cluster from Raw %i in C%iR%i: ", raw->index, mini, minr);
      fold = b->bi[i]->birth_fold;
      Rprintf(" p*=%.3e, n/E(n)=%.1e\n", pA, fold);
    } */
    
    bi_assign_center(b->bi[i]);
    return i;
  }

  // No significant abundance or singleton pval
  if(verbose) { Rprintf("\nNo significant pval, no new cluster.\n"); }
  return 0;
}
