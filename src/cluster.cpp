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
    raw->qual = (float *) malloc(raw->length * sizeof(float)); //E
    if (raw->qual == NULL)  Rcpp::stop("Memory allocation failed.");
    for(size_t i=0;i<raw->length;i++) { raw->qual[i] = (float) qual[i]; }
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
  bi->update_lambda = true;
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

// Initializion B object by placing all Raws into one Bi
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
// Flags update_lambda.
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
  // Assign bi->seq and flag
  if(bi->center) { strcpy(bi->seq, bi->center->seq); }
  bi->update_lambda = true;
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
  Raw *center = b->bi[i]->center;
  uint8_t *center_kmer8 = center->kmer8; // Store original kmer8 pointer
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
      b->bi[i]->comp_index.insert(std::make_pair(index, cind++));
    }
    sub_free(sub);
  }
}

/***********************************************
 * SUPPLANTED BY RcppParallel IMPLEMENTATION
 ***********************************************

typedef struct {
  unsigned int start;
  unsigned int end;
  B *b;
  unsigned int i;
  bool use_kmers;
  double kdist_cutoff;
  Comparison *comps;
  unsigned int ncol;
  double *err_mat;
} tc_input;

void *t_compare(void *tc_in) {
  tc_input *inp = (tc_input *) tc_in;
  unsigned int index, j;
  Raw *raw;
  Sub *sub;
  B *b = inp->b;
  unsigned int i = inp->i;
  Comparison *comps = inp->comps;

  for(index=inp->start, j=0; index<inp->end; index++, j++) {
    // Make sub object
    raw = b->raw[index];
    sub = sub_new(b->bi[i]->center, raw, b->score, b->gap_pen, b->homo_gap_pen, inp->use_kmers, inp->kdist_cutoff, b->band_size, b->vectorized_alignment, SSE);

    // Make comparison object
    comps[j].i = i;
    comps[j].index = index;
    comps[j].lambda = compute_lambda_ts(raw, sub, inp->ncol, inp->err_mat, b->use_quals);
    if(sub) {
      comps[j].hamming = sub->nsubs;
    } else {
      comps[j].hamming = -1;
    }

    // Free sub
    sub_free(sub);
  }
  
  return(NULL);
}

// compare_threaded:
// Performs alignments and computes lambda for all raws to the specified Bi
// Stores only those that can possibly be recruited to this Bi
// MULTITHREADED

void b_compare_threaded(B *b, unsigned int i, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat, unsigned int nthreads, bool verbose, int SSE) {
  unsigned int index, cind, thr, row, col, ncol;
  double lambda;
  Raw *raw;
  Comparison comp;
  pthread_t threads[nthreads-1];
  tc_input inp[nthreads];
  Comparison *comps = (Comparison *) malloc(sizeof(Comparison) * b->nraw);
  if(comps==NULL) Rcpp::stop("Memory allocation failed.");
  
  // Make thread-safe simple array for the error rate matrix
  double *err_mat = (double *) malloc(sizeof(double) * errMat.ncol() * errMat.nrow());
  if(err_mat==NULL) Rcpp::stop("Memory allocation failed.");
  ncol = errMat.ncol();
  if(errMat.nrow() != 16) { Rcpp::stop("Error matrix doesn't have 16 rows."); }
  for(row=0;row<errMat.nrow();row++) {
    for(col=0;col<errMat.ncol();col++) {
      err_mat[row*ncol + col] = errMat(row, col);
    }
  }
  
  // Threaded loop to perform all comparisons
  for(thr=0;thr<nthreads;thr++) {
    inp[thr].start = (b->nraw*thr)/nthreads;
    inp[thr].end = (b->nraw*(thr+1))/nthreads;
    inp[thr].b = b;
    inp[thr].i = i;
    inp[thr].use_kmers = use_kmers;
    inp[thr].kdist_cutoff = kdist_cutoff;
    inp[thr].comps = &comps[inp[thr].start];
    inp[thr].ncol = ncol;
    inp[thr].err_mat = err_mat;
    if(thr == nthreads-1) { break; } // run in this thread
    pthread_create(&threads[thr], NULL, t_compare, (void *) &inp[thr]);
  }
  t_compare((void *) &inp[nthreads-1]);
  
  // Join the threads into a final vector of comparisons (and free stuff)  
  for(thr=0; thr<(nthreads-1); thr++) {
    pthread_join(threads[thr], NULL);
  }
  
  // Selectively store
  for(index=0, cind=0; index<b->nraw; index++) {
    b->nalign++; ///t
    raw = b->raw[index];
    comp = comps[index];
    lambda = comp.lambda;
    
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
      b->bi[i]->comp_index.insert(std::make_pair(index, cind++));
    }
  }
  free(err_mat);
  free(comps);
}
***********************************************/
 
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
      b->bi[i]->comp_index.insert(std::make_pair(index, cind++));
    }
  }
  free(err_mat);
  free(comps);
}


/*
void b_e_update(B *b) {
  unsigned int i, index;

  for(i=0; i < b->nclust; i++) {
    if(b->bi[i]->update_e) {
      for(index=0; index < b->nraw; index++) {
        b->bi[i]->e[index] = b->bi[i]->lambda[index]*b->bi[i]->reads;
      }
      b->bi[i]->shuffle = true;
    }
    b->bi[i]->update_e = false;
  }
} */

/* b_shuffle2:
 move each sequence to the bi that produces the highest expected
 number of that sequence. The center of a Bi cannot leave.
*/
bool b_shuffle2(B *b) {
  unsigned int i, j, index;
  bool shuffled = false;
  Raw *raw;
  
  double *emax = (double *) malloc(b->nraw * sizeof(double));
  unsigned int *imax = (unsigned int *) malloc(b->nraw * sizeof(unsigned int));
  if(emax==NULL || imax==NULL) Rcpp::stop("Memory allocation failed.");

  // Initialize emax/imax off of cluster 0
  // Comparisons to all raws exist in cluster 0, in index order
  for(index=0;index<b->nraw;index++) {
    emax[index] = b->bi[0]->comp[index].lambda * b->bi[0]->reads;
    imax[index] = 0;
  }
  
  // Iterate over remaining comparisons, find best E/i for each raw
  for(i=1;i<b->nclust;i++) {
    for(j=0;j<b->bi[i]->comp.size();j++) {
      if(b->bi[i]->comp[j].lambda * b->bi[i]->reads > emax[b->bi[i]->comp[j].index]) { // better E
        emax[b->bi[i]->comp[j].index] = b->bi[i]->comp[j].lambda * b->bi[i]->reads;
        imax[b->bi[i]->comp[j].index] = b->bi[i]->comp[j].i;
      }
    }
  }
  
  // Iterate over raws, if best i different than current, move
  for(i=0;i<b->nclust;i++) {
    // IMPORTANT TO ITERATE BACKWARDS DUE TO BI_POP_RAW!!!!!!
    for(int r=b->bi[i]->nraw-1; r>=0; r--) {
      raw = b->bi[i]->raw[r];
      // If a better cluster was found, move the raw to the new bi
      if(imax[raw->index] != i) {
        if(raw->index == b->bi[i]->center->index) {  // Check if center
          if(VERBOSE) { Rprintf("Warning: Shuffle blocked the center of a Bi from leaving."); }
          continue;
        }
        // Moving raw
        bi_pop_raw(b->bi[i], r);
        bi_add_raw(b->bi[imax[raw->index]], raw);
        shuffled = true;  
      }  
    } // for(r=0;r<b->bi[i]->nraw;r++)
  }

  free(emax);
  free(imax);
  
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
    for(r=0;r<b->bi[i]->nraw;r++) {
      raw = b->bi[i]->raw[r];
      raw->p = get_pA(raw, b->bi[i]);
    } // for(r=0;r<b->bi[i]->nraw;r++)
  } // for(i=0;i<b->nclust;i++)
}

/* b_bud:
 Finds the minimum p-value. If significant, creates a new cluster and moves the
 raws from the raw with the minimum p-value to the new cluster.
 Returns index of new cluster, or 0 if no new cluster added.
*/

int b_bud(B *b, double min_fold, int min_hamming, bool verbose) {
  int i, r, ci;
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
      ci = b->bi[i]->comp_index[raw->index];
      hamming = b->bi[i]->comp[ci].hamming;
      lambda = b->bi[i]->comp[ci].lambda;

      // Calculate the fold over-abundance and the hamming distance to this raw
      if(hamming >= min_hamming) { // Only those passing the hamming/fold screens can be budded
        if(min_fold <= 1 || ((double) raw->reads) >= min_fold * lambda * b->bi[i]->reads) {  
          if((raw->p < minp) ||
            ((raw->p == minp && raw->reads > minreads))) { // Most significant
            mini = i; minr = r;
            minp = raw->p;
            mine = lambda * b->bi[i]->reads;
            mincomp = b->bi[i]->comp[ci];
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

/* Bi_make_consensus:
  Uses the alignments to the cluster center to construct a consensus sequence, which
  is then assigned to bi->seq

void bi_make_consensus(Bi *bi, bool use_quals) {
  int i,r,s,nti0,nti1,pos;
  double counts[4][SEQLEN] = {{0.0}};
  Sub *sub = NULL;
  Raw *raw;
  
  if(!bi->center) { return; }
  int len = bi->center->length;
  if(len >= SEQLEN) { return; }
  
  for(r=0;r<bi->nraw;r++) {
    raw = bi->raw[r];
    sub = bi->sub[raw->index];
    // Assign all counts to center sequence
    for(i=0;i<len;i++) {
      nti0 = ((int) bi->center->seq[i]) - 1;
      if(use_quals) {
        counts[nti0][i] += (raw->reads * (double) raw->qual[i]);  // NEED TO EXCEPT THIS IF NOT USING QUALS!!
      } else {
        counts[nti0][i] += raw->reads;
      }
    }
    // Move counts at all subs to proper nt
    if(sub) {
      for(s=0;s<sub->nsubs;s++) {
        nti0 = ((int) sub->nt0[s]) - 1;
        nti1 = ((int) sub->nt1[s]) - 1;
        pos = sub->pos[s];
        if(use_quals) {
          counts[nti0][pos] -= (raw->reads * (double) raw->qual[pos]);
          counts[nti1][pos] += (raw->reads * (double) raw->qual[pos]);
        } else {
          counts[nti0][pos] -= raw->reads;
          counts[nti1][pos] += raw->reads;
        }
      }
    }
  } // for(r=0;r<bi->nraw;r++)
  
  // Use counts to write consensus
  strcpy(bi->seq, bi->center->seq);
  for(i=0;i<len;i++) {
    if(counts[0][i] > counts[1][i] && counts[0][i] > counts[2][i] && counts[0][i] > counts[3][i]) {
      bi->seq[i] = 1;
    } 
    else if(counts[1][i] > counts[0][i] && counts[1][i] > counts[2][i] && counts[1][i] > counts[3][i]) {
      bi->seq[i] = 2;
    }
    else if(counts[2][i] > counts[0][i] && counts[2][i] > counts[1][i] && counts[2][i] > counts[3][i]) {
      bi->seq[i] = 3;
    }
    else if(counts[3][i] > counts[0][i] && counts[3][i] > counts[1][i] && counts[3][i] > counts[2][i]) {
      bi->seq[i] = 4;
    } else { // TIES! These are not properly dealt with currently!
      bi->seq[i] = bi->center->seq[i];
    }
  }
  bi->seq[len] = '\0';
} */

/* B_make_consensus:
  Makes consensus sequences for all Bi.

void b_make_consensus(B *b) {
  for (int i=0; i<b->nclust; i++) {
    bi_make_consensus(b->bi[i], b->use_quals);
  }
}
*/
