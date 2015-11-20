#include <Rcpp.h>
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
  raw->kmer = get_kmer(seq, KMER_SIZE);
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
  free(raw->kmer);
  free(raw);  
}

// The constructor for the Bi object.
Bi *bi_new(unsigned int totraw) {
//  Bi *bi = (Bi *) malloc(sizeof(Bi)); //E
  Bi *bi = new Bi;
  if (bi == NULL)  Rcpp::stop("Memory allocation failed!\n");
  bi->raw = (Raw **) malloc(RAWBUF * sizeof(Raw *)); //E
  if (bi->raw == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->maxraw = RAWBUF;
  bi->sub = (Sub **) malloc(totraw * sizeof(Sub *)); //E
  if (bi->sub == NULL)  Rcpp::stop("Memory allocation failed.");
  for(int i=0;i<totraw;i++) { bi->sub[i] = NULL; }   // Init to null pointers
  bi->lambda = (double *) malloc(totraw * sizeof(double)); //E
  if (bi->lambda == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->e = (double *) malloc(totraw * sizeof(double)); //E
  if (bi->e == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->totraw = totraw;
  bi->center = NULL;
  strcpy(bi->seq, "");
  bi->update_lambda = true;
  bi->update_e = true;
  bi->shuffle = true;
  
  bi->reads = 0;
  bi->nraw = 0;
  bi->birth_sub=NULL;
  return bi;
}

// The destructor for the Bi object
void bi_free(Bi *bi) {
  free(bi->raw);
  for(size_t index=0;index<bi->totraw;index++) { sub_free(bi->sub[index]); }
  sub_free(bi->birth_sub);
  free(bi->sub);
  free(bi->lambda);
  free(bi->e);
  delete bi;
}

// The constructor for the B object. Takes in array of Raws.
B *b_new(Raw **raws, unsigned int nraw, int score[4][4], int gap_pen, double omegaA, int band_size, bool vectorized_alignment, bool use_quals) {
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
  b->bi[0]->birth_sub = NULL;
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
    sub_free(bi->sub[pop->index]);
    bi->sub[pop->index] = NULL;
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
  unsigned int r, max_reads = 0;
  
  // Assign the raw with the most reads as the center
  bi->center = NULL;
  for(r=0;r<bi->nraw;r++) {
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
 Performs alignments and lambda of all raws to the specified Bi
*/
void b_compare(B *b, unsigned int i, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat, bool verbose) {
  unsigned int index;
  double lambda;
  Raw *raw;
  Sub *sub;
  Comparison comp;
  
  // align all raws to this sequence and compute corresponding lambda
  if(verbose) { Rprintf("C%iLU:", i); }
  for(index=0; index<b->nraw; index++) {
    raw = b->raw[index];
    // get sub object
    sub = sub_new(b->bi[i]->center, raw, b->score, b->gap_pen, use_kmers, kdist_cutoff, b->band_size, b->vectorized_alignment);
    b->nalign++;
    if(!sub) { b->nshroud++; }

    // Store sub in the cluster object Bi
    sub_free(b->bi[i]->sub[index]);
    b->bi[i]->sub[index] = sub;
    
    // Calculate lambda for that sub
    lambda = compute_lambda(raw, sub, errMat, b->use_quals);
    
    // Store lambda and set self
    b->bi[i]->lambda[index] = lambda;
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
    }
  }
  b->bi[i]->update_e = true;
  b_e_update(b);
}

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
}

/* b_shuffle:
 move each sequence to the bi that produces the highest expected
 number of that sequence. The center of a Bi cannot leave.
*/
bool b_shuffle(B *b) {
  int i, r, j, foo;
  int ibest, index;
  bool shuffled = false;
  Raw *raw;
  double e, maxe;
  
  // Make list of clusters that have shuffle flags
  std::vector<int> shuffled_bis;
  for(i=0;i<b->nclust;i++) {
    if(b->bi[i]->shuffle) {
      shuffled_bis.push_back(i);
    }
  }
  
  // Iterate over raws via clusters
  for(i=0; i<b->nclust; i++) {
    // IMPORTANT TO ITERATE BACKWARDS DUE TO BI_POP_RAW!!!!!!
    for(r=b->bi[i]->nraw-1; r>=0; r--) {
      // Find cluster with best e for this raw
      maxe = 0.0; ibest=-99;
      index = b->bi[i]->raw[r]->index;
      
      if(b->bi[i]->shuffle) { // E's for this cluster have changed, compare to all others
        for(j=0;j<b->nclust; j++) {
          e = b->bi[j]->e[index];
          if(e > maxe) {
            maxe = e;
            ibest = j;
          }
        }
      } else { // Compare just to other clusters that have changed E's
        for(foo=0;foo<shuffled_bis.size();foo++) {
          j = shuffled_bis[foo];
          e = b->bi[j]->e[index];
          if(e > maxe) {
            maxe = e;
            ibest = j;
          }
        }
      }
      
      // Check if no cluster assigned
      if(ibest == -99) {
        ibest=i;
      }
        
      // If a better cluster was found, move the raw to the new bi
      if(maxe > b->bi[i]->e[index]) {
        if(index == b->bi[i]->center->index) {  // Check if center
          if(VERBOSE) {
            Rprintf("Warning: Shuffle blocked the center of a Bi from leaving.\n");
            Rprintf("Attempted: Raw %i from C%i to C%i (%.4e (lam=%.2e,n=%i) -> %.4e (%s: lam=%.2e,n=%i))\n", \
                index, i, ibest, \
                b->bi[i]->e[index], b->bi[i]->lambda[index], b->bi[i]->reads, \
                b->bi[ibest]->e[index], b->bi[ibest]->sub[index]->key, b->bi[ibest]->lambda[index], b->bi[ibest]->reads);
          }
        } else { // Moving raw
          raw = bi_pop_raw(b->bi[i], r);
          bi_add_raw(b->bi[ibest], raw);
          shuffled = true;
        }  
      }
    } //End loop(s) over raws (r).
  } // End loop over clusters (i).
  
  for(i=0; i<b->nclust; i++) { b->bi[i]->shuffle = false; }
  return shuffled;
}

/* b_shuffle2:
 move each sequence to the bi that produces the highest expected
 number of that sequence. The center of a Bi cannot leave.
*/
bool b_shuffle2(B *b) {
  unsigned int i, j, index;
  bool shuffled = false;
  Raw *raw, *pop;
  
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

  return shuffled;
}

/* bi_free_absent_subs:
   Frees subs of raws that are not currently in this cluster.
   */
void bi_free_absent_subs(Bi *bi, unsigned int nraw) {
  unsigned int r, index;
  bool *keep = (bool *) malloc(nraw * sizeof(bool)); //E
  if (keep == NULL)  Rcpp::stop("Memory allocation failed.");
  for(index=0;index<nraw;index++) { keep[index] = false; }
  
  for(r=0;r<bi->nraw;r++) {
    keep[bi->raw[r]->index] = true;
  }
  for(index=0;index<nraw;index++) {
    if(!keep[index]) {
      sub_free(bi->sub[index]);
      bi->sub[index] = NULL;
    }
  }
  free(keep);
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
      raw->p = get_pA(raw, b->bi[i]->sub[raw->index], b->bi[i]->lambda[raw->index], b->bi[i]);      
    } // for(r=0;r<b->bi[i]->nraw;r++)
  } // for(i=0;i<b->nclust;i++)
}

/* b_bud:
 Finds the minimum p-value. If significant, creates a new cluster and moves the
 raws from the raw with the minimum p-value to the new cluster.
 Returns index of new cluster, or 0 if no new cluster added.
*/

int b_bud(B *b, double min_fold, int min_hamming, bool verbose) {
  int i, r;
  int mini, minr, minreads;
  double minp = 1.0;
  double pA=1.0;
  double fold;
  Raw *raw;
  Sub *sub, *birth_sub;

  // Find i, r indices and value of minimum pval.
  mini=-999; minr=-999; minreads=0; minp=1.0;
  for(i=0;i<b->nclust;i++) {
    for(r=0; r<b->bi[i]->nraw; r++) {
      raw = b->bi[i]->raw[r];
      sub = b->bi[i]->sub[raw->index];
      
      if(!sub) {  // This shouldn't happen
        Rprintf("Warning: Raw has null sub in b_bud.\n");
        continue; 
      }

      // Calculate the fold over-abundnace and the hamming distance to this raw
      if(sub->nsubs >= min_hamming) { // Only those passing the hamming/fold screens can be budded
        if(min_fold <= 1 || ((double) raw->reads) >= min_fold * b->bi[i]->e[raw->index]) {  
          if((raw->p < minp) ||
            ((raw->p == minp && raw->reads > minreads))) { // Most significant
            mini = i; minr = r;
            minp = raw->p;
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
    birth_sub = sub_copy(b->bi[mini]->sub[b->bi[mini]->raw[minr]->index]); // do this before bi_pop_raw destructs the sub
    raw = bi_pop_raw(b->bi[mini], minr);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "A");
    b->bi[i]->birth_pval = pA;
    b->bi[i]->birth_fold = raw->reads/b->bi[mini]->e[raw->index];
    b->bi[i]->birth_e = b->bi[mini]->e[raw->index];
    b->bi[i]->birth_sub = birth_sub;
    
    // Add raw to new cluster.
    bi_add_raw(b->bi[i], raw);

    if(verbose) { 
      double qave = 0.0;
      if(raw->qual) {
        for(int s=0;s<birth_sub->nsubs;s++) {
          qave += raw->qual[birth_sub->pos[s]];
        }
        qave = qave/((double)birth_sub->nsubs);
      }
      Rprintf("\nNew cluster from Raw %i in C%iR%i: ", raw->index, mini, minr);
      fold = b->bi[i]->birth_fold;
      Rprintf(" p*=%.3e, n/E(n)=%.1e (%.1e fold per sub)\n", pA, fold, fold/birth_sub->nsubs);
      Rprintf("Reads: %i, E: %.2e, Nsubs: %i, Ave Qsub:%.1f\n", raw->reads, b->bi[mini]->e[raw->index], birth_sub->nsubs, qave);
      Rprintf("%s\n", ntstr(birth_sub->key));
    }
    
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
*/
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
}

/* B_make_consensus:
  Makes consensus sequences for all Bi.
*/
void b_make_consensus(B *b) {
  for (int i=0; i<b->nclust; i++) {
    bi_make_consensus(b->bi[i], b->use_quals);
  }
}
