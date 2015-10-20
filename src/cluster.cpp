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
#define FAMBUF 50
#define CLUSTBUF 50

/* private function declarations */
Raw *raw_new(char *seq, int reads);
Fam *fam_new();
Bi *bi_new(int totraw);
void fam_free(Fam *fam);
void bi_free(Bi *bi);

int b_add_bi(B *b, Bi *bi);
Bi *bi_new(int totraw);
Raw *bi_pop_raw(Bi *bi, int f, int r);
void bi_shove_raw(Bi *bi, Raw *raw);

int fam_add_raw(Fam *fam, Raw *raw);
int bi_add_fam(Bi *bi, Fam *fam);
int b_add_bi(B *b, Bi *bi);
void bi_shove_raw(Bi *bi, Raw *raw);
void bi_add_raw(Bi *bi, Raw *raw);   // NOT WORKING, MAY BE DELETED
Raw *fam_pop_raw(Fam *fam, int r);
Raw *bi_pop_raw(Bi *bi, int f, int r);
Fam *bi_pop_fam(Bi *bi, int f);

void bi_census(Bi *bi);
void bi_center_update(Bi *bi);
void fam_center_update(Fam *fam);
void bi_fam_update(Bi *bi, int score[4][4], int gap_pen, int band_size, bool use_quals, bool verbose);
void bi_make_consensus(Bi *bi, bool use_quals);

/*
 raw_new:
 The constructor for the Raw object.
 */
Raw *raw_new(char *seq, int reads) {
  Raw *raw = (Raw *) malloc(sizeof(Raw)); //E
  if (raw == NULL)  Rcpp::stop("Memory allocation failed.");
  raw->seq = (char *) malloc(strlen(seq)+1); //E
  if (raw->seq == NULL)  Rcpp::stop("Memory allocation failed.");
  strcpy(raw->seq, seq);
  raw->length = strlen(seq);
  raw->qual = NULL;
  raw->kmer = get_kmer(seq, KMER_SIZE);
  raw->reads = reads;
  return raw;
}

// raw_qual_new: A constructor for Raw objects with quals
Raw *raw_qual_new(char *seq, double *qual, int reads) {
  Raw *raw = raw_new(seq, reads);
  if(qual) {
    raw->qual = (float *) malloc(raw->length * sizeof(float)); //E
    if (raw->qual == NULL)  Rcpp::stop("Memory allocation failed.");
    for(int i=0;i<raw->length;i++) { raw->qual[i] = qual[i]; }
  } else {
    Rcpp::stop("Error: NULL qual provided to raw_qual_new constructor.");
  }
  return raw;
}

void raw_free(Raw *raw) {
  free(raw->seq);
  if(raw->qual) { free(raw->qual); }
  free(raw->kmer);
  free(raw);  
}
/*
 fam_new:
 The constructor for the Fam object.
 */
Fam *fam_new() {
  Fam *fam = (Fam *) malloc(sizeof(Fam)); //E
  if (fam == NULL)  Rcpp::stop("Memory allocation failed.");
  fam->seq = (char *) malloc(SEQLEN); //E
  if (fam->seq == NULL)  Rcpp::stop("Memory allocation failed.");
  fam->raw = (Raw **) malloc(RAWBUF * sizeof(Raw *)); //E
  if (fam->raw == NULL)  Rcpp::stop("Memory allocation failed.");
  fam->maxraw = RAWBUF;
  fam->nraw = 0;
  fam->reads = 0;
  return fam;
}

/*
 fam_free:
 The destructor for the Fam object.
 */
void fam_free(Fam *fam) {
  free(fam->seq);
  free(fam->raw);
  free(fam);
}

/*
 fam_add_raw:
 Add a raw to a fam object.
 Update reads/nraw.
 Return index to raw in fam.
 */
int fam_add_raw(Fam *fam, Raw *raw) {
  if(fam->nraw >= fam->maxraw) {    // Extend Raw* buffer
    fam->raw = (Raw **) realloc(fam->raw, (fam->maxraw+RAWBUF) * sizeof(Raw *)); //E
    if (fam->raw == NULL)  Rcpp::stop("Memory allocation failed.");
    fam->maxraw+=RAWBUF;
  }

  fam->raw[fam->nraw] = raw;
  fam->reads += raw->reads;
  return(fam->nraw++);
}

/*
 fam_pop_raw:
 Removes a raw from a fam object.
 Update fam reads/nraw.
 Returns pointer to that raw.
 THIS REORDERS THE FAMS LIST OF RAWS. CAN BREAK INDICES.
 THIS DOES NOT UPDATE READS OR NRAW IN THE PARENT BI.
 */
Raw *fam_pop_raw(Fam *fam, int r) {
  Raw *pop;
  
  if(r<fam->nraw) {
    pop = fam->raw[r];
    fam->raw[r] = fam->raw[fam->nraw-1];
    fam->raw[fam->nraw-1] = NULL;
    fam->nraw--;
    fam->reads -= pop->reads;
  } else {  // Trying to pop an out of range raw
    Rprintf("fam_pop_raw: Not enough raws %i (%i)\n", r, fam->nraw);
    pop = NULL;
  }
  return pop;
}

/*
 bi_pop_raw:
 Removes a raw from the child object.
 Update bi reads/nraw.
 Returns pointer to that raw.
 SEE NOTES ON FAM_POP_RAW
 */
Raw *bi_pop_raw(Bi *bi, int f, int r) {
  Raw *pop;
  
  if(f<bi->nfam) {
    pop = fam_pop_raw(bi->fam[f], r);
    if(pop != NULL) {
      bi->nraw--;
      bi->reads -= pop->reads;
      bi->update_fam = true;
      bi->update_e = true;
      sub_free(bi->sub[pop->index]);  ///! Free the sub
      bi->sub[pop->index] = NULL;     ///! Free the sub
    }
  } else {  // Trying to pop an out of range raw
    Rprintf("bi_pop_raw: not enough fams %i (%i)\n", f, bi->nfam);
    pop = NULL;
  }
  return pop;
}

/*
 bi_pop_fam:
 Removes a fam from cluster.
 Update bi reads/nraw.
 Returns pointer to the removed fam.
 */
Fam *bi_pop_fam(Bi *bi, int f) {
  Fam *pop;
  
  if(f<bi->nfam) {
    pop = bi->fam[f];
    bi->fam[f] = bi->fam[bi->nfam-1];
    bi->fam[bi->nfam-1] = NULL;
    bi->nfam--;
    bi->nraw -= pop->nraw;
    for(int r=0;r<pop->nraw;r++) {  // Doesn't rely on fam reads
      bi->reads -= pop->raw[r]->reads;
      sub_free(bi->sub[pop->raw[r]->index]);   ///! Free the sub
      bi->sub[pop->raw[r]->index] = NULL;      ///! Free the sub
      bi->update_e = true;
    }
  } else {  // Trying to pop an out of range fam
    Rprintf("bi_pop_fam: not enough fams %i (%i)\n", f, bi->nfam);
    pop = NULL;
  }
  return pop;
}

/*
 bi_shove_raw:
 Hackily just adds raw to fam[0].
 fam_update and e_update are needed on this cluster afterwards.
 */
void bi_shove_raw(Bi *bi, Raw *raw) {
  if(bi->nfam == 0) {
    bi_add_fam(bi, fam_new());
  }
  fam_add_raw(bi->fam[0], raw);
  bi->nraw++;
  bi->reads += raw->reads;
  bi->update_fam = true;
  bi->update_e = true;
}

/*
 bi_add_raw:
 Adds raw to appropriate fam. Makes new fam if necessary.
 ISSUE WITH ALIGN, NEED ERR. UNUSED FOR NOW.
 UNUSED UNUSED UNUSED UNUSED UNUSED UNUSED
 */
void bi_add_raw(Bi *bi, Raw *raw) {
//  Rprintf("bi_add_raw\n");
  Sub *sub;
  int foo, result;
  char buf[10];  // holds sprintf("%i", fam_index)
  
  sub = bi->sub[raw->index];
  result = sm_exists(bi->sm, sub->key);
  
  if (result == 0) {                  // Handle value not found
    foo = bi_add_fam(bi, fam_new());
    fam_add_raw(bi->fam[foo], raw);
    bi->fam[foo]->sub = sub;                   // Sub set on new fam formation.
    sprintf(buf, "%i", foo);
    sm_put(bi->sm, sub->key, buf);
  } else {                            // Joining existing family
    sm_get(bi->sm, sub->key, buf, sizeof(buf));
    foo = atoi(buf);
    fam_add_raw(bi->fam[foo], raw);
  }
  
  bi->nraw++;
  bi->reads += raw->reads;
}

/* bi_new:
 The constructor for the Bi object.
 */
Bi *bi_new(int totraw) {
  Bi *bi = (Bi *) malloc(sizeof(Bi)); //E
  if (bi == NULL)  Rcpp::stop("Memory allocation failed!\n");
  strcpy(bi->seq, "");
  bi->fam = (Fam **) malloc(FAMBUF * sizeof(Fam *)); //E
  if (bi->fam == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->maxfam = FAMBUF;
  bi->sub = (Sub **) malloc(totraw * sizeof(Sub *)); //E
  if (bi->sub == NULL)  Rcpp::stop("Memory allocation failed.");
  for(int i=0;i<totraw;i++) { bi->sub[i] = NULL; }   // Init to null pointers
  bi->lambda = (double *) malloc(totraw * sizeof(double)); //E
  if (bi->lambda == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->e = (double *) malloc(totraw * sizeof(double)); //E
  if (bi->e == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->totraw = totraw;
  bi->sm = sm_new(MIN_BUCKETS); //E
  if (bi->sm == NULL)  Rcpp::stop("Memory allocation failed.");
  bi->center = NULL;
  bi->update_lambda = true;
  bi->update_fam = true;
  bi->update_e = true;
  bi->shuffle = true;
  
  bi->nfam = 0;
  bi->reads = 0;
  bi->nraw = 0;
  bi->birth_sub=NULL;
  return bi;
}

/* bi_free:
 Destructs bi objects.
*/
void bi_free(Bi *bi) {
  for(int f=0;f<bi->nfam;f++) { fam_free(bi->fam[f]); }
  free(bi->fam);
  for(int i=0;i<bi->totraw;i++) { sub_free(bi->sub[i]); }
  sub_free(bi->birth_sub);
  free(bi->sub);
  free(bi->lambda);
  free(bi->e);
  sm_delete(bi->sm);
  free(bi);
}

/* bi_add_fam:
 Adds a fam to an existing Bi. Returns its index.
 */
int bi_add_fam(Bi *bi, Fam *fam) {
  if(bi->nfam >= bi->maxfam) {    // Extend Fam* buffer
    bi->fam = (Fam **) realloc(bi->fam, (bi->maxfam+FAMBUF) * sizeof(Fam *)); //E
    if (bi->fam == NULL)  Rcpp::stop("Memory allocation failed.");
    bi->maxfam+=FAMBUF;
  }

  bi->fam[bi->nfam] = fam;
  bi->reads += fam->reads;
  return(bi->nfam++);
}

void bi_census(Bi *bi) {
  int f, r;
  int reads=0, nraw=0;
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      reads += bi->fam[f]->raw[r]->reads;
      nraw++;
    }
  }
  if(reads != bi->reads) {
    bi->reads = reads;
    bi->update_e = true;
  }
  bi->nraw = nraw;
}

/* B_new:
 The constructor for the B object. Takes in a Uniques object.
 Places all sequences into the same family within one cluster.
*/
B *b_new(Raw **raws, int nraw, int score[4][4], int gap_pen, double omegaA, bool use_singletons, double omegaS, int band_size, bool use_quals) {
  int i, j, nti;
  size_t index;

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
  b->use_singletons = use_singletons;
  b->omegaS = omegaS;
  b->band_size = band_size;
  b->use_quals = use_quals;
  
  // Copy the score matrix
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      b->score[i][j] = score[i][j];
    }
  }

  double tot_nnt[] = {0.0,0.0,0.0,0.0};
  b->raw = raws;
  for (index = 0; index < b->nraw; index++) {
    for(i=0;i<b->raw[index]->length;i++) {
      nti = ((int) b->raw[index]->seq[i]) - 1;
      if(nti == 0 || nti == 1 || nti ==2 || nti == 3) {
        tot_nnt[nti]++;
      }
    }
    b->raw[index]->index = index;
    b->reads += b->raw[index]->reads;
  }

  // Create the (for now) solitary lookup table from lambda -> pS
  // Use the average sequence composition
  // Pval lookup depends on sequence composition when the error matrix
  //   is not uniform, but this dependence is not super-strong normally.
  // However, we very much need the sequences to be all close to the same
  //   length for this to be valid!!!!
  if(b->use_singletons) {
    b_make_pS_lookup(b);
  }

  // Initialize with one cluster/one-family containing all the raws.
  b_init(b);
  return b;
}

/* b_init:
 Initialized b with all raws in one cluster, with center, flagged to update lambda and fam.
*/
void b_init(B *b) {
  // Destruct existing clusters/fams
  for(int i=0; i<b->nclust; i++) {
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
  for (size_t index=0; index<b->nraw; index++) {
    bi_shove_raw(b->bi[0], b->raw[index]);
  }

  bi_census(b->bi[0]);
  bi_center_update(b->bi[0]); // Makes cluster center sequence
}

/* b_free:
  Destruct the B object.
*/
void b_free(B *b) {
  for(int i=0;i<b->nclust;i++) { bi_free(b->bi[i]); }
  free(b->bi);
  free(b);
}

/* b_add_bi:
 Adds a new cluster Bi to the clustering B. Returns its index.
 */
int b_add_bi(B *b, Bi *bi) {
  if(b->nclust >= b->maxclust) {    // Extend Bi* buffer
    b->bi = (Bi **) realloc(b->bi, (b->maxclust+CLUSTBUF) * sizeof(Bi *)); //E
    if (b->bi == NULL)  Rcpp::stop("Memory allocation failed.");
    b->maxclust+=CLUSTBUF;
  }
  b->bi[b->nclust] = bi;
  bi->i = b->nclust;
  return(b->nclust++);
}

/*
 lambda_update:
 updates the alignments and lambda of all raws to Bi with
 changed center sequences.
*/
void b_lambda_update(B *b, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat, bool verbose) {
  int i;
  size_t index;
  double lambda;
  Sub *sub;
  
  for (i = 0; i < b->nclust; i++) {
    if(b->bi[i]->update_lambda) {   // center sequence for Bi[i] has changed
      // update alignments and lambda of all raws to this sequence
      if(verbose) { Rprintf("C%iLU:", i); }
      for(index=0; index<b->nraw; index++) {
        // get sub object
        sub = sub_new(b->bi[i]->center, b->raw[index], b->score, b->gap_pen, use_kmers, kdist_cutoff, b->band_size);
        b->nalign++;
        if(!sub) { b->nshroud++; }
  
        // Store sub in the cluster object Bi
        sub_free(b->bi[i]->sub[index]);
        b->bi[i]->sub[index] = sub;
        
        // Calculate lambda for that sub
        /// bi->self = get_self(bi->seq, err);
        /// lambda = compute_lambda(sub, b->bi[i]->self, b->err, b->use_quals);  // SELF USED TO BE SET IN BI_CONSENSUS_UPDATE
        lambda = compute_lambda3(b->raw[index], sub, errMat, b->use_quals);
        
        // Store lambda and set self
        b->bi[i]->lambda[index] = lambda;
        if(index == b->bi[i]->center->index) { b->bi[i]->self = lambda; }
        if(index == TARGET_RAW) { Rprintf("lam(TARG)=%.2e; ", b->bi[i]->lambda[index]); }
      }
      b->bi[i]->update_lambda = false;
      b->bi[i]->update_fam = true;
      b->bi[i]->update_e = true;
    } // if(b->bi[i]->update_lambda)
  }
  b_e_update(b);
}

/* bi_fam_update(Bi *bi, ...):
   Creates and allocates new fams based on the sub between each raw and the
   cluster center.
  Currently completely destructive of old fams.
   */
void bi_fam_update(Bi *bi, int score[4][4], int gap_pen, int band_size, bool use_quals, bool verbose) {
  int f, r, result, r_c;
  Sub *sub;
  char buf[10];
  
  sm_delete(bi->sm);
  bi->sm = sm_new(MIN_BUCKETS + (int) (BUCKET_SCALE * bi->nraw));  // n_buckets scales with # of raws
  // More buckets = more memory, but less collision (and therefore less costly strcmp w/in buckets).
  bi_census(bi);
  
  // Make list of pointers to the raws
  Raw **raws = (Raw **) malloc(bi->nraw * sizeof(Raw *)); //E
  if (raws == NULL)  Rcpp::stop("Memory allocation failed.");
  r_c=0;
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      raws[r_c] = bi->fam[f]->raw[r];
      r_c++;
    }
  }
  
  // Destruct old fams
  for(f=0;f<bi->nfam;f++) { fam_free(bi->fam[f]); }
  bi->nfam = 0;
  
  // Guarantee that non-NULL sub objects exist for all raws to the cluster center
  // NULL subs can occur (rarely) with kmers if cluster center changes (in shuffle) and 
  //    now exceeds kmerdist to another raw in the cluster
  for(r_c=0;r_c<bi->nraw;r_c++) {
    if(!(bi->sub[raws[r_c]->index])) { // Protect from and replace null subs
      bi->sub[raws[r_c]->index] = sub_new(bi->center, raws[r_c], score, gap_pen, false, 1., band_size);
      if(verbose) Rprintf("F");
      // DOESN'T UPDATE LAMBDA HERE, THAT ONLY HAPPENS IN COMPUTE_LAMBDA
    }
  }
  
  // Construct hash from raw->sub->key to new fam index.
  // Make fams, that contain all raws with the same substitution pattern.
  for(r_c=0;r_c<bi->nraw;r_c++) {
    // Make hash key from substitutions
    sub = bi->sub[raws[r_c]->index];
    result = sm_exists(bi->sm, sub->key);

    // Place raw in fams.
    if (result == 0) {                  // Handle value not found
      f = bi_add_fam(bi, fam_new());
      fam_add_raw(bi->fam[f], raws[r_c]);
      bi->fam[f]->sub = sub;                   // Sub set on new fam formation.
      sprintf(buf, "%i", f); // strmap only takes strings as values
      sm_put(bi->sm, sub->key, buf);
    } else {                            // Joining existing family
      sm_get(bi->sm, sub->key, buf, sizeof(buf));
      f = atoi(buf);
      fam_add_raw(bi->fam[f], raws[r_c]);
    }
  }
  
  // center_update/align the fams
  for(f=0;f<bi->nfam;f++) {
    fam_center_update(bi->fam[f]);
    sub = bi->sub[bi->fam[f]->center->index];

    if(!sub) { // Protect from null subs, but this should never arise...
      Rcpp::stop("Error: bi_fam_update hit a null sub. THIS SHOULDNT HAPPEN!!!!!\n");
    }
    
    bi->fam[f]->sub = sub;
    bi->fam[f]->lambda = bi->lambda[bi->fam[f]->center->index];
  }
  free(raws);
  bi->update_fam = false;
}

void b_fam_update(B *b, bool verbose) {
  for (int i=0; i<b->nclust; i++) {
    if(b->bi[i]->update_fam) {  // center has changed or a raw has been shoved EITHER DIRECTION breaking the fam structure
      if(verbose) { Rprintf("C%iFU:", i); }
      bi_fam_update(b->bi[i], b->score, b->gap_pen, b->band_size, b->use_quals, verbose);
    }
  }
}

/* bi_free_absent_subs(Bi *bi, int nraw):
   Frees subs of raws that are not currently in this cluster.
   */
void bi_free_absent_subs(Bi *bi, int nraw) {
  int f, r, index;
  bool *keep = (bool *) malloc(nraw * sizeof(bool)); //E
  if (keep == NULL)  Rcpp::stop("Memory allocation failed.");
  for(index=0;index<nraw;index++) { keep[index] = false; }
  
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      keep[bi->fam[f]->raw[r]->index] = true;
    }
  }
  for(index=0;index<nraw;index++) {
    if(!keep[index]) {
      sub_free(bi->sub[index]);
      bi->sub[index] = NULL;
    }
  }
  free(keep);
}

void b_e_update(B *b) {
  int i;
  size_t index;

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
  int i, f, r, j, foo;
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
  
  // Iterate over raws via clusters/fams
  for(i=0; i<b->nclust; i++) {
    for(f=0; f<b->bi[i]->nfam; f++) {
      // IMPORTANT TO ITERATE BACKWARDS DUE TO FAM_POP_RAW!!!!!!
      for(r=b->bi[i]->fam[f]->nraw-1; r>=0; r--) {
        // Find cluster with best e for this raw
        maxe = 0.0; ibest=-99;
        index = b->bi[i]->fam[f]->raw[r]->index;
        
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
            raw = bi_pop_raw(b->bi[i], f, r);
            bi_shove_raw(b->bi[ibest], raw);
            shuffled = true;
          }  
        }
      } //End loop(s) over raws (r).
    } // End loop over fams (f).
  } // End loop over clusters (i).
  
  for(i=0; i<b->nclust; i++) { b->bi[i]->shuffle = false; }
  return shuffled;
}

/* b_shuffle:
 Move sequences to the newest bi if its E-value is higher than for the
 cluster they are currently in. Cluster centers cannot move.
*/
bool b_shuffle_oneway(B *b) {
  int i, f, r, inew;
  int ibest, index;
  bool shuffled = false;
  Raw *raw;
  double newe;
  
  inew = b->nclust-1;
  // Iterate over raws via clusters/fams
  for(i=0; i<inew; i++) {
    for(f=0; f<b->bi[i]->nfam; f++) {
      // IMPORTANT TO ITERATE BACKWARDS DUE TO FAM_POP_RAW!!!!!!
      for(r=b->bi[i]->fam[f]->nraw-1; r>=0; r--) {
        // Is e better in the new cluster?
        index = b->bi[i]->fam[f]->raw[r]->index;
        newe = b->bi[inew]->e[index];
                
        // If new cluster is better, move the raw to the new bi
        if(newe > b->bi[i]->e[index]) {
          if(index == b->bi[i]->center->index) {  // Check if center
            if(VERBOSE) {
              Rprintf("Warning: Shuffle blocked the center of a Bi from leaving.\n");
              Rprintf("Attempted: Raw %i from C%i to C%i (%.4e (lam=%.2e,n=%i) -> %.4e (%s: lam=%.2e,n=%i))\n", \
                  index, i, ibest, \
                  b->bi[i]->e[index], b->bi[i]->lambda[index], b->bi[i]->reads, \
                  b->bi[ibest]->e[index], b->bi[ibest]->sub[index]->key, b->bi[ibest]->lambda[index], b->bi[ibest]->reads);
            }
          } else { // Moving raw
            raw = bi_pop_raw(b->bi[i], f, r);
            bi_shove_raw(b->bi[inew], raw);
            shuffled = true;
          }  
        }
      } //End loop(s) over raws (r).
    } // End loop over fams (f).
  } // End loop over clusters (i).
  
  for(i=0; i<b->nclust; i++) { b->bi[i]->shuffle = false; }
  return shuffled;
}

/* b_p_update:
 Calculates the abundance p-value for each family in the clustering.
 Depends on the lambda between the fam and its cluster, and the reads of each.
*/
void b_p_update(B *b) {
  int i, f;
  Fam *fam;
  for(i=0;i<b->nclust;i++) {
    for(f=0;f<b->bi[i]->nfam;f++) {
      fam = b->bi[i]->fam[f];
      fam->p =  get_pA(fam, b->bi[i]);
      
      if(b->use_singletons) { // Calculate singleton pval (pS) from cluster lookup
        fam->pS = get_pS(fam, b->bi[i], b);
      }
      
    } // for(f=0;f<b->bi[i]->nfam;f++)
  } // for(i=0;i<b->nclust;i++)
}

/* b_bud:
 Finds the minimum p-value. If significant, creates a new cluster and moves the
 raws from the fam with the minimum p-value to the new cluster.
 Returns index of new cluster, or 0 if no new cluster added.
*/

int b_bud(B *b, double min_fold, int min_hamming, bool verbose) {
  int i, f, r;
  int mini, minf, totfams, minreads;
  double minp = 1.0, minlam=1.0;
  double pA=1.0, pS=1.0;
  double fold;
  int hamming;
  Fam *fam;
  Sub *birth_sub;

  // Find i, f indices and value of minimum pval.
  mini=-999; minf=-999; minreads=0; minp=1.0; totfams=0;
  for(i=0;i<b->nclust;i++) {
    for(f=0; f<b->bi[i]->nfam; f++) {
      totfams++;
      fam = b->bi[i]->fam[f];
      
      if(!fam->sub) {  // This shouldn't happen
        Rprintf("Warning: Fam has null sub in b_bud.\n");
        continue; 
      }

      // Calculate the fold over-abundnace and the hamming distance to this fam
      if(fam->sub->nsubs >= min_hamming) { // Only those passing the hamming/fold screens can be budded
        if(min_fold <= 1 || ((double) fam->reads) >= min_fold * b->bi[i]->e[fam->center->index]) {  
          if((fam->p < minp) ||
            ((fam->p == minp && fam->reads > minreads))) { // Most significant
            mini = i; minf = f;
            minp = fam->p;
            minreads = fam->reads;
          }
        }
      }
    }
  }
  
  // Bonferoni correct the abundance pval by the number of fams and compare to OmegaA
  // (quite conservative, although probably unimportant given the abundance model issues)
  pA = minp*totfams;
  if(pA < b->omegaA && mini >= 0 && minf >= 0) {  // A significant abundance pval
    birth_sub = sub_copy(b->bi[mini]->fam[minf]->sub); // do this before bi_pop_fam destructs the sub
    fam = bi_pop_fam(b->bi[mini], minf);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "A");
    b->bi[i]->birth_pval = pA;
    b->bi[i]->birth_fold = fam->reads/b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_e = b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_sub = birth_sub;
    
    // Move raws into new cluster, could be more elegant but this works.
    for(r=0;r<fam->nraw;r++) {
      bi_shove_raw(b->bi[i], fam->raw[r]);
    }

    if(verbose) { 
      double qave = 0.0;
      if(fam->center->qual) {
        for(int s=0;s<birth_sub->nsubs;s++) {
          for(int r=0;r<fam->nraw;r++) {
            qave += (fam->raw[r]->reads * fam->raw[r]->qual[birth_sub->pos[s]]);
          }
        }
        qave = qave/((double)birth_sub->nsubs * fam->reads);
      }
      Rprintf("\nNew cluster from Raw %i in C%iF%i: ", fam->center->index, mini, minf);
      fold = ((double) fam->reads)/b->bi[mini]->e[fam->center->index];
      Rprintf(" p*=%.3e, n/E(n)=%.1e (%.1e fold per sub)\n", pA, fold, fold/birth_sub->nsubs);
      Rprintf("Reads: %i, E: %.2e, Nsubs: %i, Ave Qsub:%.1f\n", fam->reads, b->bi[mini]->e[fam->center->index], birth_sub->nsubs, qave);
      Rprintf("%s\n", ntstr(birth_sub->key));
    }
    
    bi_center_update(b->bi[i]);
    fam_free(fam);
    return i;
  }

  // No significant abundance pval
  if(!b->use_singletons) {
    return 0;
  }
  
  // Using singletons, so find minimum pS
  mini=-999; minf=-999; minp=1.0; minlam = 1.0;
  for(i=0;i<b->nclust;i++) {
    for(f=0; f<b->bi[i]->nfam; f++) {
      // Calculate the fold over-abundance and the hamming distance to this fam
      fam = b->bi[i]->fam[f];
      fold = ((double) fam->reads)/b->bi[i]->e[fam->center->index];
      if(fam->sub) {
        hamming = fam->sub->nsubs;
      } else { 
        Rprintf("Warning: Fam has null sub in b_bud.\n");
        hamming = 1; // An unmotivated default here
      }
      
      if(fold >= min_fold && hamming >= min_hamming) {  // Hamming/fold screen
        if((fam->pS < minp) ||
          ((fam->pS == minp) && (fam->lambda < minlam))){ // Most significant
          mini = i; minf = f;
          minp = fam->pS;
          minlam = fam->lambda;
        } 
      }
    }
  }

  // Bonferoni correct the singleton pval by the number of reads and compare to OmegaS
  // DADA_ML MATCH: minp*b->nclust*b->bi[i]->reads < b->omegaS
  pS = minp * b->reads;
  if(pS < b->omegaS && mini >= 0 && minf >= 0) {  // A significant singleton pval
    birth_sub = sub_copy(b->bi[mini]->fam[minf]->sub); // do this before bi_pop_fam destructs the sub
    fam = bi_pop_fam(b->bi[mini], minf);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "S");
    b->bi[i]->birth_pval = pS;
    b->bi[i]->birth_fold = fam->reads/b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_e = b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_sub = birth_sub;
    
    // Move raws into new cluster, could be more elegant but this works.
    for(r=0;r<fam->nraw;r++) {
      bi_shove_raw(b->bi[i], fam->raw[r]);
    }

    if(verbose) { 
      Rprintf("\nNew cluster from Raw %i in C%iF%i: p*=%.3e (SINGLETON: lam=%.3e)\n", fam->center->index, mini, minf, pS, fam->lambda);
    }
    
    bi_center_update(b->bi[i]);
    fam_free(fam);
    return i;
  }

  // No significant abundance or singleton pval
  if(verbose) { Rprintf("\nNo significant pval, no new cluster.\n"); }
  return 0;
}

/* fam_center_update:
 Takes a fam object, and calculates and assigns its center sequence.
 Currently this is done trivially by choosing the most abundant sequence.
 NOTE: IT IS NOT CLEAR THAT THIS IS THE SAME AS THE MATLAB FUNCTIONALITY
 */
void fam_center_update(Fam *fam) {
  int max_reads = 0;
  int r, maxr= -99;
  
  for(r=0;r<fam->nraw;r++) {
    if(fam->raw[r]->reads > max_reads) {
      maxr=r;
      max_reads = fam->raw[r]->reads;
    }
  }
  
  if(maxr>=0) {       // Make sure there is at least one raw before strcpy
    fam->center = fam->raw[maxr];
    strcpy(fam->seq, fam->raw[maxr]->seq);
  } else {            // No raw, make empty string.
    fam->center = NULL;
    fam->seq[0] = '\0';
  }
}

/* Bi_center_update:
 Takes a Bi object, and calculates and assigns its center sequence.
 Currently this is done trivially by choosing the most abundant sequence.
 If center changes flags update_fam and update_lambda
*/
void bi_center_update(Bi *bi) {
  int max_reads = 0;
  int f, r;
  int maxf = 0, maxr= 0;
  
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      if(bi->fam[f]->raw[r]->reads > max_reads) { // Most abundant
         maxf=f; maxr=r;
         max_reads = bi->fam[f]->raw[r]->reads;
      }
    }
  }
  
  if(strcmp(bi->seq, bi->fam[maxf]->raw[maxr]->seq) != 0) {  // strings differ
    bi->update_lambda = true;
    bi->update_fam = true;
    if(VERBOSE) {
      Rprintf("bi_c_u: New center raw %i (%i)\n", bi->fam[maxf]->raw[maxr]->index, bi->fam[maxf]->raw[maxr]->reads);
      Rprintf("\tNew(%i): %s\n", (int) bi->fam[maxf]->raw[maxr]->length, ntstr(bi->fam[maxf]->raw[maxr]->seq));
      Rprintf("\tOld(%i): %s\n", (int) strlen(bi->seq), ntstr(bi->seq));
    }
    bi->center = bi->fam[maxf]->raw[maxr];
    strcpy(bi->seq,bi->fam[maxf]->raw[maxr]->seq);
  }
}

/* B_center_update:
 Updates all Bi consensi.
 */
void b_center_update(B *b) {
  for (int i=0; i<b->nclust; i++) {
    bi_center_update(b->bi[i]);
  }
}

/* Bi_make_consensus:
  Uses the alignments to the cluster center to construct a consensus sequence, which
  is then assigned to bi->seq
*/
void bi_make_consensus(Bi *bi, bool use_quals) {
  int i,f,r,s,nti0,nti1,pos;
  double counts[4][SEQLEN] = {{0.0}};
  Sub *sub = NULL;
  Raw *raw;
  
  if(!bi->center) { return; }
  int len = bi->center->length;
  if(len >= SEQLEN) { return; }
  
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      raw = bi->fam[f]->raw[r];
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
    } // for(r=0;r<bi->fam[f]->nraw;r++)
  } // for(f=0;f<bi->nfam;f++)
  
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
