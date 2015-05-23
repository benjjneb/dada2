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
void bi_consensus_update(Bi *bi);
void fam_consensus_update(Fam *fam);
void bi_fam_update(Bi *bi, double score[4][4], double gap_pen, int band_size, bool use_quals);

/*
 raw_new:
 The constructor for the Raw object.
 */
Raw *raw_new(char *seq, int reads) {
  Raw *raw = (Raw *) malloc(sizeof(Raw)); //E
  if (raw == NULL)  Rcpp::stop("Memory allocation failed!\n");
  raw->seq = (char *) malloc(strlen(seq)+1); //E
  if (raw->seq == NULL)  Rcpp::stop("Memory allocation failed!\n");
  strcpy(raw->seq, seq);
  raw->qual = NULL;
  raw->kmer = get_kmer(seq, KMER_SIZE);
  raw->reads = reads;
  return raw;
}

// raw_qual_new: A constructor for Raw objects with quals
Raw *raw_qual_new(char *seq, double *qual, int reads) {
  int seqlen = strlen(seq);
  Raw *raw = (Raw *) malloc(sizeof(Raw)); //E
  if (raw == NULL)  Rcpp::stop("Memory allocation failed!\n");
  raw->seq = (char *) malloc(seqlen+1); //E
  if (raw->seq == NULL)  Rcpp::stop("Memory allocation failed!\n");
  strcpy(raw->seq, seq);
  raw->qual = NULL;
  if(qual) {
    raw->qual = (float *) malloc(seqlen * sizeof(float)); //E
    if (raw->qual == NULL)  Rcpp::stop("Memory allocation failed!\n");
    for(int i=0;i<seqlen;i++) { raw->qual[i] = qual[i]; }
  } else {
    Rcpp::stop("Error: NULL qual provided to raw_qual_new constructor.\n");
  }
  raw->kmer = get_kmer(seq, KMER_SIZE);
  raw->reads = reads;
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
  if (fam == NULL)  Rcpp::stop("Memory allocation failed!\n");
  fam->seq = (char *) malloc(SEQLEN); //E
  if (fam->seq == NULL)  Rcpp::stop("Memory allocation failed!\n");
  fam->raw = (Raw **) malloc(RAWBUF * sizeof(Raw *)); //E
  if (fam->raw == NULL)  Rcpp::stop("Memory allocation failed!\n");
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
    if (fam->raw == NULL)  Rcpp::stop("Memory allocation failed!\n");
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
    printf("fam_pop_raw: Not enough raws %i (%i)\n", r, fam->nraw);
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
    }
  } else {  // Trying to pop an out of range raw
    printf("bi_pop_raw: not enough fams %i (%i)\n", f, bi->nfam);
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
    }
  } else {  // Trying to pop an out of range fam
    printf("bi_pop_fam: not enough fams %i (%i)\n", f, bi->nfam);
    pop = NULL;
  }
  return pop;
}

/*
 bi_shove_raw:
 Hackily just adds raw to fam[0].
 consensus_update and fam_update are needed on this cluster afterwards.
 */
void bi_shove_raw(Bi *bi, Raw *raw) {
  if(bi->nfam == 0) {
    bi_add_fam(bi, fam_new());
  }
  fam_add_raw(bi->fam[0], raw);
  bi->nraw++;
  bi->reads += raw->reads;
  bi->update_fam = true;
}

/*
 bi_add_raw:
 Adds raw to appropriate fam. Makes new fam if necessary.
 ISSUE WITH ALIGN, NEED ERR. UNUSED FOR NOW.
 UNUSED UNUSED UNUSED UNUSED UNUSED UNUSED
 */
void bi_add_raw(Bi *bi, Raw *raw) {
//  printf("bi_add_raw\n");
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
  bi->seq = (char *) malloc(SEQLEN); //E
  if (bi->seq == NULL)  Rcpp::stop("Memory allocation failed!\n");
  strcpy(bi->seq, "");
  bi->fam = (Fam **) malloc(FAMBUF * sizeof(Fam *)); //E
  if (bi->fam == NULL)  Rcpp::stop("Memory allocation failed!\n");
  bi->maxfam = FAMBUF;
  bi->sub = (Sub **) malloc(totraw * sizeof(Sub *)); //E
  if (bi->sub == NULL)  Rcpp::stop("Memory allocation failed!\n");
  for(int i=0;i<totraw;i++) { bi->sub[i] = NULL; }   // Init to null pointers
  bi->lambda = (double *) malloc(totraw * sizeof(double)); //E
  if (bi->lambda == NULL)  Rcpp::stop("Memory allocation failed!\n");
  bi->e = (double *) malloc(totraw * sizeof(double)); //E
  if (bi->e == NULL)  Rcpp::stop("Memory allocation failed!\n");
  bi->totraw = totraw;
  bi->sm = sm_new(MIN_BUCKETS); //E
  if (bi->sm == NULL)  Rcpp::stop("Memory allocation failed!\n");
  bi->update_lambda = true;
  bi->update_fam = true;
  bi->shuffle = true;
  
  bi->nfam = 0;
  bi->reads = 0;
  bi->nraw = 0;
  return bi;
}

/* bi_free:
 Destructs bi objects.
*/
void bi_free(Bi *bi) {
  for(int f=0;f<bi->nfam;f++) { fam_free(bi->fam[f]); }
  free(bi->fam);
  free(bi->seq);
  for(int i=0;i<bi->totraw;i++) { sub_free(bi->sub[i]); }
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
    if (bi->fam == NULL)  Rcpp::stop("Memory allocation failed!\n");
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
  bi->reads = reads;
  bi->nraw = nraw;
}

/* B_new:
 The constructor for the B object. Takes in a Uniques object.
 Places all sequences into the same family within one cluster.
*/
B *b_new(Raw **raws, int nraw, double score[4][4], double gap_pen, double omegaA, bool use_singletons, double omegaS, int band_size, bool use_quals) {
  int i, j, nti;
  size_t index;

  // Allocate memory
  B *b = (B *) malloc(sizeof(B)); //E
  if (b == NULL)  Rcpp::stop("Memory allocation failed!\n");
  b->bi = (Bi **) malloc(CLUSTBUF * sizeof(Bi *)); //E
  if (b->bi == NULL)  Rcpp::stop("Memory allocation failed!\n");
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
    for(i=0;i<strlen(b->raw[index]->seq);i++) {
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
 Initialized b with all raws in one cluster, with consensus, flagged to update lambda and fam.
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
  b->nalign;
  b->nshroud;

  // Add all raws to that cluster
  for (size_t index=0; index<b->nraw; index++) {
    bi_shove_raw(b->bi[0], b->raw[index]);
  }

  bi_census(b->bi[0]);
  bi_consensus_update(b->bi[0]); // Makes cluster consensus sequence
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
    if (b->bi == NULL)  Rcpp::stop("Memory allocation failed!\n");
    b->maxclust+=CLUSTBUF;
  }
  b->bi[b->nclust] = bi;
  bi->i = b->nclust;
  return(b->nclust++);
}

/*
 lambda_update:
 updates the alignments and lambda of all raws to Bi with
 updated consensus sequences. 
*/
void b_lambda_update(B *b, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat) {
  int i;
  size_t index;
  double lambda;
  Sub *sub;
  
  for (i = 0; i < b->nclust; i++) {
    if(b->bi[i]->update_lambda) {   // consensus sequence for Bi[i] has changed
      // update alignments and lambda of all raws to this sequence
      if(tVERBOSE) { printf("C%iLU:", i); }
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
        if(index == TARGET_RAW) { printf("lam(TARG)=%.2e; ", b->bi[i]->lambda[index]); }
      }
      b->bi[i]->update_lambda = false;
      if(!b->bi[i]->update_fam) {
        printf("Warning: Lambda updated but update_fam flag not set in C%i\n", i);
        b->bi[i]->update_fam = true;
      }
    } // if(b->bi[i]->update_lambda)
  }
  b_e_update(b);
}

/* bi_fam_update(Bi *bi, ...):
   Creates and allocates new fams based on the sub between each raw and the
   cluster consensus.
  Currently completely destructive of old fams.
   */
void bi_fam_update(Bi *bi, double score[4][4], double gap_pen, int band_size, bool use_quals) {
  int f, r, result, r_c;
  Sub *sub;
  char buf[10];
  
  sm_delete(bi->sm);
  bi->sm = sm_new(MIN_BUCKETS + (int) (BUCKET_SCALE * bi->nraw/2));  // n_buckets scales with # of raws
  // More buckets = more memory, but less collision (and therefore less costly strcmp w/in buckets).
  bi_census(bi);
  
  // Make list of pointers to the raws
  Raw **raws = (Raw **) malloc(bi->nraw * sizeof(Raw *)); //E
  if (raws == NULL)  Rcpp::stop("Memory allocation failed!\n");
  r_c=0;
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      raws[r_c] = bi->fam[f]->raw[r];
      r_c++;
    }
  }
  
  if(r_c != bi->nraw) {
    printf("Warning: bi_fam_update --- nraw inconsistent (%i, %i).\n", r_c, bi->nraw);
  }
  
  // Destruct old fams
  for(f=0;f<bi->nfam;f++) { fam_free(bi->fam[f]); }
  bi->nfam = 0;
  
  // Guarantee that non-NULL sub objects exist for all raws to the cluster center
  // NULL subs can occur (rarely) with kmers if cluster center changes (in shuffle) and 
  //    now exceeds kmerdist to another raw in the cluster
  for(r_c=0;r_c<bi->nraw;r_c++) {
    if(!(bi->sub[raws[r_c]->index])) { // Protect from and replace null subs
      if(tVERBOSE) { printf("Warning: bi_fam_update hit a null sub.\n"); }
      bi->sub[raws[r_c]->index] = sub_new(bi->center, raws[r_c], score, gap_pen, false, 1., band_size);
      // DOESN'T UPDATE LAMBDA HERE, which seems fine, as its most "fair" that it keeps its 0 lambda from being outside kmerdist
      ///! Check that lambda is in fact zero, warn otherwise.
      if(bi->lambda[raws[r_c]->index] != 0) {
        printf("Warning: Unexepected non-zero lambda %.2f for raw outside kmerdist in fam_update.\n", bi->lambda[raws[r_c]->index]);
      }
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
  
  // consensus_update/align the fams
  for(f=0;f<bi->nfam;f++) {
    fam_consensus_update(bi->fam[f]);
    sub = bi->sub[bi->fam[f]->center->index];

    if(!sub) { // Protect from null subs, but this should never arise...
      Rcpp::stop("Error: bi_fam_update hit a null sub. THIS SHOULDNT HAPPEN!!!!!\n");
    }
    
    bi->fam[f]->sub = sub;
    bi->fam[f]->lambda = bi->lambda[bi->fam[f]->center->index];
  }
  if(tVERBOSE) printf("(nraw=%d,nfam=%d), ", bi->nraw, bi->nfam);
  free(raws);
  bi->update_fam = false;
}

void b_fam_update(B *b) {
  for (int i=0; i<b->nclust; i++) {
    if(b->bi[i]->update_fam) {  // Consensus has changed or a raw has been shoved EITHER DIRECTION breaking the fam structure
      if(tVERBOSE) { printf("C%iFU:", i); }
      bi_fam_update(b->bi[i], b->score, b->gap_pen, b->band_size, b->use_quals);
    }
  }
}

void b_e_update(B *b) {
  double e;
  for(int i=0; i < b->nclust; i++) {
    for(int index=0; index < b->nraw; index++) {
      e = b->bi[i]->lambda[index]*b->bi[i]->reads;
      if(e != b->bi[i]->e[index]) { // E changed
        b->bi[i]->e[index] = e;
        b->bi[i]->shuffle = true;
      }
    }
  }
}

/* b_shuffle:
 move each sequence to the bi that produces the highest expected
 number of that sequence. The center of a Bi cannot leave.
*/
bool b_shuffle(B *b) {
  int i, f, r, j;
  int ibest, index;
  bool shuffled = false;
  Raw *raw;
  double e, maxe;
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
          for(j=0;j<b->nclust; j++) {
            if(b->bi[j]->shuffle) {
              e = b->bi[j]->e[index];
              if(e > maxe) {
                maxe = e;
                ibest = j;
              }
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
            if(tVERBOSE) {
              printf("Warning: Shuffle blocked the center of a Bi from leaving.\n");
              printf("Attempted: Raw %i from C%i to C%i (%.4e (lam=%.2e,n=%i) -> %.4e (%s: lam=%.2e,n=%i))\n", \
                  index, i, ibest, \
                  b->bi[i]->e[index], b->bi[i]->lambda[index], b->bi[i]->reads, \
                  b->bi[ibest]->e[index], b->bi[ibest]->sub[index]->key, b->bi[ibest]->lambda[index], b->bi[ibest]->reads);
            }
          } else { // Moving raw
            raw = bi_pop_raw(b->bi[i], f, r);
            bi_shove_raw(b->bi[ibest], raw);
            b->bi[i]->update_fam = true;
            b->bi[ibest]->update_fam = true;  // DUPLICATIVE FLAGGING FROM SHOVE_RAW FUNCTION
            shuffled = true;
//            if(VERBOSE) { printf("shuffle: Raw %i from C%i to C%i (%.4e -> %.4e)\n", index, i, ibest, b->bi[i]->e[index], b->bi[ibest]->e[index]); }
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

int b_bud(B *b, double min_fold, int min_hamming) {
  int i, f, r;
  int mini, minf, totfams, minreads;
  double minp = 1.0, minlam=1.0;
  double pA=1.0, pS=1.0;
  double fold;
  int hamming;
  Fam *fam;

  // Find i, f indices and value of minimum pval.
  mini=-999; minf=-999; minreads=0; minp=1.0; totfams=0;
  for(i=0;i<b->nclust;i++) {
    for(f=0; f<b->bi[i]->nfam; f++) {
      totfams++;
      fam = b->bi[i]->fam[f];
      
      // Calculate the fold over-abundnace and the hamming distance to this fam
      fold = ((double) fam->reads)/b->bi[i]->e[fam->center->index];
      if(fam->sub) {
        hamming = fam->sub->nsubs;
      } else { 
        printf("Warning: Fam has null sub in b_bud.\n");
        hamming = 1; // An unmotivated default here
      }
      
      if(fold >= min_fold && hamming >= min_hamming) {  // Only those passing the hamming/fold screens can be budded
        if((fam->p < minp) ||
          ((fam->p == minp && fam->reads > minreads))) { // Most significant
          mini = i; minf = f;
          minp = fam->p;
          minreads = fam->reads;
        } 
      }
    }
  }
  
  // Bonferoni correct the abundance pval by the number of fams and compare to OmegaA
  // (quite conservative, although probably unimportant given the abundance model issues)
  pA = minp*totfams;
  if(pA < b->omegaA && mini >= 0 && minf >= 0) {  // A significant abundance pval
    fam = bi_pop_fam(b->bi[mini], minf);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "A");
    b->bi[i]->birth_pval = pA;
    b->bi[i]->birth_fold = fam->reads/b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_e = b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_sub = fam->sub;
    
    // Move raws into new cluster, could be more elegant but this works.
    for(r=0;r<fam->nraw;r++) {
      bi_shove_raw(b->bi[i], fam->raw[r]);
    }

    if(tVERBOSE) { 
      double qave = 0.0;
      if(fam->center->qual) {
        for(int s=0;s<fam->sub->nsubs;s++) {
          for(int r=0;r<fam->nraw;r++) {
            qave += (fam->raw[r]->reads * fam->raw[r]->qual[fam->sub->pos[s]]);
          }
        }
        qave = qave/((double)fam->sub->nsubs * fam->reads);
      }
      printf("\nNew cluster from Raw %i in C%iF%i: ", fam->center->index, mini, minf);
      fold = ((double) fam->reads)/b->bi[mini]->e[fam->center->index];
      printf(" p*=%.3e, n/E(n)=%.1e (%.1e fold per sub)\n", pA, fold, fold/fam->sub->nsubs);
      printf("Reads: %i, E: %.2e, Nsubs: %i, Ave Qsub:%.1f\n", fam->reads, b->bi[mini]->e[fam->center->index], fam->sub->nsubs, qave);
      printf("%s\n", ntstr(fam->sub->key));
    }
    
    bi_consensus_update(b->bi[i]);
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
        printf("Warning: Fam has null sub in b_bud.\n");
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
    fam = bi_pop_fam(b->bi[mini], minf);
    i = b_add_bi(b, bi_new(b->nraw));
    strcpy(b->bi[i]->birth_type, "S");
    b->bi[i]->birth_pval = pS;
    b->bi[i]->birth_fold = fam->reads/b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_e = b->bi[mini]->e[fam->center->index];
    b->bi[i]->birth_sub = fam->sub;
    
    // Move raws into new cluster, could be more elegant but this works.
    for(r=0;r<fam->nraw;r++) {
      bi_shove_raw(b->bi[i], fam->raw[r]);
    }

    if(tVERBOSE) { 
      printf("\nNew cluster from Raw %i in C%iF%i: p*=%.3e (SINGLETON: lam=%.3e)\n", fam->center->index, mini, minf, pS, fam->lambda);
    }
    
    bi_consensus_update(b->bi[i]);
    fam_free(fam);
    return i;
  }

  // No significant abundance or singleton pval
  if(tVERBOSE) { printf("\nNo significant pval, no new cluster.\n"); }
  return 0;
}

/* fam_consensus_update:
 Takes a fam object, and calculates and assigns its consensus sequence.
 Currently this is done trivially by choosing the most abundant sequence.
 NOTE: IT IS NOT CLEAR THAT THIS IS THE SAME AS THE MATLAB FUNCTIONALITY
 */
void fam_consensus_update(Fam *fam) {
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

/* Bi_consensus_update:
 Takes a Bi object, and calculates and assigns its consensus sequence.
 Currently this is done trivially by choosing the most abundant sequence.
 If consensus changes...
    Updates .center and .self
    Flags update_fam and update_lambda
*/
void bi_consensus_update(Bi *bi) {
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
      printf("bi_c_u: New Consensus from raw %i (%i)\n", bi->fam[maxf]->raw[maxr]->index, bi->fam[maxf]->raw[maxr]->reads);
      printf("\tNew(%i): %s\n", (int) strlen(bi->fam[maxf]->raw[maxr]->seq), ntstr(bi->fam[maxf]->raw[maxr]->seq));
      printf("\tOld(%i): %s\n", (int) strlen(bi->seq), ntstr(bi->seq));
    }
    bi->center = bi->fam[maxf]->raw[maxr];
    strcpy(bi->seq,bi->fam[maxf]->raw[maxr]->seq);
  }
}

/* B_consensus_update:
 Updates all Bi consensi.
 */
void b_consensus_update(B *b) {
  for (int i=0; i<b->nclust; i++) {
    bi_consensus_update(b->bi[i]);
  }
}

