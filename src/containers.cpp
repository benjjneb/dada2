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

/********* CONSTRUCTORS AND DESTRUCTORS *********/

// The constructor for the Raw object.
Raw *raw_new(char *seq, double *qual, unsigned int reads, bool prior) {
  // Allocate
  Raw *raw = (Raw *) malloc(sizeof(Raw)); //E
  if (raw == NULL)  Rcpp::stop("Memory allocation failed.");
  raw->seq = (char *) malloc(strlen(seq)+1); //E
  if (raw->seq == NULL)  Rcpp::stop("Memory allocation failed.");
  // Assign sequence and associated properties
  strcpy(raw->seq, seq);
  raw->length = strlen(seq);
  raw->reads = reads;
  raw->prior = prior;
  // Allocate and copy quals (quals downgraded to uint8_t here for memory savings)
  if(qual) { 
    raw->qual = (uint8_t *) malloc(raw->length * sizeof(uint8_t)); //E
    if (raw->qual == NULL)  Rcpp::stop("Memory allocation failed.");
    for(size_t i=0;i<raw->length;i++) { raw->qual[i] = (uint8_t) round(qual[i]); }
  } else {
    raw->qual = NULL;
  }
  raw->p = 0.0;
  raw->E_minmax = -999.0;
  raw->lock = false;
  raw->correct = true;
  return raw;
}

// The destructor for the Raw object.
void raw_free(Raw *raw) {
  free(raw->seq);
  if(raw->qual) { free(raw->qual); }
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
  bi->check_locks = true;
  
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
B *b_new(Raw **raws, unsigned int nraw, double omegaA, double omegaP, bool use_quals) {
  unsigned int index;
  
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
  b->omegaA = omegaA;
  b->omegaP = omegaP;
  b->use_quals = use_quals;
  
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
  b->bi[0]->birth_from = 0;
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
