#include <string.h>
#include <stdlib.h>
#include "dada.h"
// [[Rcpp::interfaces(cpp)]]

/************* KMERS *****************
 * Current Kmer implementation assumes A/C/G/T only.
 */

double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  int i;
  int n_kmer = 1 << (2*k); // 4^k kmers
  uint16_t dotsum = 0;
  double dot = 0.0;
  
  for(i=0;i<n_kmer; i++) {
    dotsum += (kv1[i] < kv2[i] ? kv1[i] : kv2[i]);
  }

  dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}

uint16_t *get_kmer(char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  int len = strlen(seq);
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  uint16_t *kvec = (uint16_t *) malloc(n_kmers * sizeof(uint16_t)); //E
  if (kvec == NULL)  Rcpp::stop("Memory allocation failed.");
  for(kmer=0;kmer<n_kmers;kmer++) { kvec[kmer] = 0; }

  if(len <=0 || len > SEQLEN) {
    Rcpp::stop("Unexpected sequence length.");
  }

  for(i=0; i<len-k; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        Rcpp::stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }
    
    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; } 
    else if(kmer >= n_kmers) {
      Rcpp::stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
  return kvec;
}

/************* ALIGNMENT *****************
 * Banded Needleman Wunsch
 */

char **raw_align(Raw *raw1, Raw *raw2, int score[4][4], int gap_p, int homo_gap_p, bool use_kmers, double kdist_cutoff, int band, bool vectorized_alignment) {
  char **al;
  double kdist;
  
  if(use_kmers) {
    kdist = kmer_dist(raw1->kmer, raw1->length, raw2->kmer, raw2->length, KMER_SIZE);
  }
  
  if(use_kmers && kdist > kdist_cutoff) {
    al = NULL;
  } else if(vectorized_alignment) { // ASSUMES SCORE MATRIX REDUCES TO MATCH/MISMATCH
//    al = nwalign_endsfree_vectorized(raw1->seq, raw2->seq, (int16_t) score[0][0], (int16_t) score[0][1], (int16_t) gap_p, band);
    al = nwalign_vectorized2(raw1->seq, raw2->seq, (int16_t) score[0][0], (int16_t) score[0][1], (int16_t) gap_p, 0, band);
  } else if(homo_gap_p != gap_p && homo_gap_p <= 0) {
    al = nwalign_endsfree_homo(raw1->seq, raw2->seq, score, gap_p, homo_gap_p, band);
  } else {
    al = nwalign_endsfree(raw1->seq, raw2->seq, score, gap_p, band);
  }

  return al;
}

/* note: input sequence must end with string termination character, '\0' */
char **nwalign_endsfree(const char *s1, const char *s2, int score[4][4], int gap_p, int band) {
  static size_t nnw = 0;
  int i, j;
  int l, r;
  unsigned int len1 = strlen(s1);
  unsigned int len2 = strlen(s2);
  int diag, left, up;
  
  unsigned int nrow = len1+1;
  unsigned int ncol = len2+1;
  int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
  int *p = (int *) malloc(nrow * ncol * sizeof(int)); //E
  if(d == NULL || p == NULL) Rcpp::stop("Memory allocation failed.");
  
  // Fill out left columns of d, p.
  for (i = 0; i <= len1; i++) {
    d[i*ncol] = 0; // ends-free gap
    p[i*ncol] = 3;
  }
  
  // Fill out top rows of d, p.
  for (j = 0; j <= len2; j++) {
    d[j] = 0; // ends-free gap
    p[j] = 2;
  }
  
  // Calculate left/right-bands in case of different lengths
  int lband, rband;
  if(len2 > len1) {
    lband = band;
    rband = band+len2-len1;
  } else if(len1 > len2) {
    lband = band+len1-len2;
    rband = band;
  } else {
    lband = band;
    rband = band;
  }
  
  // Fill out band boundaries of d.
  if(band>=0 && (band<len1 || band<len2)) {
    for(i=0;i<=len1;i++) {
      if(i-lband-1 >= 0) { d[i*ncol + i-lband-1] = -9999; }
      if(i+rband+1 <= len2) { d[i*ncol + i+rband+1] = -9999; }
    }
  }
  
  // Fill out the body of the DP matrix.
  for (i = 1; i <= len1; i++) {
    if(band>=0) {
      l = i-lband; if(l < 1) { l = 1; }
      r = i+rband; if(r>len2) { r = len2; }
    } else { l=1; r=len2; }

    for (j = l; j <= r; j++) {
      // Score for the left move.
      if (i == len1) {
        left = d[i*ncol + j-1]; // Ends-free gap.
      } else {
        left = d[i*ncol + j-1] + gap_p;
      }
      
      // Score for the up move.
      if (j == len2) {
        up = d[(i-1)*ncol + j]; // Ends-free gap.
      } else {
        up = d[(i-1)*ncol + j] + gap_p;
      }

      // Score for the diagonal move.
      diag = d[(i-1)*ncol + j-1] + score[s1[i-1]-1][s2[j-1]-1];
      
      // Break ties and fill in d,p.
      if (up >= diag && up >= left) {
        d[i*ncol + j] = up;
        p[i*ncol + j] = 3;
      } else if (left >= diag) {
        d[i*ncol + j] = left;
        p[i*ncol + j] = 2;
      } else {
        d[i*ncol + j] = diag;
        p[i*ncol + j] = 1;
      }
    }
  }
    
  char *al0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *al1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");

  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;  

  while ( i > 0 || j > 0 ) {
//    Rprintf("(%i, %i): p=%i, d=%i\n", i, j, p[i*ncol + j], d[i*ncol + j]);
    switch ( p[i*ncol + j] ) {
    case 1:
      al0[len_al] = s1[--i];
      al1[len_al] = s2[--j];
      break;
    case 2:
      al0[len_al] = 6;
      al1[len_al] = s2[--j];
      break;
    case 3:
      al0[len_al] = s1[--i];
      al1[len_al] = 6;
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  
  
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");

  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  // Free allocated memory
  free(d);
  free(p);
  free(al0);
  free(al1);
  
  nnw++;
  return al;
}

/* note: input sequence must end with string termination character, '\0' */
/* 08-17-15: MJR homopolymer free gapping version of ends-free alignment */
char **nwalign_endsfree_homo(const char *s1, const char *s2, int score[4][4], int gap_p, int homo_gap_p, int band) {
  static size_t nnw = 0;
  int i, j, k;
  int l, r;
  unsigned int len1 = strlen(s1);
  unsigned int len2 = strlen(s2);
  int diag, left, up;
  
  //find locations where s1 has homopolymer and put 1s in homo1
  unsigned char *homo1 = (unsigned char *) malloc(len1*sizeof(unsigned char)); //E
  unsigned char *homo2 = (unsigned char *) malloc(len2*sizeof(unsigned char)); //E
  if(homo1 == NULL || homo2 == NULL) Rcpp::stop("Memory allocation failed.");
  for (i=0,j=0;j<len1;j++) {
    if (j==len1-1 || s1[j]!=s1[j+1]) {
      for(k=i;k<=j;k++) {
        if (j-i>=2) {//min homopolymer length = 3
          homo1[k] = 1;
        } else {
          homo1[k] = 0;
        }
      }
      i = j+1;
    }
  }
  
  //find locations where s2 has homopolymer and put 1s in homo2
  for (i=0,j=0;j<len2;j++) {
    if (j==len2-1 || s2[j]!=s2[j+1]) {
      for(k=i;k<=j;k++) {
        if (j-i>=2) { //min homopolymer length = 3
          homo2[k] = 1;
        } else {
          homo2[k] = 0;
        }
      }
      i = j+1;
    }
  }

  unsigned int nrow = len1+1;
  unsigned int ncol = len2+1;
  int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
  int *p = (int *) malloc(nrow * ncol * sizeof(int)); //E
  if(d == NULL || p == NULL) Rcpp::stop("Memory allocation failed.");
  
  // Fill out left columns of d, p.
  for (i = 0; i <= len1; i++) {
    d[i*ncol] = 0; // ends-free gap
    p[i*ncol] = 3;
  }
  
  // Fill out top rows of d, p.
  for (j = 0; j <= len2; j++) {
    d[j] = 0; // ends-free gap
    p[j] = 2;
  }
  
  // Calculate left/right-bands in case of different lengths
  int lband, rband;
  if(len2 > len1) {
    lband = band;
    rband = band+len2-len1;
  } else if(len1 > len2) {
    lband = band+len1-len2;
    rband = band;
  } else {
    lband = band;
    rband = band;
  }
  
  // Fill out band boundaries of d.
  if(band>=0 && (band<len1 || band<len2)) {
    for(i=0;i<=len1;i++) {
      if(i-lband-1 >= 0) { d[i*ncol + i-lband-1] = -9999; }
      if(i+rband+1 <= len2) { d[i*ncol + i+rband+1] = -9999; }
    }
  }
  
  // Fill out the body of the DP matrix.
  for (i = 1; i <= len1; i++) {
    if(band>=0) {
      l = i-lband; if(l < 1) { l = 1; }
      r = i+rband; if(r>len2) { r = len2; }
    } else { l=1; r=len2; }
    
    for (j = l; j <= r; j++) {
      // Score for the left move.
      if (i == len1) {
        left = d[i*ncol + j-1]; // Ends-free gap.
      } else if (homo2[j-1]) {
        left = d[i*ncol + j-1] + homo_gap_p; //Homopolymer gap
      } else {
        left = d[i*ncol + j-1] + gap_p;
      }
      
      // Score for the up move.
      if (j == len2) {
        up = d[(i-1)*ncol + j]; // Ends-free gap.
      } else if (homo1[i-1]) {
          up = d[(i-1)*ncol + j] + homo_gap_p; //Homopolymer gap
      } else {
        up = d[(i-1)*ncol + j] + gap_p;
      }
      
      // Score for the diagonal move.
      diag = d[(i-1)*ncol + j-1] + score[s1[i-1]-1][s2[j-1]-1];
      
      // Break ties and fill in d,p.
      if (up >= diag && up >= left) {
        d[i*ncol + j] = up;
        p[i*ncol + j] = 3;
      } else if (left >= diag) {
        d[i*ncol + j] = left;
        p[i*ncol + j] = 2;
      } else {
        d[i*ncol + j] = diag;
        p[i*ncol + j] = 1;
      }
    }
  }
  
  char *al0 = (char *) malloc((len1+len2+1) * sizeof(char)); //E
  char *al1 = (char *) malloc((len1+len2+1) * sizeof(char)); //E
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;  
  
  while ( i > 0 || j > 0 ) {
    switch ( p[i*ncol + j] ) {
    case 1:
      al0[len_al] = s1[--i];
      al1[len_al] = s2[--j];
      break;
    case 2:
      al0[len_al] = 6;
      al1[len_al] = s2[--j];
      break;
    case 3:
      al0[len_al] = s1[--i];
      al1[len_al] = 6;
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  
  
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");
  
  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  // Free allocated memory
  free(d);
  free(p);
  free(homo1);
  free(homo2);
  free(al0);
  free(al1);
  
  nnw++;
  return al;
}


// Provided for the R nwalign function if endsfree=FALSE
// Not used within the dada method
// Separate function to avoid if statement within performance critical nwalign_endsfree
/* note: input sequence must end with string termination character, '\0' */
char **nwalign(const char *s1, const char *s2, int score[4][4], int gap_p, int band) {
  static size_t nnw = 0;
  int i, j;
  int l, r;
  unsigned int len1 = strlen(s1);
  unsigned int len2 = strlen(s2);
  int diag, left, up;
  
  unsigned int nrow = len1+1;
  unsigned int ncol = len2+1;
  int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
  int *p = (int *) malloc(nrow * ncol * sizeof(int)); //E
  if(d == NULL || p == NULL) Rcpp::stop("Memory allocation failed.");
  
  d[0] = 0;
  p[0] = 0; // Should never be queried
  
  // Fill out left columns of d, p.
  for (i = 1; i <= len1; i++) {
    d[i*ncol] = d[(i-1)*ncol] + gap_p;
    p[i*ncol] = 3;
  }
  
  // Fill out top rows of d, p.
  for (j = 1; j <= len2; j++) {
    d[j] = d[j-1] + gap_p;
    p[j] = 2;
  }
  
  // Calculate left/right-bands in case of different lengths
  int lband, rband;
  if(len2 > len1) {
    lband = band;
    rband = band+len2-len1;
  } else if(len1 > len2) {
    lband = band+len1-len2;
    rband = band;
  } else {
    lband = band;
    rband = band;
  }
  
  // Fill out band boundaries of d.
  if(band>=0 && (band<len1 || band<len2)) {
    for(i=0;i<=len1;i++) {
      if(i-lband-1 >= 0) { d[i*ncol + i-lband-1] = -9999; }
      if(i+rband+1 <= len2) { d[i*ncol + i+rband+1] = -9999; }
    }
  }
  
  // Fill out the body of the DP matrix.
  for (i = 1; i <= len1; i++) {
    if(band>=0) {
      l = i-lband; if(l < 1) { l = 1; }
      r = i+rband; if(r>len2) { r = len2; }
    } else { l=1; r=len2; }

    for (j = l; j <= r; j++) {
      // Score for the left move.
      left = d[i*ncol + j-1] + gap_p;
      
      // Score for the up move.
      up = d[(i-1)*ncol + j] + gap_p;

      // Score for the diagonal move.
      diag = d[(i-1)*ncol + j-1] + score[s1[i-1]-1][s2[j-1]-1];
      
      // Break ties and fill in d,p.
      if (up >= diag && up >= left) {
        d[i*ncol + j] = up;
        p[i*ncol + j] = 3;
      } else if (left >= diag) {
        d[i*ncol + j] = left;
        p[i*ncol + j] = 2;
      } else {
        d[i*ncol + j] = diag;
        p[i*ncol + j] = 1;
      }
    }
  }
    
  char *al0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *al1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");

  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;  

  while ( i > 0 || j > 0 ) {
    switch ( p[i*ncol + j] ) {
    case 1:
      al0[len_al] = s1[--i];
      al1[len_al] = s2[--j];
      break;
    case 2:
      al0[len_al] = 6;
      al1[len_al] = s2[--j];
      break;
    case 3:
      al0[len_al] = s1[--i];
      al1[len_al] = 6;
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  
  
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");

  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  // Free allocated memory
  free(d);
  free(p);
  free(al0);
  free(al1);
  
  nnw++;
  return al;
}

/************* SUBS *****************
 * Compressed storage for alignment.
 * Keeps only substitutions.
 * (could also consider CIGAR format)
 */

/*
 al2subs:
 takes in an alignment represented as a char ** al. creates
 a Sub object from the substitutions of al[1] relative to al[0].
 that is, the identity of al[1] is stored at positions where it
 differs from al[0]
 */
Sub *al2subs(char **al) {
  int i, i0, i1, align_length, len0, nsubs;
  bool is_nt0, is_nt1;
  char *al0, *al1; // dummy pointers to the sequences in the alignment
  
  // define dummy pointer to key string
  char *pkey;
  
  if(!al) { // Null alignment (outside kmer thresh) -> Null sub
    Sub *sub = NULL;
    return sub;
  }
  
  // create Sub obect and initialize memory
  Sub *sub = (Sub *) malloc(sizeof(Sub)); //E
  if (sub == NULL)  Rcpp::stop("Memory allocation failed.");

  // traverse alignment and find length of sq0 and nsubs for memory allocation
  len0 = 0; nsubs = 0;
  al0 = al[0]; al1 = al[1];
  align_length = strlen(al[0]);
  for(i=0;i<align_length;i++) {
    is_nt0 = ((al0[i] == 1) || (al0[i] == 2) || (al0[i] == 3) || (al0[i] == 4) || (al0[i] == 5)); // A/C/G/T/N (non-gap) in seq0
    is_nt1 = ((al1[i] == 1) || (al1[i] == 2) || (al1[i] == 3) || (al1[i] == 4) || (al1[i] == 5)); // A/C/G/T/N (non-gap) in seq1
    if(is_nt0) { len0++; }
    
    if(is_nt0 && is_nt1) { // Possible sub
      if((al0[i] != al1[i]) && (al0[i] != 5) && (al1[i] != 5)) { // Ns don't make subs
        nsubs++;
      }
    }
  }

  sub->len0 = len0;
  sub->map = (uint16_t *) malloc(len0 * sizeof(uint16_t)); //E
  sub->pos = (uint16_t *) malloc(nsubs * sizeof(uint16_t)); //E
  sub->nt0 = (char *) malloc(nsubs); //E
  sub->nt1 = (char *) malloc(nsubs); //E
  sub->key = (char *) malloc((6*nsubs) + 1); // Must be modified if SEQLEN > 1000 //E
  if (sub->map == NULL || sub->pos == NULL || sub->nt0 == NULL || sub->nt1 == NULL || sub->key == NULL) {
    Rcpp::stop("Memory allocation failed.");
  }
  sub->nsubs=0;
    
  // traverse the alignment and record substitutions while building the hash key
  pkey = sub->key;
  
  i0 = -1; i1 = -1;
  al0 = al[0]; al1 = al[1];
  for(i=0;i<align_length;i++) {
    is_nt0 = ((al0[i] == 1) || (al0[i] == 2) || (al0[i] == 3) || (al0[i] == 4) || (al0[i] == 5)); // A/C/G/T/N (non-gap) in seq0
    is_nt1 = ((al1[i] == 1) || (al1[i] == 2) || (al1[i] == 3) || (al1[i] == 4) || (al1[i] == 5)); // A/C/G/T/N (non-gap) in seq1
    if(is_nt0) { i0++; }
    if(is_nt1) { i1++; }

    // Record to map
    if(is_nt0) {
      if(is_nt1) { 
        sub->map[i0] = i1; // Assigning a signed into to a uint16_t...
      } else {
        sub->map[i0] = GAP_GLYPH; // Indicates gap
      }
    }

    if(is_nt0 && is_nt1) { // Possible sub
      if((al0[i] != al1[i]) && (al0[i] != 5) && (al1[i] != 5)) { // Ns don't make subs
        sub->pos[sub->nsubs] = i0;
        sub->nt0[sub->nsubs] = al0[i];
        sub->nt1[sub->nsubs] = al1[i];
        
        // Assuming space is available
        // Assuming sequences are <1000 nts long (3 digit positions)
        *pkey++ = al0[i];
        *pkey++ = '0' + i0/100;
        *pkey++ = '0' + (i0 % 100)/10;
        *pkey++ = '0' + (i0 % 10);
        *pkey++ = al1[i];
        *pkey++ = ',';
        sub->nsubs++;
      }
    }
  } // for(i=0;i<align_length;i++)
  *pkey = '\0';

  return sub;
}

// Wrapper for al2subs(raw_align(...)) that manages memory and qualities
Sub *sub_new(Raw *raw0, Raw *raw1, int score[4][4], int gap_p, int homo_gap_p, bool use_kmers, double kdist_cutoff, int band, bool vectorized_alignment) {
  int s;
  char **al;
  Sub *sub;

  al = raw_align(raw0, raw1, score, gap_p, homo_gap_p, use_kmers, kdist_cutoff, band, vectorized_alignment);
  sub = al2subs(al);

  if(sub) {
    sub->q0 = NULL;
    sub->q1 = NULL;
    if(raw0->qual && raw1->qual) {
      sub->q0 = (double *) malloc(sub->nsubs * sizeof(double)); //E
      sub->q1 = (double *) malloc(sub->nsubs * sizeof(double)); //E
      if (sub->q0 == NULL || sub->q1 == NULL) { Rcpp::stop("Memory allocation failed."); }
      
      for(s=0;s<sub->nsubs;s++) {
        sub->q0[s] = raw0->qual[sub->pos[s]];
        sub->q1[s] = raw1->qual[sub->map[sub->pos[s]]];
      }
    }
  }

  if(al) { // not a NULL align
    free(al[0]);
    free(al[1]);
    free(al);
  }

  return sub;
}

// Copies the given sub into a newly allocated sub object
Sub *sub_copy(Sub *sub) {
  int nsubs, len0;
  
  if(sub == NULL) { return(NULL); }
  nsubs = sub->nsubs;
  len0 = sub->len0;
  
  Sub *rsub = (Sub *) malloc(sizeof(Sub)); //E
  if (rsub == NULL)  Rcpp::stop("Memory allocation failed.");
  rsub->map = (uint16_t *) malloc(len0 * sizeof(uint16_t)); //E
  rsub->pos = (uint16_t *) malloc(nsubs * sizeof(uint16_t)); //E
  rsub->nt0 = (char *) malloc(nsubs); //E
  rsub->nt1 = (char *) malloc(nsubs); //E
  rsub->key = (char *) malloc((6*nsubs) + 1); // Must be modified if SEQLEN > 1000 //E
  if (rsub->map == NULL || rsub->pos == NULL || rsub->nt0 == NULL || rsub->nt1 == NULL || rsub->key == NULL) {
    Rcpp::stop("Memory allocation failed.");
  }
  
  rsub->nsubs = sub->nsubs;
  rsub->len0 = sub->len0;
  memcpy(rsub->map, sub->map, len0 * sizeof(uint16_t));
  memcpy(rsub->pos, sub->pos, nsubs * sizeof(uint16_t));
  memcpy(rsub->nt0, sub->nt0, nsubs);
  memcpy(rsub->nt1, sub->nt1, nsubs);
  memcpy(rsub->key, sub->key, (6*nsubs) + 1);

  if(sub->q0 && sub->q1) {
    rsub->q0 = (double *) malloc(nsubs * sizeof(double)); //E
    rsub->q1 = (double *) malloc(nsubs * sizeof(double)); //E
    if (rsub->q0 == NULL || rsub->q1 == NULL) { Rcpp::stop("Memory allocation failed."); }
    memcpy(rsub->q0, sub->q0, nsubs * sizeof(double));
    memcpy(rsub->q1, sub->q1, nsubs * sizeof(double));
  } else {
    rsub->q0 = NULL;
    rsub->q1 = NULL;
  }

  return rsub;
}

// Destructor for sub object
void sub_free(Sub *sub) {
  if(sub) { // not a NULL sub
    free(sub->key);
    free(sub->nt1);
    free(sub->nt0);
    free(sub->pos);
    free(sub->map);
    if(sub->q0) { free(sub->q0); }
    if(sub->q1) { free(sub->q1); }
    free(sub);
  }
}
