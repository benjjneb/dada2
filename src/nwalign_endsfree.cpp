#include <string.h>
#include <stdlib.h>
#include "dada.h"
// [[Rcpp::interfaces(cpp)]]

/************* KMERS *****************
 * Current Kmer implementation assumes A/C/G/T only.
 */

double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  int i;
  int n_kmer = 2 << (2*k); // 4^k kmers
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
  size_t n_kmers = (2 << (2*k));  // 4^k kmers
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

char **raw_align(Raw *raw1, Raw *raw2, int score[4][4], int gap_p, bool use_kmers, double kdist_cutoff, int band) {
  char **al;
  double kdist;
  
  if(use_kmers) {
    kdist = kmer_dist(raw1->kmer, raw1->length, raw2->kmer, raw2->length, KMER_SIZE);
  }
  
  if(use_kmers && kdist > kdist_cutoff) {
    al = NULL;
  } else {
    al = nwalign_endsfree(raw1->seq, raw2->seq, score, gap_p, band);
  }

  return al;
}

/* note: input sequence must end with string termination character, '\0' */
char **nwalign_endsfree(char *s1, char *s2, int score[4][4], int gap_p, int band) {
  static size_t nnw = 0;
  int i, j;
  int l, r;
  int len1 = strlen(s1);
  int len2 = strlen(s2);
  int d[len1 + 1][len2 + 1]; // d: DP matrix
  int p[len1 + 1][len2 + 1]; // backpointer matrix with 1 for diagonal, 2 for left, 3 for up.
  int diag, left, up;
      
  // Fill out left columns of d, p.
  for (i = 0; i <= len1; i++) {
    d[i][0] = 0; // ends-free gap
    p[i][0] = 3;
  }
  
  // Fill out top rows of d, p.
  for (j = 0; j <= len2; j++) {
    d[0][j] = 0; // ends-free gap
    p[0][j] = 2;
  }
  
  for (j = 0; j <= len2; j++) {
    d[0][j] = 0; // ends-free gap
    p[0][j] = 2;
  }

  // Fill out band boundaries of d.
  if(band) {
    for(i=0;i<=len1;i++) {
      if(i-band-1 >= 0) { d[i][i-band-1] = -999; }
      if(i+band+1 <= len2) { d[i][i+band+1] = -999; }
    }
  }
  
  // Fill out the body of the DP matrix.
  for (i = 1; i <= len1; i++) {
    if(band) {
      l = i-band; if(l < 1) { l = 1; }
      r = i+band; if(r>len2) { r = len2; }
    } else { l=1; r=len2; }

    for (j = l; j <= r; j++) {
      // Score for the left move.
      if (i == len1)
        left = d[i][j-1]; // Ends-free gap.
      else
        left = d[i][j-1] + gap_p;
      
      // Score for the up move.
      if (j == len2) 
        up = d[i-1][j]; // Ends-free gap.
      else
        up = d[i-1][j] + gap_p;

      // Score for the diagonal move.
      diag = d[i-1][j-1] + score[s1[i-1]-1][s2[j-1]-1];
      
      // Break ties and fill in d,p.
      if (up >= diag && up >= left) {
        d[i][j] = up;
        p[i][j] = 3;
      } else if (left >= diag) {
        d[i][j] = left;
        p[i][j] = 2;
      } else {
        d[i][j] = diag;
        p[i][j] = 1;
      }
    }
  }
    
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) calloc( len1 + len2 + 1, sizeof(char)); //initialize to max possible length //E
  al[1] = (char *) calloc( len1 + len2 + 1, sizeof(char)); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");

  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;  

  while ( i > 0 || j > 0 ) {
    switch ( p[i][j] ) {
    case 1:
      al[0][len_al] = s1[--i];
      al[1][len_al] = s2[--j];
      break;
    case 2:
      al[0][len_al] = 6;
      al[1][len_al] = s2[--j];
      break;
    case 3:
      al[0][len_al] = s1[--i];
      al[1][len_al] = 6;
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  
  // Remove empty tails of strings.
  al[0] = (char *) realloc(al[0],len_al+1); //E
  al[1] = (char *) realloc(al[1],len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");
  
  // Reverse the alignment strings (since traced backwards).
  char temp;
  for (i = 0, j = len_al - 1; i <= j; i++, j--) {
    temp = al[0][i];
    al[0][i] = al[0][j];
    al[0][j] = temp;
    temp = al[1][i];
    al[1][i] = al[1][j];
    al[1][j] = temp;
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
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
  int i, i0, i1, bytes, align_length, len0, nsubs;
  bool is_nt0, is_nt1;
  char *al0, *al1; // dummy pointers to the sequences in the alignment
  
  // define dummy pointer to key string
  char *pkey;
  int key_size = 0; // key buffer bytes written
  
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
Sub *sub_new(Raw *raw0, Raw *raw1, int score[4][4], int gap_p, bool use_kmers, double kdist_cutoff, int band) {
  int s;
  char **al;
  Sub *sub;

  al = raw_align(raw0, raw1, score, gap_p, use_kmers, kdist_cutoff, band);
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
