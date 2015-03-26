#include <string.h>
#include <stdlib.h>
#include "dada.h"
// [[Rcpp::interfaces(cpp)]]

double kmer_dist(int *kv1, int len1, int *kv2, int len2, int k) {
  double dot = 0;
  
  for(int i=0;i<(2 << (2*k)); i++) {
    dot += (kv1[i] < kv2[i] ? kv1[i] : kv2[i]);
  }
  dot = dot/((len1 < len2 ? len1 : len2) - k + 1.);
  
  if(dot < 0 || dot > 1) { printf("Bad dot: %.4e\n", dot); }

  return (1. - dot);
}

int *get_kmer(char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  int len = strlen(seq);
  size_t kmer = 0;
  size_t n_kmers = (2 << (2*k));  // 4^k kmers
  int *kvec = (int *) calloc(n_kmers,  sizeof(int));

  if(len <=0 || len > SEQLEN) {
    printf("Unexpected sequence length: %d\n", len);
  }

  for(i=0; i<len-k; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        printf("Unexpected nucleotide: %d\n", seq[j]);
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
      printf("Kmer index out of range!\n");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
  return kvec;
}

char **raw_align(Raw *raw1, Raw *raw2, double score[4][4], double gap_p, bool use_kmers, double kdist_cutoff, int band) {
  static long nalign=0;
  static long nshroud=0;
  char **al;
  double kdist;
  
  if(use_kmers) {
    kdist = kmer_dist(raw1->kmer, strlen(raw1->seq), raw2->kmer, strlen(raw2->seq), KMER_SIZE);
  } else { kdist = 0; }
  
  if(use_kmers && kdist > kdist_cutoff) {
    al = NULL;
    nshroud++;
  } else {
    al = nwalign_endsfree(raw1->seq, raw2->seq, score, gap_p, band);
  }
  nalign++;
  
  if(tVERBOSE && nalign%ALIGN_SQUAWK==0) { printf("%d alignments, %d shrouded\n", (int) nalign, (int) nshroud); }

  return al;
}

/* note: input sequence must end with string termination character, '\0' */
char **nwalign_endsfree(char *s1, char *s2, double score[4][4], double gap_p, int band) {
  static size_t nnw = 0;
  int i, j;
  int l, r;
  size_t len1 = strlen(s1);
  size_t len2 = strlen(s2);
  double d[len1 + 1][len2 + 1]; /* d: DP matrix */
  int p[len1 + 1][len2 + 1]; /* backpointer matrix with 1 for diagonal, 2 for left, 3 for up. */  
  double diag, left, up;
  
  int is_nt = 0;
  if(VERBOSE) {
    for(i=0;i<strlen(s1);i++) {
      if(s1[i] > 6) { is_nt=1; }
    }
    if(is_nt) { printf("NT STR IN ALIGN\n"); }
  }
  
  /* Fill out left columns of d, p. */
  for (i = 0; i <= len1; i++) {
    d[i][0] = 0; /* ends-free gap */
    p[i][0] = 3;
  }
  
  /* Fill out top rows of d, p. */
  for (j = 0; j <= len2; j++) {
    d[0][j] = 0; /* ends-free gap */
    p[0][j] = 2;
  }
  
  for (j = 0; j <= len2; j++) {
    d[0][j] = 0; /* ends-free gap */
    p[0][j] = 2;
  }

  /* Fill out band boundaries of d. */
  if(band) {
    for(i=0;i<=len1;i++) {
      if(i-band-1 >= 0) { d[i][i-band-1] = -999; }
      if(i+band+1 <= len2) { d[i][i+band+1] = -999; }
    }
  }
  
  /* Fill out the body of the DP matrix. */
  for (i = 1; i <= len1; i++) {
    if(band) {
      l = i-band; if(l < 1) { l = 1; }
      r = i+band; if(r>len2) { r = len2; }
    } else { l=1; r=len2; }

    for (j = l; j <= r; j++) {
      /* Score for the left move. */
      if (i == len1)
        left = d[i][j-1]; /* Ends-free gap. */
      else
        left = d[i][j-1] + gap_p;
      
      /* Score for the up move. */
      if (j == len2) 
        up = d[i-1][j]; /* Ends-free gap. */
      else
        up = d[i-1][j] + gap_p;

      /* Score for the diagonal move. */
      diag = d[i-1][j-1] + score[s1[i-1]-1][s2[j-1]-1];
//      printf("(%i, %i) --- diag:%.2f = %.2f + %.2f (%i, %i)\n", i, j, diag, d[i-1][j-1], score[s1[i-1]-1][s2[j-1]-1], s1[i-1]-1, s2[j-1]);
      
//      printf("(%i, %i) --- l:%.2f, u:%.2f, d:%.2f\n", i, j, left, up, diag);
      /* Break ties and fill in d,p. */
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
    
  /* Allocate memory to alignment strings. */
  char **al = (char **) malloc( 2 * sizeof(char *) );
  al[0] = (char *) calloc( len1 + len2 + 1, sizeof(char)); //initialize to max possible length
  al[1] = (char *) calloc( len1 + len2 + 1, sizeof(char));

  /* Trace back over p to form the alignment. */
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
      printf("NWALIGN: OOR p[%i][%i] = %i\n", i, j, p[i][j]);
    }
    len_al++;
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  
  /* Remove empty tails of strings. */
  al[0] = (char *) realloc(al[0],len_al+1);
  al[1] = (char *) realloc(al[1],len_al+1);
  
  // Reverse the alignment strings. // why?
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

/*
 al2subs:
 takes in an alignment represented as a char ** al. creates
 a Sub object from the substitutions of al[1] relative to al[0].
 that is, the identity of al[1] is stored at positions where it
 differs from al[0]
 */
Sub *al2subs(char **al) {
  int i, i0, i1, bytes;
  int align_length;
  int is_nt0, is_nt1;
  char *al0, *al1; /* dummy pointers to the sequences in the alignment */
  
  /* define dummy pointers to positions/nucleotides/keys */
//  int *ppos;
//  char *pnt0;
//  char *pnt1;
  char *pkey;
  int key_size = 0; /* key buffer bytes written*/
  
  if(!al) { // Null alignment (outside kmer thresh) -> Null sub
    Sub *sub = NULL;
    return sub;
  }
  
  // create Sub obect and initialize memory
  Sub *sub = (Sub *) malloc(sizeof(Sub));
  sub->map = (int *) malloc(strlen(al[0]) * sizeof(int));
  sub->pos = (int *) malloc(strlen(al[0]) * sizeof(int));
  sub->nt0 = (char *) malloc(strlen(al[0]));
  sub->nt1 = (char *) malloc(strlen(al[0]));
  sub->key = (char *) malloc(6*strlen(al[0])); // Will break if seqs get into the many thousands
  sub->nsubs=0;
  
  /* traverse the alignment and record substitutions while building the hash key*/
  pkey = sub->key;
  
  i0 = -1; i1 = -1;
  al0 = al[0]; al1 = al[1];
  align_length = strlen(al[0]);
  for(i=0;i<align_length;i++) {
    is_nt0 = ((al0[i] == 1) || (al0[i] == 2) || (al0[i] == 3) || (al0[i] == 4) || (al0[i] == 5)); // A/C/G/T/N (non-gap) in seq0
    is_nt1 = ((al1[i] == 1) || (al1[i] == 2) || (al1[i] == 3) || (al1[i] == 4) || (al1[i] == 5)); // A/C/G/T/N (non-gap) in seq1
    if(is_nt0) { i0++; }
    if(is_nt1) { i1++; }

    // Record to map
    if(is_nt0) {
      if(is_nt1) { 
        sub->map[i0] = i1;
      } else {
        sub->map[i0] = -999; // Indicates gap
      }
    }

    if(is_nt0 && is_nt1) { // Possible sub
      if((al0[i] != al1[i]) && (al0[i] != 5) && (al1[i] != 5)) { // Ns don't make subs
        sub->pos[sub->nsubs] = i0;
        sub->nt0[sub->nsubs] = al0[i];
        sub->nt1[sub->nsubs] = al1[i];
        
//        *ppos = i0;
//        *pnt0 = al0[i];
//        *pnt1 = al1[i];
//        ppos++; pnt0++; pnt1++;

        bytes = sprintf(pkey,"%c%d%c,",al0[i],i0,al1[i]);
        key_size += bytes;
        pkey += bytes;
        sub->nsubs++;
      }
    }
  }
  *pkey = '\0';
  int2nt(sub->key, sub->key);  // store sub keys as readable strings
    
  /* give back unused memory */
  sub->pos = (int *) realloc(sub->pos, sub->nsubs * sizeof(int));
  sub->nt0 = (char *) realloc(sub->nt0, sub->nsubs);
  sub->nt1 = (char *) realloc(sub->nt1, sub->nsubs);
  sub->key = (char *) realloc(sub->key, strlen(sub->key) + 1); /* +1 necessary for the "\0" */
//  printf("Key %i\t%i\t%s\n", key_size, strlen(sub->key), ntstr(sub->key));
  
  if(VERBOSE && sub->nsubs>150) {
    printf("AL2SUBS: Many subs (%i) \n", sub->nsubs);
  }
  return sub;
}

// Wrapper for al2subs(raw_align(...)) that manages memory and quality scores
Sub *sub_new(Raw *raw0, Raw *raw1, double score[4][4], double gap_p, bool use_kmers, double kdist_cutoff, int band, bool use_quals) {
  int s;
  char **al;
  Sub *sub;

  al = raw_align(raw0, raw1, score, gap_p, use_kmers, kdist_cutoff, band);
  sub = al2subs(al);

  if(sub) {
    sub->q0 = NULL;
    sub->q1 = NULL;
    if(use_quals) {
      sub->q0 = (double *) malloc(sub->nsubs * sizeof(double));
      sub->q1 = (double *) malloc(sub->nsubs * sizeof(double));
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
