#include <Rcpp.h>
#include "dada.h"
using namespace Rcpp;

char **nwalign_endsfree_vectorized(char *s1, char *s2, int16_t score[4][4], int16_t gap_p, size_t band) {
  size_t row, col;
  size_t i,j;
  size_t index;
  size_t len1 = strlen(s1);
  size_t len2 = strlen(s2);
  int16_t d[800][800]; // d: DP matrix
  int16_t p[800][800]; // backpointer matrix with 1 for diagonal, 2 for left, 3 for up.
  // Too big for the Effing stack at 2*SEQLEN+1 x 2*SEQLEN+1
  int16_t dbuf[2*SEQLEN+1];
  int16_t pbuf[2*SEQLEN+1];
  int16_t diag_buf[SEQLEN];
  int16_t diag, left, up, entry, pentry;
  size_t center;
  size_t col_min, col_max;
  const int16_t *ptr_left, *ptr_diag, *ptr_up;
  
  // Simplifying to match/mismatch alignment
  int16_t match = score[0][0];
  int16_t mismatch = score[0][2];

  // assuming len1=len2 for now
  // assuming band is even for now
  if(band>=0 && band<len1) {
    center=1 + band/2;
  } else {
    center=1 + len1/2;
  }
  
  // Initialize top to -9999
  for(row=0;row<band;row++) {
    for(col=0;col<=center+band+1;col++) {
      d[row][col] = -9999;
      p[row][col] = 0; // Should never be queried
    }
  }
  
  // Fill out starting point
  d[0][center] = 0;
  p[row][col] = 0; // Should never be queried
  
  // Fill out "left" "column" of d, p.
  for(row=1, col=center-1; col>0; row+=2, col--) {
    d[row][col] = 0; // ends-free gap
    p[row][col] = 3;
    d[row+1][col] = 0; // ends-free gap
    p[row+1][col] = 3;
  }
  
  // Fill out "top" "row" of d, p.
  for (row=1,col=center+1; col<=center+band/2; row+=2, col++) {
    d[row][col-1] = 0; // ends-free gap
    p[row][col-1] = 2;
    d[row+1][col] = 0; // ends-free gap
    p[row+1][col] = 2;
  }
    
  // Fill out band boundaries
  if(band>=0 && band<len1) {
    for(row=0;row<=len1+len2;row++) {
      d[row][0] = -9999;
    }
  }
  if(band>=0 && band<len2) {
    for(row=0;row<=len1+len2;row+=2) {
      d[row][center+band/2+1] = -9999;
      d[row+1][center+band/2] = -9999;
      d[row+1][center+band/2+1] = -9999;
    }
  }
  
  // Fill out top wedge (Row 1 taken care of by ends-free)
  for(row=2;row<=band;row+=2) { // If band odd, could broach band in last row?
    // Fill out even row
    col_min = center-row/2+1;  // avoid ends-free part
    col_max = center+row/2-1;  // avoid ends-free part
    for(col=col_min,i=row-2,j=0;col<1+col_max;col++,i--,j++) {
      left = d[row-1][col-1] + gap_p;
      diag = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
      up = d[row-1][col] + gap_p;
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      d[row][col] = entry; // Vectorizes if this points to another array
      p[row][col] = pentry;
    }
    
    // Fill out odd row (row+1)
    col_min--; // Because avoided ends-free previously
    
    for(col=col_min,i=row-1,j=0;col<1+col_max;col++,i--,j++) { // avoid ends-free part
      left = d[row][col] + gap_p;
      diag = d[row-1][col] + (s1[i] == s2[j] ? match : mismatch);
      up = d[row][col+1] + gap_p;
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      d[row+1][col] = entry; // Vectorizes if this points to another array
      p[row+1][col] = pentry;
    }
  }

  // Fill out banded body
  for(row=band;row<=len1+len2-band;row+=2) { // Using fact that len1+len2=even, overlapping by 1 row with next part
    // Fill out even row
    col_min = 1;
    col_max = center+band/2;
    ptr_left = &d[row-1][col_min-1];
    ptr_diag = &d[row-2][col_min];
    ptr_up = &d[row-1][col_min];
//    #pragma clang loop vectorize(enable)
    for(col=col_min,i=row/2+band/2-1,j=row/2-band/2-1;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
    }
    for(col=col_min;col<1+col_max;col++) {
      left = (*ptr_left++) + gap_p;
      diag = diag_buf[col];
      up = (*ptr_up++) + gap_p;
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      dbuf[col] = entry; // Vectorizes if this points to another array
      pbuf[col] = pentry;
    }
    memcpy(&d[row][1], &dbuf[1], (band+1)*sizeof(int16_t));
    memcpy(&p[row][1], &pbuf[1], (band+1)*sizeof(int16_t));
    
    // Fill out odd row
    col_max--;
    ptr_left = &d[row][col_min];
    ptr_diag = &d[row-1][col_min];
    ptr_up = &d[row][col_min+1];

    for(col=col_min,i=row/2+band/2-1,j=row/2-band/2;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[row-1][col] + (s1[i] == s2[j] ? match : mismatch);
    }
    for(col=col_min;col<1+col_max;col++) {
      left = (*ptr_left++) + gap_p;
      diag = diag_buf[col];
      up = (*ptr_up++) + gap_p;
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      dbuf[col] = entry; // Vectorizes if this points to another array
      pbuf[col] = pentry;
    }
    memcpy(&d[row+1][1], &dbuf[1], band*sizeof(int16_t));
    memcpy(&p[row+1][1], &pbuf[1], band*sizeof(int16_t));
  } // for(row=band;row<=len1+len2;row+=2) (Banded body)

  // Fill out bottom wedge
  for(row=len1+len2-band+1;row<=len1+len2;row+=2) { // Using fact that len1+len2=even, starting on odd row now
    // Fill out odd row
    col_min = center-(1+len1+len2-row)/2;  // includes ends-free part
    col_max = center+(1+len1+len2-row)/2-1;  // includes ends-free part
    
    for(col=col_min,i=len1-1,j=row-len1-1;col<1+col_max;col++,i--,j++) { // includes ends-free part, could make limit based on j<=len2
      left = d[row-1][col] + (i==(len1-1) ? 0 : gap_p);
      diag = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
      up = d[row-1][col+1] + (j==(len2-1) ? 0 : gap_p);
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      d[row][col] = entry; // Vectorizes if this points to another array
      p[row][col] = pentry;
    }
    
    col_min++; // Because avoided ends-free previously
    // Fill out even row (row+1)
    for(col=col_min,i=len1-1,j=row-len1;col<1+col_max;col++,i--,j++) {
      left = d[row][col-1] + (i==(len1-1) ? 0 : gap_p);
      diag = d[row-1][col] + (s1[i] == s2[j] ? match : mismatch);
      up = d[row][col] + (j==(len2-1) ? 0 : gap_p);
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      d[row+1][col] = entry; // Vectorizes if this points to another array
      p[row+1][col] = pentry;
    }
    
    
  }

/*
  for(row=len1+len2-10;row<=len1+len2;row++) {
    for(col=0;col<=center+band/2+1;col++) {
      Rprintf("%i,",p[row][col]);
    }
    Rprintf("\n");
  }
  Rprintf("\n");

  for(row=len1+len2-10;row<=len1+len2;row++) {
    for(col=0;col<=center+band/2+1;col++) {
      Rprintf("%+04i,",d[row][col]);
    }
    Rprintf("\n");
  }
  Rprintf("\n");
*/

  char al0[2*SEQLEN+1];
  char al1[2*SEQLEN+1];
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;
  
  while ( i > 0 || j > 0 ) {
    ///v Rprintf("ij=(%i,%i), rc=(%i,%i), p[][]=%i\n", i,j,i+j,center-(i-j), p[i+j][center-(i-j)]);
    switch ( p[i+j][1 + (band-i+j)/2] ) {
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
        ///v Rprintf("ij=(%i,%i), rc=(%i,%i), p[][]=%i\n", i,j,i+j,center-(i-j), p[i+j][center-(i-j)]);
        Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Failed memory allocation.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Failed memory allocation.");
  
  // Reverse the alignment strings (since traced backwards).
  for (index=0;index<len_al;index++) {
    al[0][index] = al0[len_al-index-1];
    al[1][index] = al1[len_al-index-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  
  return al;
}

// [[Rcpp::export]]
Rcpp::CharacterVector C_nwvec(std::string s1, std::string s2, Rcpp::NumericMatrix score, int gap_p, int band) {
  size_t i, j;
  char *seq1, *seq2;
  char **al;

  seq1 = (char *) malloc(s1.size()+1); //E
  seq2 = (char *) malloc(s2.size()+1); //E
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  strcpy(seq1, s1.c_str());
  strcpy(seq2, s2.c_str());
  nt2int(seq1, seq1);
  nt2int(seq2, seq2);
  
  int16_t c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = (int16_t) score(i,j);
    }
  }
  
  al = nwalign_endsfree_vectorized(seq1, seq2, c_score, (int16_t) gap_p, (size_t) band);

  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);

  Rcpp::CharacterVector rval;
  rval.push_back(std::string(al[0]));
  rval.push_back(std::string(al[1]));
  
  free(seq1);
  free(seq2);
  free(al[0]);
  free(al[1]);
  return(rval);
}
