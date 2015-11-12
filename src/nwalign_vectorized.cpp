#include <Rcpp.h>
#include "dada.h"
using namespace Rcpp;

void dploop_vec(int16_t *__restrict__ ptr_left, int16_t *__restrict__ ptr_diag, int16_t *__restrict__ ptr_up, int16_t *__restrict__ d, int16_t *__restrict__ p, int16_t gap_p, size_t n) {
  int16_t left, diag, up, entry, pentry;
  size_t i = 0;
  
  while(i<n) {
    left = *ptr_left + gap_p;
    diag = *ptr_diag;
    up = *ptr_up + gap_p;
    
    entry = up >= left ? up : left;
    pentry = up >= left ? 3 : 2;
    pentry = entry >= diag ? pentry : 1;
    entry = entry >= diag ? entry : diag;

    *d = entry; // Vectorizes if this points to another array
    *p = pentry;
    
    ptr_left++;
    ptr_diag++;
    ptr_up++;
    d++;
    p++;
    i++;
  }
}

char **nwalign_endsfree_vectorized(char *s1, char *s2, int16_t score[4][4], int16_t gap_p, size_t band) {
  size_t row, col;
  size_t i,j;
  size_t len1 = strlen(s1);
  size_t len2 = strlen(s2);
  int16_t d[800][800]; // d: DP matrix
  int16_t p[800][800]; // backpointer matrix with 1 for diagonal, 2 for left, 3 for up.
  // Too big for the Effing stack at 2*SEQLEN+1 x 2*SEQLEN+1
  int16_t diag_buf[SEQLEN];
  int16_t diag, left, up, entry, pentry, d_free;
  size_t center;
  size_t col_min, col_max;
  size_t i_max, j_min;
  int16_t *ptr_left, *ptr_diag, *ptr_up;
  
  // Simplifying to match/mismatch alignment
  int16_t match = score[0][0];
  int16_t mismatch = score[0][2];

  // assuming len1=len2 for now
  center=1 + (band+1)/2;
  
  // Initialize top to -9999
  for(row=0;row<band;row++) {
    for(col=0;col<center+band+2;col++) {
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
  if(band%2 == 0) { // even band
    for(row=0;row<len1+len2+1;row++) {
      d[row][0] = -9999;
    }
  
    for(row=0;row<len1+len2+1;row+=2) {
      d[row][band+2] = -9999;
      d[row+1][band+1] = -9999;
      d[row+1][band+2] = -9999;
    }
  } else { // odd band
    for(row=0;row<len1+len2+1;row+=2) {
      d[row][0] = -9999;
      d[row][1] = -9999;
      d[row+1][0] = -9999;
    }    
    for(row=0;row<len1+len2+1;row++) {
      d[row][band+2] = -9999;
    }
  }
  
  // Fill out top wedge (Row 1 taken care of by ends-free)
  row = 2;
  col_min = center; // Do not fill out the ends-free cells
  col_max = center;
  i_max = 0;
  j_min = 0;
  while(row <= band) {
    // Fill out even row
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_left = &d[row-1][col_min-1];
    ptr_diag = &diag_buf[col_min];
    ptr_up = &d[row-1][col_min];
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row][col_min], &p[row][col_min], gap_p, col_max-col_min+1);
    
    col_min--;
    i_max++;
    if(++row > band) break;
    
    // Fill out odd row
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_left = &d[row-1][col_min];
    ptr_diag = &diag_buf[col_min];
    ptr_up = &d[row-1][col_min+1];
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row][col_min], &p[row][col_min], gap_p, col_max-col_min+1);
    
    col_max++;
    i_max++;
    row++;
  }

  // ----- FILL OUT BANDED BODY ------
  // Initialize indexing values
  row=band+1;
  if(band % 2 == 0) { // even
    col_min=1;
    col_max=band;
  } else {
    col_min=2;
    col_max=band+1;
  }
  j_min = 0; // First row out of band, can still address first (0th) element of string 
  i_max = band - 1; // Remember, i+j = row - 2 (because of zero-indexing of string)
  
  // Loop over banded region
  while(row<=len1+len2-band) {// Using fact that len1+len2=even. This loop covers an even number of rows (band+1...len1+len2-band)
    // Fill out short row (different parity than band, n_elements=band)
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_diag = &diag_buf[col_min];
    if(row%2 == 0) { // even row
      ptr_left = &d[row-1][col_min-1];
      ptr_up = &d[row-1][col_min];
    } else { // odd row
      ptr_left = &d[row-1][col_min];
      ptr_up = &d[row-1][col_min+1];
    }
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row][col_min], &p[row][col_min], gap_p, col_max-col_min+1);
    
    if(band%2 == 0) { // even band (and odd row)
      col_max++;
    } else { // odd band (and even row)
      col_min--;
    }
    i_max++;
    row++;
    
    // Fill out long row (same parity as band, n_elements=band+1)  
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_diag = &diag_buf[col_min];
    if(row%2 == 0) { // even row
      ptr_left = &d[row-1][col_min-1];
      ptr_up = &d[row-1][col_min];
    } else { // odd row
      ptr_left = &d[row-1][col_min];
      ptr_up = &d[row-1][col_min+1];
    }
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row][col_min], &p[row][col_min], gap_p, col_max-col_min+1);
    
    if(band%2 == 0) { // even band (and even row)
      col_max--;
    } else { // odd band (and odd row)
      col_min++;
    }
    j_min++;
    row++;
  }
  
  // ----- FILL OUT BOTTOM WEDGE ------
  // Initialize starting values
  row = len1+len2-band+1; // Always opposite parity of band (len1=len2)
  if(band%2 == 0) { // even band, odd row
    col_min = 1;
    col_max = band;
  } else { // odd band, even row
    col_min = 2;
    col_max = band+1;
  }
  i_max = len1-1; // reached end of seq1
  j_min = row - i_max - 2; // i+j=row-2

  // Loop over lower wedge.......
  while(row <= len1+len2) {
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[row-2][col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_diag = &diag_buf[col_min];
    if(row%2 == 0) { // even row
      ptr_left = &d[row-1][col_min-1];
      ptr_up = &d[row-1][col_min];
    } else { // odd row
      ptr_left = &d[row-1][col_min];
      ptr_up = &d[row-1][col_min+1];
    }
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row][col_min], &p[row][col_min], gap_p, col_max-col_min+1);
    
    // Redo the ends-free cells
    // Left column
    if(row%2==0) {
      d_free = d[row-1][col_min-1];
    } else {
      d_free = d[row-1][col_min];
    }
    if(d_free > d[row][col_min]) { // ends-free gap is better
      d[row][col_min] = d_free;
      p[row][col_min] = 2;
    } else if(d_free == d[row][col_min] && p[row][col_min] == 1) { // left gap takes precedence over diagonal move (for consistency)
      p[row][col_min] = 2;      
    }
    // Right column
    if(row%2==0) {
      d_free = d[row-1][col_max];
    } else {
      d_free = d[row-1][col_max+1];
    }
    if(d_free > d[row][col_max]) { // ends-free gap is better
      d[row][col_max] = d_free;
      p[row][col_max] = 3;
    } else if(d_free == d[row][col_max] && p[row][col_max] != 3) { // up gap takes precedence over left or diagonal move (for consistency)
      p[row][col_max] = 3;
    }
    
    // Update indices
    if(row%2 == 0) { // even row
      col_max--;
    } else { // odd row
      col_min++;      
    }
    j_min++;
    row++;
  }

/*  for(row=len1+len2-band;row<=len1+len2;row++) {
    for(col=0;col<=2*band;col++) {
      Rprintf("%04i ", d[row][col]);
    }
    Rprintf("\n");
  }

  for(row=len1+len2-band;row<=len1+len2;row++) {
    for(col=0;col<=2*band;col++) {
      Rprintf("%i ", p[row][col]);
    }
    Rprintf("\n");
  }*/

  char al0[2*SEQLEN+1];
  char al1[2*SEQLEN+1];
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;
  
  while ( i > 0 || j > 0 ) {
    ///v Rprintf("ij=(%i,%i), rc=(%i,%i), p[][]=%i\n", i,j,i+j,center-(i-j), p[i+j][center-(i-j)]);
///    switch ( p[i+j][1 + (band-i+j)/2] ) {
    switch ( p[i+j][(2*center-i+j)/2] ) {
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
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
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
  
  int c_score_int[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = (int) score(i,j);
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
