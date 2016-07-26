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

void parr(int16_t *arr, int nrow, int ncol) {
  int col, row;
  for(row=0;row<nrow;row++) {
    for(col=0;col<ncol;col++) {
      Rprintf("%05d\t", arr[row*ncol+col]);
    }
    Rprintf("\n");
  }
}

char **nwalign_vectorized2(char *s1, char *s2, int16_t match, int16_t mismatch, int16_t gap_p, int16_t end_gap_p, int band) {
  size_t row, col, ncol, nrow, foo;
  size_t i,j;
  size_t len1, len2;
  int16_t d_free;
  size_t start_col, end_col;
  size_t col_min, col_max, even;
  size_t i_max, j_min;
  int16_t *ptr_left, *ptr_diag, *ptr_up, *ptr_index, *ptr_d, *ptr_p;
  bool swap = false;
  char *ptr_char;

  len1 = strlen(s1);
  len2 = strlen(s2);
  if(len1 > len2) { // Ensure s1 is the shorter sequences
    ptr_char = s1;
    s1 = s2;
    s2 = ptr_char;
    swap = true;
    foo = len1;
    len1 = len2;
    len2 = foo;
  }
  if(band < 0) { band = len2; }
  
  // Allocate the DP matrices
  start_col = 1 + (1+(band<len1 ? band : len1))/2;
  end_col = start_col + (len2-len1)/2;
  //  ncol = 3 + (len1+len2+1)/2; // 3 = left boundary + center + right boundary !!!
  ncol = 2 + start_col + ((len2-len1+band)<len2 ? (len2-len1+band) : len2)/2;
  nrow = len1 + len2 + 1;
  int16_t *d = (int16_t *) malloc(ncol * nrow * sizeof(int16_t));
  int16_t *p = (int16_t *) malloc(ncol * nrow * sizeof(int16_t));
  int16_t *diag_buf = (int16_t *) malloc(ncol * sizeof(int16_t));
  if (d == NULL || p == NULL || diag_buf == NULL)  Rcpp::stop("Memory allocation failed.");
  
  for(row=0;row<nrow;row++) {
    for(col=0;col<ncol;col++) {
      d[row*ncol+col] = -9999;
      p[row*ncol+col] = 0;
    }
  }
  
  // Fill out starting point
  d[start_col] = 0;
  p[start_col] = 0; // Should never be queried
  
  // Fill out "left" "column" of d, p.
  row=1;
  col=start_col-1;
  ptr_d=&d[ncol]; // start of 1st row
  ptr_p=&p[ncol];
  d_free = end_gap_p; // end-gap score
  while(row < (1 + (band < len1 ? band : len1))) {
    ptr_d[col] = d_free; // end gap
    ptr_p[col] = 3;
    if(row%2==0) {
      col--;
    }
    row++;
    d_free += end_gap_p;
    ptr_d += ncol;
    ptr_p += ncol;
  }
  
  // Fill out "top" "row" of d, p.
  row=1;
  col=start_col;
  ptr_d=&d[ncol]; // start of 1st row
  ptr_p=&p[ncol];
  d_free = end_gap_p; // end-gap score
  while(row < (1+(band+len2-len1 < len2 ? band+len2-len1 : len2))) {
    ptr_d[col] = d_free;
    ptr_p[col] = 2;
    if(row%2 == 1) {
      col++;
    }
    row++;
    d_free += end_gap_p;
    ptr_d += ncol;
    ptr_p += ncol;
  }
  
/*  Rprintf("D MATRIX:.....\n");
  parr(d, len2+1, ncol);
  Rprintf("\n");
  Rprintf("P MATRIX:.....\n");
  parr(p, len2+1, ncol);
  Rprintf("\n");*/
  
    // Fill out DP matrix (Row 0/1 taken care of by ends-free)
  row = 2;
  col_min = start_col; // Do not fill out the ends-free cells
  col_max = start_col;
  i_max = 0; // 1st nt
  j_min = 0; // 1st nt
  even = 1; // TRUE
  
  while(row <= (len1+len2)) {
      // Fill out row
    for(col=col_min,i=i_max,j=j_min;col<(1+col_max);col++,i--,j++) {
      diag_buf[col] = d[(row-2)*ncol + col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_left = &d[(row-1)*ncol + col_min-even];
    ptr_diag = &diag_buf[col_min];
    ptr_up = &d[(row-1)*ncol + col_min+1-even];
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    
    // Recalculate ends-free for lower boundary
    // Redo the ends-free cells
    // Left column
    if(row > len1 && (end_gap_p > gap_p)) {
      d_free = ptr_left[0] + end_gap_p; // first cell of left array
      if(d_free > d[row*ncol + col_min]) { // ends-free gap is better
        d[row*ncol + col_min] = d_free;
        p[row*ncol + col_min] = 2;
      } else if(d_free == d[row*ncol + col_min] && p[row*ncol + col_min] == 1) { // left gap takes precedence over diagonal move (for consistency)
        p[row*ncol + col_min] = 2;
      }
    }
    // Right column
    if(row > len2 && (end_gap_p > gap_p)) {
      d_free = ptr_up[col_max-col_min] + end_gap_p; // last cell of up array
      if(d_free > d[row*ncol + col_max]) { // ends-free gap is better
        d[row*ncol + col_max] = d_free;
        p[row*ncol + col_max] = 3;
      } else if(d_free == d[row*ncol + col_max] && p[row*ncol + col_max] != 3) { // up gap takes precedence over left or diagonal move (for consistency)
        p[row*ncol + col_max] = 3;
      }
    }
    
    // Update indices
    if(row==(band<len1 ? band : len1)) { // Offset to the boundary, didn't do in top wedge cause pre-filled
      col_min--;
      i_max++;
      j_min--;
    }
    if(row < band && row < len1) { // upper tri for seq1
      if(even) { col_min--; }
      i_max++;
    } else if(i_max < len1-1) { // banded area
      if(even) { col_min--; i_max++; }
      else { col_min++; j_min++; }
    } else { // lower tri for seq1
      if(!even) { col_min++; }
      j_min++;
    }

    if(row == (band+len2-len1 < len2 ? band+len2-len1 : len2)) {
      // Offset to the boundary, didn't do in top wedge cause pre-filled
      col_max++;
    }
    if(row<(band+len2-len1 < len2 ? band+len2-len1 : len2)) {
      if(!even) { col_max++; }
    } else if((row+1)/2 + col_max - start_col < len2) { // "j_max" (1-index) < len2
      if(even) { col_max--; }
      else { col_max++; }
    } else {
      if(even) { col_max--; }
    }
    
    row++;
    even = 1 - even;
  }

/*  Rprintf("D MATRIX:.....\n");
  parr(d, len1+len2+1, ncol);
  Rprintf("\n");
  Rprintf("P MATRIX:.....\n");
  parr(p, len1+len2+1, ncol);
  Rprintf("\n");*/

  char *al0 = (char *) malloc((nrow+1) * sizeof(char));
  char *al1 = (char *) malloc((nrow+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;
  
  while ( i > 0 || j > 0 ) {
    switch ( p[(i+j)*ncol + (2*start_col+j-i)/2] ) {
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
  
  // Free DP objects
  free(d);
  free(p);
  free(diag_buf);
  
  // Return to input ordering
  if(swap) {
    ptr_char = al0;
    al0 = al1;
    al1 = ptr_char;
  }
  
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
  
  free(al0);
  free(al1);
  
  return al;
}

char **nwalign_endsfree_vectorized(char *s1, char *s2, int16_t match, int16_t mismatch, int16_t gap_p, int band) {
  size_t row, col, ncol, nrow;
  size_t i,j;
  size_t len1 = strlen(s1);
  size_t len2 = strlen(s2);
  int16_t d_free;
  size_t center;
  size_t col_min, col_max;
  size_t i_max, j_min;
  int16_t *ptr_left, *ptr_diag, *ptr_up, *ptr_index;
  
  // Require same length sequences
  if(len1 != len2) {
    Rcpp::stop("Vectorized alignment currently only implemented for same length sequences.");
  }
  
  // Deal with band possibilities
  if(band == 0) {
    Rcpp::stop("Vectorized alignment not currently implemented for band = 0.");
  }
  if(band<0 || band>len1) {
    band = len1;
  }
  
  // Allocate the DP matrices
  center=1 + (band+1)/2;
  ncol = 3 + band;
  nrow = len1 + len2 + 2;
  int16_t *d = (int16_t *) malloc(ncol * nrow * sizeof(int16_t));
  int16_t *p = (int16_t *) malloc(ncol * nrow * sizeof(int16_t));
  int16_t *diag_buf = (int16_t *) malloc(ncol * sizeof(int16_t));
  if (d == NULL || p == NULL || diag_buf == NULL)  Rcpp::stop("Memory allocation failed.");
  
  for(row=0;row<nrow;row++) {
    for(col=0;col<ncol;col++) {
      d[row*ncol+col] = -9999;
      p[row*ncol+col] = 0;
    }
  }
  
  // Fill out starting point
  d[center] = 0;
  p[center] = 0; // Should never be queried
  
  // Fill out "left" "column" of d, p.
  row=1;
  col=center-1;
  while(row<band+1) {
    d[row*ncol + col] = 0; // ends-free gap
    p[row*ncol + col] = 3;
    if(row%2==0) {
      col--;
    }
    row++;
  }
  
  // Fill out "top" "row" of d, p.
  row=1;
  col=center;
  while(row<band+1) {
    d[row*ncol + col] = 0;
    p[row*ncol + col] = 2;
    if(row%2 == 1) {
      col++;
    }
    row++;
  }
  
/*  Rprintf("D MATRIX:.....\n");
  parr(d, len2+1, ncol);
  Rprintf("\n");
  Rprintf("P MATRIX:.....\n");
  parr(p, len2+1, ncol);
  Rprintf("\n"); */
  
  // Fill out band boundaries
  if(band%2 == 0) { // even band
    for(row=0,ptr_index=&d[0];row<len1+len2+1;row++,ptr_index+=ncol) {
      *ptr_index = -9999;
    }
    
    ptr_index = &d[band+2];
    for(row=0;row<len1+len2+1;row+=2) {
      *ptr_index = -9999; // even row
      ptr_index += (ncol-1);
      *ptr_index = -9999; // odd row
      ptr_index += (ncol+1);
    }
  } else { // odd band
    ptr_index = &d[1];
    for(row=0;row<len1+len2+1;row+=2) {
      *ptr_index = -9999;
      ptr_index += (ncol-1);
      *ptr_index = -9999;
      ptr_index += (ncol+1);
    }
    for(row=0,ptr_index=&d[band+2];row<len1+len2+1;row++,ptr_index+=ncol) {
      *ptr_index = -9999;
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
      diag_buf[col] = d[(row-2)*ncol + col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_left = &d[(row-1)*ncol + col_min-1];
    ptr_diag = &diag_buf[col_min];
    ptr_up = &d[(row-1)*ncol + col_min];
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    
    col_min--;
    i_max++;
    if(++row > band) break;
    
    // Fill out odd row
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = d[(row-2)*ncol + col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_left = &d[(row-1)*ncol + col_min];
    ptr_diag = &diag_buf[col_min];
    ptr_up = &d[(row-1)*ncol + col_min+1];
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    
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
  while(row<len1+len2-band+1) {// Using fact that len1+len2=even. This loop covers an even number of rows (band+1...len1+len2-band)
    // Fill out short row (different parity than band, n_elements=band)
    ptr_diag = &d[(row-2)*ncol];
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = ptr_diag[col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_diag = &diag_buf[col_min];
    ptr_left = &d[(row-1)*ncol + col_min-1]; // even row, "left" is up-left
    if(row%2) { // odd row, "left" is straight up
      ptr_left++;
    }
    ptr_up = ptr_left + 1;
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    
    if(band%2 == 0) { // even band (and odd row)
      col_max++;
    } else { // odd band (and even row)
      col_min--;
    }
    i_max++;
    row++;
    
    // Fill out long row (same parity as band, n_elements=band+1)  
    ptr_diag = &d[(row-2)*ncol];
    for(col=col_min,i=i_max,j=j_min;col<1+col_max;col++,i--,j++) {
      diag_buf[col] = ptr_diag[col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_diag = &diag_buf[col_min];
    ptr_left = &d[(row-1)*ncol + col_min-1]; // even row, "left" is up-left
    if(row%2) { // odd row, "left" is straight up
      ptr_left++;
    }
    ptr_up = ptr_left + 1;
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    
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
      diag_buf[col] = d[(row-2)*ncol + col] + (s1[i] == s2[j] ? match : mismatch);
    }
    ptr_diag = &diag_buf[col_min];
    if(row%2 == 0) { // even row
      ptr_left = &d[(row-1)*ncol + col_min-1];
      ptr_up = &d[(row-1)*ncol + col_min];
    } else { // odd row
      ptr_left = &d[(row-1)*ncol + col_min];
      ptr_up = &d[(row-1)*ncol + col_min+1];
    }
    dploop_vec(ptr_left, ptr_diag, ptr_up, &d[row*ncol + col_min], &p[row*ncol + col_min], gap_p, col_max-col_min+1);
    
    // Redo the ends-free cells
    // Left column
    if(row%2==0) {
      d_free = d[(row-1)*ncol + col_min-1];
    } else {
      d_free = d[(row-1)*ncol + col_min];
    }
    if(d_free > d[row*ncol + col_min]) { // ends-free gap is better
      d[row*ncol + col_min] = d_free;
      p[row*ncol + col_min] = 2;
    } else if(d_free == d[row*ncol + col_min] && p[row*ncol + col_min] == 1) { // left gap takes precedence over diagonal move (for consistency)
      p[row*ncol + col_min] = 2;
    }
    // Right column
    if(row%2==0) {
      d_free = d[(row-1)*ncol + col_max];
    } else {
      d_free = d[(row-1)*ncol + col_max+1];
    }
    if(d_free > d[row*ncol + col_max]) { // ends-free gap is better
      d[row*ncol + col_max] = d_free;
      p[row*ncol + col_max] = 3;
    } else if(d_free == d[row*ncol + col_max] && p[row*ncol + col_max] != 3) { // up gap takes precedence over left or diagonal move (for consistency)
      p[row*ncol + col_max] = 3;
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
  
/*  Rprintf("D MATRIX:.....\n");
  parr(d, len1+len2+1, ncol);
  Rprintf("\n");
  Rprintf("P MATRIX:.....\n");
  parr(p, len1+len2+1, ncol);
  Rprintf("\n"); */
  
  char *al0 = (char *) malloc((nrow+1) * sizeof(char));
  char *al1 = (char *) malloc((nrow+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;
  
  while ( i > 0 || j > 0 ) {
    switch ( p[(i+j)*ncol + (2*center-i+j)/2] ) {
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
  
  // Free DP objects
  free(d);
  free(p);
  free(diag_buf);
  
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
  
  free(al0);
  free(al1);
  
  return al;
}

// [[Rcpp::export]]
Rcpp::CharacterVector C_nwvec(std::string s1, std::string s2, int16_t match, int16_t mismatch, int16_t gap_p, int band, bool endsfree, bool test) {
  char *seq1, *seq2;
  char **al;

  seq1 = (char *) malloc(s1.size()+1); //E
  seq2 = (char *) malloc(s2.size()+1); //E
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  strcpy(seq1, s1.c_str());
  strcpy(seq2, s2.c_str());
  nt2int(seq1, seq1);
  nt2int(seq2, seq2);
  
  if(test) {
    if(endsfree) {
      al = nwalign_vectorized2(seq1, seq2, match, mismatch, gap_p, 0, (size_t) band);
    } else {
      al = nwalign_vectorized2(seq1, seq2, match, mismatch, gap_p, gap_p, (size_t) band);
    }
  } else {
    if(!endsfree) Rprintf("Standard vectorized aligner is endsfree-only.\n");
    al = nwalign_endsfree_vectorized(seq1, seq2, match, mismatch, gap_p, (size_t) band);
  }

  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);

  Rcpp::CharacterVector rval;
  rval.push_back(std::string(al[0]));
  rval.push_back(std::string(al[1]));
  
  free(seq1);
  free(seq2);
  free(al[0]);
  free(al[1]);
  free(al);
  return(rval);
}
