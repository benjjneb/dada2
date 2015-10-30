#include <Rcpp.h>
#include "dada.h"
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

char **nwalign_endsfree_vectorized(char *s1, char *s2, int score[4][4], int gap_p, int band) {
  int row, col;
  int i,j;
  int len1 = strlen(s1);
  int len2 = strlen(s2);
  int d[800][800]; // d: DP matrix
  int p[800][800]; // backpointer matrix with 1 for diagonal, 2 for left, 3 for up.
  // Too big for the Effing stack at 2*SEQLEN+1 x 2*SEQLEN+1
  // HAVE TO FIX THIS!!
  int dbuf[2*SEQLEN+1];
  int pbuf[2*SEQLEN+1];
  int diag, left, up, entry, pentry;
//  int excess1, excess2;
  
  int center;
  
//  excess1 = len1>len2 ? len1-len2 : 0;
//  excess2 = len2>len1 ? len2-len1 : 0;

  // assuming len1=len2 for now
  if(band>=0 && band<len1) {
    center=band+1;
  } else {
    center=len1+1;
  }
  
  // Fill out "left" "column" of d, p.
  for (row=0,col=center; col>=0; row++, col--) {
    d[row][col] = 0; // ends-free gap
    p[row][col] = 3;
  }
  
  // Fill out "top" "row" of d, p.
  for (row=0,col=center; col<=len2+center; row++, col++) { // perf: can cut off at band
    d[row][col] = 0; // ends-free gap
    p[row][col] = 2;
  }
    
  // Fill out band boundaries
  if(band>=0 && band<len1) {
    for(row=0;row<=len1+len2;row++) {
      d[row][0] = -9999;
    }
  }
  if(band>=0 && band<len2) {
    for(row=0;row<=len1+len2;row++) {
      d[row][center+band+1] = -9999;
    }
  }
  
  // Fill out top wedge (Row 1 taken care of by ends-free)
  for(row=2;row<center;row++) {
    for(col=center-(row-1);col<=center+(row-1);col++) { // ONLY HALF THE ENTRIES NEED TO BE LOOPED OVER
      left = d[row-1][col-1] + gap_p;
      diag = d[row-2][col] + score[s1[(row+(col-center))/2 - 1]-1][s2[(row-(col-center))/2 - 1]-1];
      up = d[row-1][col+1] + gap_p;
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      d[row][col] = entry; // Vectorizes if this points to another array
      p[row][col] = pentry;
    }
  }

  // Fill out banded body (and bottom wedge)
  ////!!! MAKING MATCH/MISMATCH HERE
  int match = score[0][0];
  int mismatch = score[0][2];
  for(row=center;row<=len1+len2;row++) {
//    #pragma clang loop vectorize(enable)
    for(col=1;col<=center+band;col++) { // ONLY HALF THE ENTRIES NEED TO BE LOOPED OVER
      left = d[row-1][col-1] + gap_p;
//      diag = d[row-2][col] + score[s1[(row+(col-center))/2 - 1]-1][s2[(row-(col-center))/2 - 1]-1];
      diag = d[row-2][col] + s1[col] == s2[col] ? match : mismatch;
      up = d[row-1][col+1] + gap_p;
      
      entry = up >= left ? up : left;
      pentry = up >= left ? 3 : 2;
      pentry = entry >= diag ? pentry : 1;
      entry = entry >= diag ? entry : diag;

      dbuf[col] = entry; // Vectorizes if this points to another array
      pbuf[col] = pentry;
    }
    memcpy(&d[row][1], &dbuf[1], (2*band+1)*sizeof(int));
    memcpy(&p[row][1], &pbuf[1], (2*band+1)*sizeof(int));
  }

/*
  for(row=0;row<10;row++) {
    for(col=center-6;col<=center+6;col++) {
      Rprintf("%i,",p[row][col]);
    }
    Rprintf("\n");
  }
  Rprintf("\n");

  for(row=0;row<10;row++) {
    for(col=center-6;col<=center+6;col++) {
      Rprintf("%i,",d[row][col]);
    }
    Rprintf("\n");
  }
  Rprintf("\n");
*/

  char al0[2*SEQLEN+1];
  char al1[2*SEQLEN+1];
  
  // Trace back over p to form the alignment.
  unsigned int len_al = 0;
  i = len1;
  j = len2;
  
  while ( i > 0 || j > 0 ) {
    ///v Rprintf("ij=(%i,%i), rc=(%i,%i), p[][]=%i\n", i,j,i+j,center-(i-j), p[i+j][center-(i-j)]);
    switch ( p[i+j][center-(i-j)] ) {
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
  int i, j;
  char *seq1, *seq2;
  char **al;

  seq1 = (char *) malloc(s1.size()+1); //E
  seq2 = (char *) malloc(s2.size()+1); //E
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  strcpy(seq1, s1.c_str());
  strcpy(seq2, s2.c_str());
  nt2int(seq1, seq1);
  nt2int(seq2, seq2);
  
  int c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = (int) score(i,j);
    }
  }
  
  al = nwalign_endsfree_vectorized(seq1, seq2, c_score, gap_p, band);

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
