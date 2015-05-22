#include <Rcpp.h>
#include "dada.h"
#include "align.h"
using namespace Rcpp;
// [[Rcpp::interfaces(r,cpp)]]

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dada_nw(std::string s1, std::string s2, Rcpp::NumericMatrix score, int gap, int band, bool FAST_ALIGN) {
  int i,j,nrow, ncol;
  char seq1[SEQLEN];
  char seq2[SEQLEN];

  long long_gap = (long) gap;
  long * nwscore;
  long * nwdiff;
  long * nwgaps;
  long * nwindels;
  long * nwalignmentlength;
  char ** nwalignment;
  nwinfo_s * nw = nw_init();

  nt2int(seq1, s1.c_str());
  nt2int(seq2, s2.c_str());

  // Copy score into a C style array
  nrow = score.nrow();
  ncol = score.ncol();
  if(nrow != 4 || ncol != 4) {
    Rcpp::Rcout << "C: Score matrix malformed:" << nrow << ", " << ncol << "\n";
    return 0;
  }
  double c_score[4][4];
  long long_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = score(i,j);
      long_score[i][j] = (long) score(i,j);
    }
  }

  if(FAST_ALIGN) {
    nw_align(seq1, seq1 + strlen(seq1), seq2, seq2 + strlen(seq2),
              long_score,
              0L, -long_gap, 0L, // q gap opens
              0L, -long_gap, 0L, // t gap opens
              0L, -long_gap, 0L, // q gap extends
              0L, -long_gap, 0L, // t gap extends
              nwscore,
              nwdiff,
              nwgaps,
              nwindels,
              nwalignmentlength,
              nwalignment,
              0L, 0L, // vsearch diagnostics
              nw);
  } else {
    nwalign_endsfree(seq1, seq2, c_score, gap, band);
  }
  return 1;
}

char ntoi(char nt) {
  char rchar = '\0';
  switch(nt) {
    case 'A':
      rchar = 1;
      break;
    case 'C':
      rchar = 2;
      break;
    case 'G':
      rchar = 3;
      break;
    case 'T':
      rchar = 4;
      break;
    case 'N':
      rchar = 5;
      break;
    case '-':
      rchar = 6;
      break;
    case '\0':
      rchar = '\0';
      break;
    default:
      rchar = nt;
      break;
  }
  return(rchar);
}

