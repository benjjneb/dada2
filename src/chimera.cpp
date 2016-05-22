#include "dada.h"
#include <Rcpp.h>
using namespace Rcpp;

int get_ham_endsfree(const char *seq1, const char *seq2, int len);

//------------------------------------------------------------------
// Determines whether sq is a perfect bimera of some combination from pars.
// 
// @param sq The query DNA sequence.
// @param pars Potential bimeric "parents".
// @param max_shift A \code{integer(1)} of the maximum alignment shift allowed.
// 
// [[Rcpp::export]]
bool C_is_bimera(std::string sq, std::vector<std::string> pars, bool allow_one_off, int min_one_off_par_dist, Rcpp::NumericMatrix score, int gap_p, int max_shift) {
  // For now only finding perfect bimeras
  int i, j, left, right, left_oo, right_oo, pos, len, max_len;
  // Make  c-style 2d score array
  int c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = (int) score(i,j);
    }
  }
  
  // Make integer-ized c-style sequence strings
  char *seq1 = (char *) malloc(sq.size()+1); //E
  // Get max length of the parent sequences
  max_len=0;
  for(i=0;i<pars.size();i++) {
    len = pars[i].size();
    if(len>max_len) { max_len = len; }
  }
  char *seq2 = (char *) malloc(max_len+1); //E
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  nt2int(seq1, sq.c_str());
  
  char **al;
  int max_left=0, max_right=0;
  int oo_max_left=0, oo_max_right=0, oo_max_left_oo=0, oo_max_right_oo=0;
  
  bool rval = false;
  for(i=0;i<pars.size() && rval==false;i++) {
    nt2int(seq2, pars[i].c_str());
    al = nwalign_endsfree(seq1, seq2, c_score, gap_p, max_shift);  // Remember, alignments must be freed!
    len = strlen(al[0]);
    
    pos=0; left=0;
    while(al[0][pos] == 6 && pos<len) {
      pos++; // Scan in until query starts
    }
    while(al[1][pos] == 6 && pos<max_shift) {
      pos++; left++; // Credit as ends-free coverage until parent starts
    }
    while(pos<len && al[0][pos] == al[1][pos]) {
      pos++; left++; // Credit as covered until a mismatch
    }
    if(allow_one_off) {
      // Step forward, and credit to one-off further matches (and this mismatch if not a gap)
      left_oo = left;
      pos++;
      if(pos<len && al[0][pos] != 6) { left_oo++; }
      while(pos<len && al[0][pos] == al[1][pos]) {
        pos++; left_oo++;
      }
    }
    
    pos=len-1; right=0;
    while(al[0][pos] == 6 && pos >= 0) { 
      pos--;
    }
    while(al[1][pos] == 6 && pos>+(len-max_shift)) {
      pos--; right++;
    }
    while(pos>=0 && al[0][pos] == al[1][pos]) {
      pos--; right++;
    }
    if(allow_one_off) {
      // Step forward, and credit to one-off further matches (and this mismatch if not a gap)
      right_oo = right;
      pos--;
      if(pos>=0 && al[0][pos] != 6) { right_oo++; }
      while(pos>=0 && al[0][pos] == al[1][pos]) {
        pos--; right_oo++;
      }
    }
    
    if((left+right) >= sq.size()) { // Toss id/pure-shift/internal-indel "parents"
      continue;
    }
    if(left > max_left) { max_left=left; }
    if(right > max_right) { max_right=right; }
    
    // Need to evaluate whether parents are allowed for one-off models
    if(allow_one_off && get_ham_endsfree(al[0], al[1], len) >= min_one_off_par_dist) {
      if(left > oo_max_left) { oo_max_left=left; }
      if(right > oo_max_right) { oo_max_right=right; }
      if(left_oo > oo_max_left_oo) { oo_max_left_oo=left_oo; }
      if(right_oo > oo_max_right_oo) { oo_max_right_oo=right_oo; }
    }
    
    // Evaluate, and break if found bimeric model
    if((max_right+max_left)>=sq.size()) {
      rval = true;
    }
    if(allow_one_off) {
      if((oo_max_left+oo_max_right_oo)>=sq.size() || (oo_max_left_oo+oo_max_right)>=sq.size()) {
        rval=true;
      }
    }
    free(al[0]);
    free(al[1]);
    free(al);
  }
  
  free(seq1);
  free(seq2);
  return(rval);
}

// Internal function to get hamming distance between aligned seqs
//  without counting end gaps
int get_ham_endsfree(const char *seq1, const char *seq2, int len) {
  int i, j, pos, ham;
  bool gap1, gap2;
  // Find start of internal part of alignment
  i=0;
  gap1 = (seq1[i]==6); 
  gap2 = (seq2[i]==6);
  while(gap1 || gap2) {
    i++;
    gap1 = (gap1 && (seq1[i]==6));
    gap2 = (gap2 && (seq2[i]==6));
  }
  // Find end of internal part of alignment
  j=len-1;
  gap1 = (seq1[j]==6); 
  gap2 = (seq2[j]==6);
  while(gap1 || gap2) {
    j--;
    gap1 = (gap1 && (seq1[j]==6));
    gap2 = (gap2 && (seq2[j]==6));
  }
  // Calculate hamming distance over internal part
  ham=0;
  for(pos=i;pos<=j;pos++) {
    if(seq1[pos] != seq2[pos]) { ham++; }
  }
  return(ham);
}

