#include "dada.h"
#include <Rcpp.h>
using namespace Rcpp;

int get_ham_endsfree(const char *seq1, const char *seq2);
void get_lr(char **al, int &left, int &right, int &left_oo, int &right_oo, bool allow_one_off, int max_shift);

//------------------------------------------------------------------
// Determines whether sq is a perfect bimera of some combination from pars.
// 
// @param sq The query DNA sequence.
// @param pars Potential bimeric "parents".
// @param max_shift A \code{integer(1)} of the maximum alignment shift allowed.
// 
// [[Rcpp::export]]
bool C_is_bimera(std::string sq, std::vector<std::string> pars, bool allow_one_off, int min_one_off_par_dist, int match, int mismatch, int gap_p, int max_shift) {
  int i, left, right, left_oo, right_oo;
  char **al;
  int max_left=0, max_right=0;
  int oo_max_left=0, oo_max_right=0, oo_max_left_oo=0, oo_max_right_oo=0;
  bool rval = false;
  
  for(i=0;i<pars.size() && rval==false;i++) {
    al = nwalign_vectorized2(sq.c_str(), pars[i].c_str(), (int16_t) match, (int16_t) mismatch, (int16_t) gap_p, 0, max_shift);  // Remember, alignments must be freed!
    get_lr(al, left, right, left_oo, right_oo, allow_one_off, max_shift);
    
    if((left+right) >= sq.size()) { // Toss id/pure-shift/internal-indel "parents"
      continue;
    }
    if(left > max_left) { max_left=left; }
    if(right > max_right) { max_right=right; }
    
    // Need to evaluate whether parents are allowed for one-off models
    if(allow_one_off && get_ham_endsfree(al[0], al[1]) >= min_one_off_par_dist) {
      if(left > oo_max_left) { oo_max_left=left; }
      if(right > oo_max_right) { oo_max_right=right; }
      if(left_oo > oo_max_left_oo) { oo_max_left_oo=left_oo; }
      if(right_oo > oo_max_right_oo) { oo_max_right_oo=right_oo; }
    }
    
    // Found bimeric model?
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
  
  return(rval);
}

// [[Rcpp::export]]
Rcpp::DataFrame C_table_bimera(Rcpp::IntegerMatrix mat, std::vector<std::string> seqs, double min_frac, int ignore_n, double min_fold, int min_abund, bool allow_one_off, int min_one_off_par_dist, int match, int mismatch, int gap_p, int max_shift) {
  int i,j,k,nsam,nflag,sqlen,left,right,left_oo,right_oo,max_left,max_right;
  int oo_max_left, oo_max_right, oo_max_left_oo, oo_max_right_oo;
  char **al;
  
  const int *vals = mat.begin(); // What happens if its not integer?
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  // matrix stored in "column-major" order (so 1st col, then 2nd col...)
  // mat(i,j) --> vals[i+j*nrow]
  
  Rcpp::LogicalVector rval(ncol, false);
  Rcpp::IntegerVector flags(ncol, 0);
  Rcpp::IntegerVector sams(ncol, 0);
  Rcpp::IntegerVector temp;
  std::vector<int> lefts(ncol);
  std::vector<int> rights(ncol);
  std::vector<int> lefts_oo(ncol);
  std::vector<int> rights_oo(ncol);
  std::vector<bool> allowed(ncol);
  
  for(j=0;j<ncol;j++) { // Evaluate each sequence
    nsam=0; nflag=0;
    sqlen = seqs[j].size();
    std::fill(lefts.begin(), lefts.end(), -1);
    std::fill(rights.begin(),rights.end(),-1);
    if(allow_one_off) {
      std::fill(lefts_oo.begin(), lefts_oo.end(), -1);
      std::fill(rights_oo.begin(),rights_oo.end(),-1);
      std::fill(allowed.begin(), allowed.end(), false); 
    }
    for(i=0;i<nrow;i++) { // Evaluate in each sample (row)
      if(vals[i+j*nrow]<=0) { continue; }
      nsam++;
      max_left=0; max_right=0;
      oo_max_left=0; oo_max_right=0; oo_max_left_oo=0; oo_max_right_oo=0;
      for(k=0;k<ncol;k++) { // Compare with all possible parents
        if(vals[i+k*nrow]>(min_fold*vals[i+j*nrow]) && vals[i+k*nrow]>=min_abund) {
          if(lefts[k]<0) { // Comparison not yet done to this potential parent
            al = nwalign_vectorized2(seqs[j].c_str(), seqs[k].c_str(), (int16_t) match, (int16_t) mismatch, (int16_t) gap_p, 0, max_shift);  // Remember, alignments must be freed!
            get_lr(al, left, right, left_oo, right_oo, allow_one_off, max_shift);
            if(allow_one_off && get_ham_endsfree(al[0], al[1]) >= min_one_off_par_dist) {
              allowed[k]=true;
            }
            
            if((left+right) < sqlen) {
              lefts[k]=left;
              rights[k]=right;
              if(allow_one_off) {
                lefts_oo[k] = left_oo;
                rights_oo[k] = right_oo;
              }
            } else {  // Ignore id/pure-shift/internal-indel "parents"
              lefts[k]=0;
              rights[k]=0;
              if(allow_one_off) {
                lefts_oo[k] = 0;
                rights_oo[k] = 0;
              }
            }
            free(al[0]);
            free(al[1]);
            free(al);
          }
          // Now compare to best parents yet found
          if(lefts[k] > max_left) { max_left=lefts[k]; }
          if(rights[k] > max_right) { max_right=rights[k]; }
          if(allow_one_off && allowed[k]) {
            if(lefts[k] > oo_max_left) { oo_max_left=lefts[k]; }
            if(rights[k] > oo_max_right) { oo_max_right=rights[k]; }
            if(lefts_oo[k] > oo_max_left_oo) { oo_max_left_oo=lefts_oo[k]; }
            if(rights_oo[k] > oo_max_right_oo) { oo_max_right_oo=rights_oo[k]; }
          }
        } // if(vals[i+k*nrow]>vals[i+j*nrow] && vals[i+k*nrow]>=2)
      } // for(k=0;k<mat.ncol();k++)
      
      // Flag if chimeric model exists
      if((max_right+max_left)>=sqlen) {
        nflag++;
      } else if(allow_one_off) {
        if((oo_max_left+oo_max_right_oo)>=sqlen || (oo_max_left_oo+oo_max_right)>=sqlen) {
          nflag++;
        }
      }
    } // for(i=0;i<mat.nrow();i++)
    
    if(nflag >= nsam || (nflag > 0 && nflag >= (nsam-ignore_n)*min_frac)) { rval[j]=true; }
    flags[j] = nflag;
    sams[j] = nsam;
  } // for(j=0;j<mat.ncol();j++)
  return(Rcpp::DataFrame::create(_["bim"]=rval,_["nflag"]=flags,_["nsam"]=sams));
}

// Internal function to get hamming distance between aligned seqs
//  without counting end gaps
int get_ham_endsfree(const char *seq1, const char *seq2) {
  int i, j, pos, ham;
  size_t len=strlen(seq2);
  bool gap1, gap2;
  // Find start of internal part of alignment
  i=0;
  gap1 = (seq1[i]=='-'); 
  gap2 = (seq2[i]=='-');
  while(gap1 || gap2) {
    i++;
    gap1 = (gap1 && (seq1[i]=='-'));
    gap2 = (gap2 && (seq2[i]=='-'));
  }
  // Find end of internal part of alignment
  j=len-1;
  gap1 = (seq1[j]=='-'); 
  gap2 = (seq2[j]=='-');
  while(gap1 || gap2) {
    j--;
    gap1 = (gap1 && (seq1[j]=='-'));
    gap2 = (gap2 && (seq2[j]=='-'));
  }
  // Calculate hamming distance over internal part
  ham=0;
  for(pos=i;pos<=j;pos++) {
    if(seq1[pos] != seq2[pos]) { ham++; }
  }
  return(ham);
}

// Internal function to get left and right overlaps, 
//   with or without allowing one mismatch
void get_lr(char **al, int &left, int &right, int &left_oo, int &right_oo, bool allow_one_off, int max_shift) {
  size_t len = strlen(al[0]);
  int pos=0; left=0;
  while(al[0][pos] == '-' && pos<len) {
    pos++; // Scan in until query starts
  }
  while(al[1][pos] == '-' && pos<max_shift) {
    pos++; left++; // Credit as ends-free coverage until parent starts
  }
  while(pos<len && al[0][pos] == al[1][pos]) {
    pos++; left++; // Credit as covered until a mismatch
  }
  if(allow_one_off) {
    // Step forward, and credit to one-off further matches (and this mismatch if not a gap)
    left_oo = left;
    pos++;
    if(pos<len && al[0][pos] != '-') { left_oo++; }
    while(pos<len && al[0][pos] == al[1][pos]) {
      pos++; left_oo++;
    }
  }
  
  pos=len-1; right=0;
  while(al[0][pos] == '-' && pos >= 0) { 
    pos--;
  }
  while(al[1][pos] == '-' && pos>+(len-max_shift)) {
    pos--; right++;
  }
  while(pos>=0 && al[0][pos] == al[1][pos]) {
    pos--; right++;
  }
  if(allow_one_off) {
    // Step forward, and credit to one-off further matches (and this mismatch if not a gap)
    right_oo = right;
    pos--;
    if(pos>=0 && al[0][pos] != '-') { right_oo++; }
    while(pos>=0 && al[0][pos] == al[1][pos]) {
      pos--; right_oo++;
    }
  }
}
