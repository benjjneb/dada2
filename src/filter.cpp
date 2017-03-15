// [[Rcpp::plugins(cpp11)]]
#include "dada.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector C_matchRef(std::vector<std::string> seqs, std::string ref,
                               unsigned int word_size, bool non_overlapping) {
  int i,j;
  std::unordered_set<std::string> phash; ///!
  Rcpp::IntegerVector rval(seqs.size());
  
  unsigned int len = ref.size();
  ref.append(ref, 0, word_size);
  
  for(i=0;i<len;i++) {
    phash.insert(ref.substr(i,word_size));
  }
  
  for(i=0;i<seqs.size();i++) {
    len=seqs[i].size();
    if(len<word_size) { continue; }
    
    for(j=0;j<=(len-word_size);j++) {
      if(phash.count(seqs[i].substr(j, word_size))) {
        rval[i]++;
        if(non_overlapping) { j+=word_size; }
      }
    }
  }
  return(rval);
}

// [[Rcpp::export]]
Rcpp::NumericVector C_matrixEE(Rcpp::IntegerMatrix inp) {
  int i,j;
  double ee;
  Rcpp::NumericVector rval(inp.nrow());

  for(i=0;i<inp.nrow();i++) {
    ee=0.0;
    for(j=0;j<inp.ncol();j++) {
      if(inp(i,j) == NA_INTEGER) { break; }
      ee += pow(10.0, -inp(i,j)/10.0);
    }
    rval(i) = ee;
  }
  return(rval);
}
