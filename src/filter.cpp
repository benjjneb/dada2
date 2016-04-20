////// [[Rcpp::plugins(cpp11)]]
#include "dada.h"
#include <Rcpp.h>
////#include <unordered_set> ///!
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector C_matchRef(std::vector<std::string> seqs, std::string ref,
                             unsigned int word_size, bool non_overlapping) {
  int i,j;
  std::set<std::string> phash; ///!
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


