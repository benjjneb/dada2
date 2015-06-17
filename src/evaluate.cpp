#include "dada.h"
#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------------------------
//' Exposes ends-free Needleman-Wunsh alignment to R.
//' @export
//' @param s1 desc
//' @param s2 desc
//' @param score desc
//' @param gap_p desc
//' @param band desc
// [[Rcpp::export]]
Rcpp::CharacterVector C_nwalign(std::string s1, std::string s2, Rcpp::NumericMatrix score, int gap_p, int band) {
  int i, j;
  char *seq1, *seq2;
  char **al;

  seq1 = (char *) malloc(s1.size()+1); //E
  seq2 = (char *) malloc(s2.size()+1); //E
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
  al = nwalign_endsfree(seq1, seq2, c_score, gap_p, band);
  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);

//  printf("SCORE:\n");
//  for(i=0;i<4;i++) { printf("{%i, %i, %i, %i}\n",c_score[i][0],c_score[i][1],c_score[i][2],c_score[i][3]); }
//  printf("\nGAP: %i, BAND: %i\n", gap_p, band);
//  printf("ALIGN:\n");
//  printf("%s\n%s\n", al[0], al[1]);

  Rcpp::CharacterVector rval;
  rval.push_back(std::string(al[0]));
  rval.push_back(std::string(al[1]));
  
  free(seq1);
  free(seq2);
  free(al[0]);
  free(al[1]);
  return(rval);
}

//------------------------------------------------------------------
//' Calculates the number of matches/mismatches/internal_indels in an alignment.
//' 
// [[Rcpp::export]]
Rcpp::IntegerVector C_eval_pair(std::string s1, std::string s2) {
  int match, mismatch, indel, start, end;
  bool s1gap, s2gap;
  if(s1.size() != s2.size()) {
    printf("Warning: Aligned strings are not the same length.\n");
    return R_NilValue;
  }

  // Find start of align (end of initial gapping)
  s1gap = s2gap = true;
  start = -1;
  do {
    start++;
    s1gap = s1gap && (s1[start] == '-');
    s2gap = s2gap && (s2[start] == '-');
  } while((s1gap || s2gap) && start<s1.size());

  // Find end of align (start of terminal gapping)
  s1gap = s2gap = true;
  end = s1.size();
  do {
    end--;
    s1gap = s1gap && (s1[end] == '-');
    s2gap = s2gap && (s2[end] == '-');
  } while((s1gap || s2gap) && end>=start);

  // Count the matches, mismatches, indels within the internal part of
  // the alignment.
  match = mismatch = indel = 0;
  for(int i=start;i<=end;i++) {
    if(s1[i]=='-' || s2[i]=='-') {
      indel++;
    } else if(s1[i] == s2[i]) {
      match++;
    } else {
      mismatch++;
    }
  }
  
  Rcpp::IntegerVector rval = Rcpp::IntegerVector::create(_["match"]=match, _["mismatch"]=mismatch, _["indel"]=indel);
  return(rval);
}

//------------------------------------------------------------------
//' Calculates the size of perfect overlap on the left and right side
//' of the child sequence (s2) to the aligned parent.
//' 
// [[Rcpp::export]]
Rcpp::IntegerVector C_get_overlaps(std::string s1, std::string s2, int max_shift) {
  int left, right, start, end, i;
  bool s1gap, s2gap, is_nt1, is_nt2;
  if(s1.size() != s2.size()) {
    printf("Warning: Aligned strings are not the same length.\n");
    return R_NilValue;
  }

  // Find start of align (end of initial gapping)
  s1gap = s2gap = true;
  start = -1;
  do {
    start++;
    s1gap = s1gap && (s1[start] == '-');
    s2gap = s2gap && (s2[start] == '-');
  } while((s1gap || s2gap) && start < max_shift && start < s1.size());

  // Find end of align (start of terminal gapping)
  s1gap = s2gap = true;
  end = s1.size();
  do {
    end--;
    s1gap = s1gap && (s1[end] == '-');
    s2gap = s2gap && (s2[end] == '-');
  } while((s1gap || s2gap) && end > s1.size()-max_shift && end > start);


  // Count the length of perfect alignment from the left and right sides
  left=0;
  for(i=0;i<s1.size();i++) {
    is_nt1 = (s1[i] == 'A' || s1[i] == 'C' || s1[i] == 'G' || s1[i] == 'T');
    is_nt2 = (s2[i] == 'A' || s2[i] == 'C' || s2[i] == 'G' || s2[i] == 'T');
    if(is_nt2) { 
      if(i >= start && !(is_nt1 && is_nt2 && s1[i]==s2[i])) { // mismatch or indel
          break;
      } else { // match or in overhang
        left++;
      }
    }
  }
  
  right=0;
  for(i=s1.size()-1;i>=0;i--) {
    is_nt1 = (s1[i] == 'A' || s1[i] == 'C' || s1[i] == 'G' || s1[i] == 'T');
    is_nt2 = (s2[i] == 'A' || s2[i] == 'C' || s2[i] == 'G' || s2[i] == 'T');
    if(is_nt2) { 
      if(i <= end && !(is_nt1 && is_nt2 && s1[i]==s2[i])) { // mismatch or indel
          break;
      } else { // match or in overhang
        right++;
      }
    }
  }
  
  Rcpp::IntegerVector rval = Rcpp::IntegerVector::create(_["left"]=left, _["right"]=right);
  return(rval);
}

//------------------------------------------------------------------
//' Calculates the consensus of two sequences (first sequence wins mismatches).
//' 
// [[Rcpp::export]]
Rcpp::CharacterVector C_pair_consensus(std::string s1, std::string s2) {
  if(s1.size() != s2.size()) {
    printf("Warning: Aligned strings are not the same length.\n");
    return R_NilValue;
  }
  
  char *oseq = (char *) malloc(s1.size()+1); //E
  for(int i=0;i<s1.size();i++) {
    if(s1[i] == s2[i]) {
      oseq[i] = s1[i];
    } else if(s2[i] == '-') {
      oseq[i] = s1[i];
    } else if(s1[i] == '-') {
      oseq[i] = s2[i];
    } else {
      oseq[i] = s1[i]; // s1 wins mismatches
    }
  }
  oseq[s1.size()] = '\0';

  std::string ostr(oseq);
  free(oseq);
  return(ostr);
}

//------------------------------------------------------------------
//' Generate the kmer-distance and the alignment distance from the
//'   given set of sequences. 
//'
//' @param seqs (Required). Character.
//'  A vector containing all unique sequences in the data set.
//'  Only A/C/G/T allowed.
//' 
//' @param score (Required). Numeric matrix (4x4).
//' The score matrix used during the alignment. Coerced to integer.
//'
//' @param gap (Required). A \code{numeric(1)} giving the gap penalty for alignment. Coerced to integer.
//'
//' @param max_aligns (Required). A \code{numeric(1)} giving the (maximum) number of
//' pairwise alignments to do.
//'
//' @return DataFrame.
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame evaluate_kmers(std::vector< std::string > seqs, int kmer_size, Rcpp::NumericMatrix score, int gap, int band, size_t max_aligns) {
  int i, j, n_iters, stride, minlen, nseqs, len1 = 0, len2 = 0;
  char *seq1, *seq2;

  int c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = (int) score(i,j);
    }
  }
  nseqs = seqs.size();
  
  // Find the kdist/align-dist for max_aligns sequence comparisons
  if(max_aligns < (nseqs * (nseqs-1)/2)) { // More potential comparisons than max
    double foo = 2 * sqrt((double) max_aligns);
    n_iters = (int) foo + 2; // n_iters * (n_iters-1)/2 > max_aligns
    stride = nseqs/n_iters;
  } else {
    max_aligns = (nseqs * (nseqs-1)/2);
    n_iters = nseqs;
    stride = 1;
  }

  size_t npairs = 0;
  Rcpp::NumericVector adist(max_aligns);
  Rcpp::NumericVector kdist(max_aligns);
  Sub *sub;
  uint16_t *kv1;
  uint16_t *kv2;

  for(i=0;i<nseqs;i=i+stride) {
    seq1 = intstr(seqs[i].c_str());
    len1 = strlen(seq1);
    kv1 = get_kmer(seq1, kmer_size);
    for(j=i+1;j<nseqs;j=j+stride) {
      seq2 = intstr(seqs[j].c_str());
      len2 = strlen(seq2);
      kv2 = get_kmer(seq2, kmer_size);

      minlen = (len1 < len2 ? len1 : len2);

      sub = al2subs(nwalign_endsfree(seq1, seq2, c_score, gap, band));
      adist[npairs] = ((double) sub->nsubs)/((double) minlen);
      
      kdist[npairs] = kmer_dist(kv1, len1, kv2, len2, kmer_size);
      npairs++;
      free(kv2);
      free(seq2);
      if(npairs >= max_aligns) { break; }
    }
    free(kv1);
    free(seq1);
    if(npairs >= max_aligns) { break; }
  }
  
  if(npairs != max_aligns) {
    Rcpp::Rcout << "Warning: Failed to reach requested number of alignments.\n";
  }
  return Rcpp::DataFrame::create(_["align"] = adist, _["kmer"] = kdist);
}


//------------------------------------------------------------------
//' Quantify the number of alignments altered by banding at the given BAND_SIZE.
//'
//' @param seqs (Required). Character.
//'  A vector containing all unique sequences in the data set.
//'  Only A/C/G/T allowed.
//' 
//' @param score (Required). Numeric matrix (4x4).
//' The score matrix used during the alignment. Coerced to integer.
//'
//' @param gap (Required). A \code{numeric(1)} giving the gap penalty for alignment. Coerced to integer.
//' 
//' @param band_size (Required). A \code{numeric(1)} giving the band size to consider.
//'
//' @param max_aligns (Required). A \code{numeric(1)} giving the (maximum) number of
//' pairwise alignments to do.
//'
//' @return DataFrame.
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame evaluate_band(std::vector< std::string > seqs, Rcpp::NumericMatrix score, int gap, int band_size, size_t max_aligns) {
  int i, j, n_iters, stride, nseqs, len1 = 0, len2 = 0;
  char *seq1, *seq2;

  int c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = (int) score(i,j);
    }
  }
  nseqs = seqs.size();

  // Find the kdist/align-dist for max_aligns sequence comparisons
  if(max_aligns < (nseqs * (nseqs-1)/2)) { // More potential comparisons than max
    double foo = 2 * sqrt((double) max_aligns);
    n_iters = (int) foo + 2; // n_iters * (n_iters-1)/2 > max_aligns
    stride = nseqs/n_iters;
  } else {
    max_aligns = (nseqs * (nseqs-1)/2);
    n_iters = nseqs;
    stride = 1;
  }

  size_t npairs = 0;
  size_t differ = 0;
  Sub *sub, *sub_band;

  for(i=0;i<nseqs;i=i+stride) {
    seq1 = intstr(seqs[i].c_str());
    len1 = strlen(seq1);
    for(j=i+1;j<nseqs;j=j+stride) {
      seq2 = intstr(seqs[j].c_str());
      len2 = strlen(seq2);

      sub = al2subs(nwalign_endsfree(seq1, seq2, c_score, gap, -1));
      sub_band = al2subs(nwalign_endsfree(seq1, seq2, c_score, gap, band_size));
      
      if(strcmp(sub->key, sub_band->key) != 0) { // different strings
        differ++;
      }
      
      npairs++;
      free(seq2);
      if(npairs >= max_aligns) { break; }
    }
    free(seq1);
    if(npairs >= max_aligns) { break; }
  }
  return Rcpp::DataFrame::create(_["nalign"] = npairs, _["ndiff"] = differ);
}


