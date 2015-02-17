#include "dada.h"
#include <Rcpp.h>
#define OMEGA_A 0.01

using namespace Rcpp;
//' @useDynLib dadac
//' @importFrom Rcpp evalCpp

B *run_dada(Uniques *uniques, double score[4][4], double err[4][4], double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS);

//------------------------------------------------------------------
//' Run DADA on the provided unique sequences/abundance pairs. 
//'
//' @param seqs (Required). Character.
//'  A vector containing all unique sequences in the data set.
//'  Only A/C/G/T/N/- allowed. Ungapped sequences recommended.
//' 
//' @param abundances (Required). Numeric.
//'  A vector of the number of reads of each unique seuqences.
//'  NAs not tolerated. Must be same length as the seqs vector.
//'
//' @param err (Required). Numeric matrix (4x4).
//'
//' @param score (Required). Numeric matrix (4x4).
//' The score matrix used during the alignment.
//'
//' @param gap (Required). A \code{numeric(1)} giving the gap penalty for alignment.
//'
//' @return List.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List dada_uniques(std::vector< std::string > seqs,  std::vector< int > abundances,
                        Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix score, Rcpp::NumericVector gap,
                        Rcpp::LogicalVector use_kmers, Rcpp::NumericVector kdist_cutoff,
                        Rcpp::NumericVector band_size,
                        Rcpp::NumericVector omegaA, 
                        Rcpp::LogicalVector use_singletons, Rcpp::NumericVector omegaS) {
  int i, j, f, s, len1, len2, nrow, ncol;
  double tote;
  Fam *fam;
  
  // Load the seqs/abundances into a Uniques struct
  len1 = seqs.size();
  len2 = abundances.size();
  if(len1 != len2) {
    Rcpp::Rcout << "C: Different input lengths:" << len1 << ", " << len2 << "\n";
    return R_NilValue;
  }
  Uniques *uniques = uniques_from_vectors(seqs, abundances);


  // Copy score into a C style array
  nrow = score.nrow();
  ncol = score.ncol();
  if(nrow != 4 || ncol != 4) {
    Rcpp::Rcout << "C: Score matrix malformed:" << nrow << ", " << ncol << "\n";
    return R_NilValue;
  }
  double c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = score(i,j);
    }
  }

  // Copy err into a C style array
  nrow = err.nrow();
  ncol = err.ncol();
  if(nrow != 4 || ncol != 4) {
    Rcpp::Rcout << "C: Error matrix malformed:" << nrow << ", " << ncol << "\n";
    return R_NilValue;
  }
  double c_err[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_err[i][j] = err(i,j);
    }
  }
  
  // Copy gap into a C double
  len1 = gap.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Gap penalty not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_gap = as<double>(gap);
  
  // Copy use kmers into a C++ bool
  len1 = use_kmers.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Use_kmers not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  bool c_use_kmers = as<bool>(use_kmers);

  // Copy kdist_cutoff into a C double
  len1 = kdist_cutoff.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Kdist cutoff not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_kdist_cutoff = as<double>(kdist_cutoff);

  len1 = band_size.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Band_size not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  int c_band_size = as<int>(band_size);

  len1 = omegaA.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: OmegaA not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_omegaA = as<double>(omegaA);

  len1 = use_singletons.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: use_singletons not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  bool c_use_singletons = as<bool>(use_singletons);

  len1 = omegaS.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: OmegaS not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_omegaS = as<double>(omegaS);

  // Run DADA
  B *bb = run_dada(uniques, c_score, c_err, c_gap, c_use_kmers, c_kdist_cutoff, c_band_size, c_omegaA, c_use_singletons, c_omegaS);
  uniques_free(uniques);

  // Extract output from Bi objects
  char **oseqs = (char **) malloc(bb->nclust * sizeof(char *));
  for(i=0;i<bb->nclust;i++) {
    oseqs[i] = (char *) malloc((strlen(bb->bi[i]->seq)+1) * sizeof(char));
    ntcpy(oseqs[i], bb->bi[i]->seq);
  }

  // Convert to R objects and return
  Rcpp::CharacterVector Rseqs;
  Rcpp::NumericVector Rabunds(bb->nclust);
  Rcpp::NumericVector Rbirth_pvals(bb->nclust);
  Rcpp::CharacterVector Rbirth_types;
  Rcpp::NumericVector Rpvals(bb->nclust);

  for(i=0;i<bb->nclust;i++) {
    Rseqs.push_back(std::string(oseqs[i]));
    Rabunds[i] = bb->bi[i]->reads;
    Rbirth_pvals[i] = bb->bi[i]->birth_pval;
    Rbirth_types.push_back(std::string(bb->bi[i]->birth_type));
    
    tote = 0.0;
    for(j=0;j<bb->nclust;j++) {
      if(i != j) {
        tote += bb->bi[j]->e[bb->bi[i]->center->index];
      }
    }
    Rpvals[i] = calc_pA(1+bb->bi[i]->reads, tote);
  }

  // Get error (or substitution) statistics
  int32_t otrans[4][4];
  b_get_trans_matrix(bb, otrans);

  Rcpp::IntegerMatrix Rtrans(4, 4);  // R INTS ARE SIGNED 32 BIT
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      Rtrans(i,j) = otrans[i][j];
    }
  }
  
  int32_t nts_by_pos[SEQLEN] = {0};
  int32_t subs_by_pos[SEQLEN] = {0};
  
  for(i=0;i<bb->nclust;i++) {
    for(j=0;j<strlen(bb->bi[i]->seq);j++) {
      nts_by_pos[j] += bb->bi[i]->reads;
    }
    
    for(f=0;f<bb->bi[i]->nfam;f++) {
      fam = bb->bi[i]->fam[f];
      if(fam->sub) { // not a NULL sub
        for(s=0;s<fam->sub->nsubs;s++) {
          subs_by_pos[fam->sub->pos[s]]+=fam->reads;
        }
      } else {  // Fams should never have NULL subs
        printf("Warning: Output fam C%iF%i had a NULL sub.\n", i, f);
      }
    } // for(f=0;f<bb->bi[i]->nfam;f++)
  }
  
  Rcpp::NumericVector Rnts_by_pos;
  Rcpp::NumericVector Rsubs_by_pos;
  for(i=0;i<SEQLEN;i++) {
    if(nts_by_pos[i]==0) { break; }
    else {
      Rnts_by_pos.push_back(nts_by_pos[i]);
      Rsubs_by_pos.push_back(subs_by_pos[i]);
    }
  }
  
  // Free memory
  for(i=0;i<bb->nclust;i++) {
    free(oseqs[i]);
  }
  free(oseqs);
  b_free(bb);
  
  // Organize return List  
  Rcpp::DataFrame df_clustering = Rcpp::DataFrame::create(_["sequence"] = Rseqs, _["abundance"]  = Rabunds, _["pval"] = Rpvals, _["birth_type"] = Rbirth_types, _["birth_pval"] = Rbirth_pvals);
  Rcpp::DataFrame df_subs = Rcpp::DataFrame::create(_["nts"] = Rnts_by_pos, _["subs"]  = Rsubs_by_pos);
  return Rcpp::List::create(_["clustering"] = df_clustering, _["substitutions"] = df_subs, _["trans"] = Rtrans);
}

B *run_dada(Uniques *uniques, double score[4][4], double err[4][4], double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS) {
  int newi = 0, round = 1;
  double inflation = 1.0;
  size_t index;
  B *bb;
  bb = b_new(uniques, err, score, gap_pen, omegaA, use_singletons, omegaS, band_size); // New cluster with all sequences in 1 bi and 1 fam
  b_fam_update(bb);     // Organizes raws into fams, makes fam consensus/lambda
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
  newi = b_bud(bb);
  
  while(newi) {
    if(tVERBOSE) printf("C: ----------- Round %i -----------\n", round++);
    b_consensus_update(bb);
    b_lambda_update(bb, use_kmers, kdist_cutoff);

    // Temporarily inflate the E's for the new cluster based on the expected number of reads from its center
    if((int) (bb->bi[newi]->center->reads/bb->bi[newi]->self) > bb->bi[newi]->reads) {
      inflation = (bb->bi[newi]->center->reads/bb->bi[newi]->self)/bb->bi[newi]->reads;
      for(index=0;index<bb->nraw;index++) {
        bb->bi[newi]->e[index] = bb->bi[newi]->e[index] * inflation;
      }
    }

    b_shuffle(bb);
    
    b_consensus_update(bb);
    b_lambda_update(bb, use_kmers, kdist_cutoff);
    b_fam_update(bb); // must have lambda_update before fam_update

    b_p_update(bb);
    newi = b_bud(bb);
  }
  return bb;
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
//' The score matrix used during the alignment.
//'
//' @param gap (Required). A \code{numeric(1)} giving the gap penalty for alignment.
//'
//' @param max_aligns (Required). A \code{numeric(1)} giving the (maximum) number of
//' pairwise alignments to do.
//'
//' @return DataFrame.
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame evaluate_kmers(std::vector< std::string > seqs, Rcpp::NumericMatrix score, Rcpp::NumericVector gap, int band, size_t max_aligns) {
  int i, j, n_iters, stride, minlen, nseqs, len1 = 0, len2 = 0;
  char *seq1, *seq2;
  double c_gap = as<double>(gap);
  double c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = score(i,j);
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
  int *kv1;
  int *kv2;

  for(i=0;i<nseqs;i=i+stride) {
    seq1 = intstr(seqs[i].c_str());
    len1 = strlen(seq1);
    kv1 = get_kmer(seq1, KMER_SIZE);
    for(j=i+1;j<nseqs;j=j+stride) {
      seq2 = intstr(seqs[j].c_str());
      len2 = strlen(seq2);
      kv2 = get_kmer(seq2, KMER_SIZE);

      minlen = (len1 < len2 ? len1 : len2);

      sub = al2subs(nwalign_endsfree(seq1, seq2, c_score, c_gap, band));
      adist[npairs] = ((double) sub->nsubs)/((double) minlen);
      
      kdist[npairs] = kmer_dist(kv1, len1, kv2, len2, KMER_SIZE);
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
//' The score matrix used during the alignment.
//'
//' @param gap (Required). A \code{numeric(1)} giving the gap penalty for alignment.
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
//  double c_gap = as<double>(gap);
  double c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = score(i,j);
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

      sub = al2subs(nwalign_endsfree(seq1, seq2, c_score, gap, 0));
      sub_band = al2subs(nwalign_endsfree(seq1, seq2, c_score, gap, band_size));
      
      if(strcmp(sub->key, sub_band->key) != 0) { // different strings
        if(tVERBOSE) { printf("\n%s\n%s\n", sub->key, sub_band->key); }
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


