#include "dada.h"
#include <Rcpp.h>
#define OMEGA_A 0.01

using namespace Rcpp;
//' @useDynLib dadac
//' @importFrom Rcpp evalCpp

B *run_dada(Uniques *uniques, double score[4][4], double err[4][4], double gap_pen);

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
//' @return DataFrame object with sequence and abundance columns,
//' corresponding to the DADA denoised sample genotypes.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List dada_uniques(std::vector< std::string > seqs,  std::vector< int > abundances, Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix score, Rcpp::NumericVector gap) {
  int i, j, len1, len2, nrow, ncol;
  
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
  double c_gap = *(gap.begin());
  
  // Run DADA
  B *bb = run_dada(uniques, c_score, c_err, c_gap);
  uniques_free(uniques);
  
  // Extract output from B object
  char **ostrs = b_get_seqs(bb);
  int *oabs = b_get_abunds(bb);
  int32_t trans[4][4];
  b_get_trans_matrix(bb, trans);
  
  // Convert to R objects and return
  Rcpp::CharacterVector oseqs;
  Rcpp::NumericVector oabunds(bb->nclust);
  for(i=0;i<bb->nclust;i++) {
    oseqs.push_back(std::string(ostrs[i]));
    oabunds[i] = oabs[i];
  }
  Rcpp::IntegerMatrix otrans(4, 4);  // R INTS ARE SIGNED 32 BIT
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      otrans(i,j) = trans[i][j];
    }
  }
  
  Rcpp::DataFrame df_genotypes = Rcpp::DataFrame::create(_["sequence"] = oseqs, _["abundance"]  = oabunds);
  return Rcpp::List::create(_["genotypes"] = df_genotypes, _["trans"] = otrans);
}

B *run_dada(Uniques *uniques, double score[4][4], double err[4][4], double gap_pen) {
  int newi, round;
//  double SCORE[4][4] = {{5, -4, -4, -4}, {-4, 5, -4, -4}, {-4, -4, 5, -4}, {-4, -4, -4, 5}};
//  double ERR[4][4] = {{0.991, 0.003, 0.003, 0.003}, {0.003, 0.991, 0.003, 0.003}, {0.003, 0.003, 0.991, 0.003}, {0.003, 0.003, 0.003, 0.991}};  
  B *bb;
  bb = b_new(uniques, err, score, gap_pen); // New cluster with all sequences in 1 bi and 1 fam
  b_fam_update(bb);     // Organizes raws into fams, makes fam consensus sequence
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
  newi = -999;
  round = 1;

  while(newi != 0) {
    newi = b_bud(bb, OMEGA_A); if(tVERBOSE) Rcout << "C: Budded\n";
    if(tVERBOSE) printf("C: ----------- Round %i (newi=%i) -----------\n", round++, newi);
    b_consensus_update(bb); if(tVERBOSE) Rcout << "C: Consensused\n";
    b_lambda_update(bb, USE_KMERS); if(tVERBOSE) Rcout << "C: Lambdad\n";
    b_shuffle(bb);  if(tVERBOSE) Rcout << "C: Shuffled\n";
    b_consensus_update(bb); if(tVERBOSE) Rcout << "C: Consensused\n";
    b_fam_update(bb); if(tVERBOSE) Rcout << "C: Fam Updated\n";
    b_p_update(bb);  if(tVERBOSE) Rcout << "C: Pvaled\n";
  }

  return bb;
}
