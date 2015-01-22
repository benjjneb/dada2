#include "dada.h"
#include <Rcpp.h>
#define OMEGA_A 0.01

using namespace Rcpp;
//' @useDynLib dadac
//' @importFrom Rcpp evalCpp

B *run_dada(Uniques *uniques);

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
//' @return DataFrame object with sequence and abundance columns,
//' corresponding to the DADA denoised sample genotypes.
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame dada_uniques( std::vector< std::string > seqs,  std::vector< int > abundances ) {
  int i;
  int len1 = seqs.size();
  int len2 = abundances.size();
  
  if(len1 != len2) {
    Rcpp::Rcout << "Different input lengths:" << len1 << ", " << len2 << "\n";
    return R_NilValue;
  }
//  for(i=0;i<len1;i++) { Rcpp::Rcout << seqs[i] << "\t" << abundances[i] << "\n"; }
  
  Uniques *uniques = uniques_from_vectors(seqs, abundances);
  B *bb = run_dada(uniques);
  uniques_free(uniques);
  
//  b_print(bb);
  
  char **ostrs = b_get_seqs(bb);
  int *oabs = b_get_abunds(bb);
  
  Rcpp::CharacterVector oseqs;
  Rcpp::NumericVector abunds(bb->nclust);
  for(i=0;i<bb->nclust;i++) {
    oseqs.push_back(std::string(ostrs[i]));
    abunds[i] = oabs[i];
  }
  
  return Rcpp::DataFrame::create(Rcpp::Named("sequence") = oseqs, 
                                Rcpp::Named("abundance")  = abunds);
}

//' Run DADA on the provided uniques filename.
//' @export
// [[Rcpp::export]]
int dada_from_file( std::string filename ) {
  B *bb;
  const char * cfn = filename.c_str();
  Uniques *uniques = uniques_from_file(cfn);
  
  bb = run_dada(uniques);
  uniques_free(uniques);
  
  b_print(bb);
  return bb->nclust;
}

B *run_dada(Uniques *uniques) {
  int newi, round;
//  double SCORE[4][4] = {{5, -4, -4, -4}, {-4, 5, -4, -4}, {-4, -4, 5, -4}, {-4, -4, -4, 5}};
  double ERR[4][4] = {{0.991, 0.003, 0.003, 0.003}, {0.003, 0.991, 0.003, 0.003}, {0.003, 0.003, 0.991, 0.003}, {0.003, 0.003, 0.003, 0.991}};  
  B *bb;
  bb = b_new(uniques, ERR, GAPPEN); // New cluster with all sequences in 1 bi and 1 fam
  b_fam_update(bb);     // Organizes raws into fams, makes fam consensus sequence
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
  newi = -1;
  if(tVERBOSE) { printf("Init Errors:\n"); err_print(ERR); }
  round = 1;

  while(newi != 0) {
    newi = b_bud(bb, OMEGA_A); if(tVERBOSE) Rcout << "Budded\n";
    printf("----------- Round %i (newi=%i) -----------\n", round++, newi);
    b_consensus_update(bb); if(tVERBOSE) Rcout << "Consensused\n";
    b_lambda_update(bb, USE_KMERS); if(tVERBOSE) Rcout << "Lambdad\n";
    b_shuffle(bb);  if(tVERBOSE) Rcout << "Shuffled\n";
    b_consensus_update(bb); if(tVERBOSE) Rcout << "Consensused\n";
    b_fam_update(bb); if(tVERBOSE) Rcout << "Fam Updated\n";
    b_p_update(bb);  if(tVERBOSE) Rcout << "Pvaled\n";
  }

  return bb;
}
