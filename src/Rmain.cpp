#include "dada.h"
#include <Rcpp.h>
#define OMEGA_A 0.01

using namespace Rcpp;
//' @useDynLib dadac
//' @importFrom Rcpp evalCpp

//int cpp_dada( std::vector< std::string > strings,  std::vector< int > abundnaces ) {
/*   int len = strings.size();
  int len2 = abundances.size();
  
  printf("cpp_dada input: %i strings and %i abundances.\n", len, len2)
  
  if(len != len2) {
    printf("Complain.\n");
  }
*/

//' Run DADA on the provided uniques filename.
//' @export
// [[Rcpp::export]]
int dada_from_file( std::string filename ) {
  int newi, round;
//  double SCORE[4][4] = {{5, -4, -4, -4}, {-4, 5, -4, -4}, {-4, -4, 5, -4}, {-4, -4, -4, 5}};
  double ERR[4][4] = {{0.991, 0.003, 0.003, 0.003}, {0.003, 0.991, 0.003, 0.003}, {0.003, 0.003, 0.991, 0.003}, {0.003, 0.003, 0.003, 0.991}};
  
  const char * cfn = filename.c_str();
  Uniques *uniques = uniques_from_file(cfn);
  
  B *bb;
  bb = b_new(uniques, ERR, GAPPEN); // New cluster with all sequences in 1 bi and 1 fam
  b_fam_update(bb);     // Organizes raws into fams, makes fam consensus sequence
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
  newi = -1;
  if(tVERBOSE) { printf("Init Errors:\n"); err_print(ERR); }
  round = 1;

  while(newi != 0) {
    newi = b_bud(bb, OMEGA_A); if(tVERBOSE) printf("Budded\n");
    printf("----------- Round %i (newi=%i) -----------\n", round++, newi);
    b_consensus_update(bb); if(tVERBOSE) printf("Consensused\n");
    b_lambda_update(bb, USE_KMERS); if(tVERBOSE) printf("Lambdad\n");
    b_shuffle(bb);  if(tVERBOSE) printf("Shuffled\n");
    b_consensus_update(bb); if(tVERBOSE) printf("Consensused\n");
    b_fam_update(bb); if(tVERBOSE) printf("Fam Updated\n");
    b_p_update(bb);  if(tVERBOSE) printf("Pvaled\n");
  }

  b_print(bb);
    
  uniques_free(uniques);
  return round;
}
