#include "dada.h"

#define TEST_FILE "group_0-27.uniques"
//#define TEST_FILE "A0014B0276_L_noN.uniques"
//#define TEST_FILE "5seq_1000read.unique"
//#define TEST_FILE "test.uniques2"
#define OMEGA_A 0.01

int main() {
  int i, j, index, len, newi;
  double var, lambda;
  Uniques *uniques = uniques_from_file(TEST_FILE);
  char *seq = (char *) malloc(BUFFER_SIZE);
  char **align;
  char **align2;
  int len1, len2;
  double SCORE[4][4] = {{5, -4, -4, -4}, {-4, 5, -4, -4}, {-4, -4, 5, -4}, {-4, -4, -4, 5}};
  double ERR[4][4] = {{0.991, 0.003, 0.003, 0.003}, {0.003, 0.991, 0.003, 0.003}, {0.003, 0.003, 0.991, 0.003}, {0.003, 0.003, 0.003, 0.991}};
  Raw *pop;
  char *seq1 = (char *) malloc(BUFFER_SIZE);
  char *seq2 = (char *) malloc(BUFFER_SIZE);
  Sub *sub;
  
  // Calculate pval from poisson cdf. Check to see if underflow case gets zero correctly.
  double mu = DBL_MIN;  // TEMPORARY SOLUTION?
  double norm = mu;
  double pval = 1 - gsl_cdf_poisson_P(2-1, mu);  // THIS MUST EQUAL ZERO WHEN MU=DBL_MIN ... DOES IT?
  printf("pval from the DBL_MIN: %.4e = %.4e/%.4e\n", pval/norm, pval, norm);
  pval = pval/norm;

  uniques_sequence(uniques, 0, seq1);
  uniques_sequence(uniques, 4, seq2);
  printf("Seq lengths: %i ... %i\n", (int) strlen(seq1), (int) strlen(seq2));
  printf("%s\n", ntstr(seq1));
  printf("%s\n", ntstr(seq2));
  
  align = b_align(seq1, seq2, SCORE, GAPPEN, FALSE);
  align_print(align);
  sub = al2subs(align);
  lambda = compute_lambda(sub, 1., ERR);
  printf("Lambda (on %i subs): %.4e\n\n", sub->nsubs, lambda);
  
  align2 = b_align(seq1, seq2, SCORE, GAPPEN, FALSE);
  align_print(align2);
  sub = al2subs(align2);
  lambda = compute_lambda(sub, 1., ERR);
  printf("Lambda (on %i subs): %.4e\n\n", sub->nsubs, lambda);
  printf("strcmp(align1[0], align2[0]) = %i, strcmp(align1[1], align2[1]) = %i\n", strcmp(align[0], align2[0]), strcmp(align[1], align2[1]));
  
  
  B *bb;
  double err[4][4];
  for(i=0;i<4;i++) { for(j=0;j<4;j++) { err[i][j] = ERR[i][j]; } }
  double prev_err[4][4];
  double delta_err = 16;
  int round, erri=0;
  
  while(delta_err > DBL_MIN && erri<10) {
    erri++;
    bb = b_new(uniques, err, GAPPEN); // New cluster with all sequences in 1 bi and 1 fam
    b_fam_update(bb);     // Organizes raws into fams, makes fam consensus sequence
    b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
    newi = -1;
    printf("Init Errors:\n"); err_print(ERR);
    round = 1;

    while(newi != 0) { //  && newi<2
      newi = b_bud(bb, OMEGA_A); if(tVERBOSE) printf("Budded\n");
      printf("----------- Round %i (newi=%i) -----------\n", round++, newi);
      b_consensus_update(bb); if(tVERBOSE) printf("Consensused\n");
      b_lambda_update(bb, USE_KMERS); if(tVERBOSE) printf("Lambdad\n");
      b_shuffle(bb);  if(tVERBOSE) printf("Shuffled\n");
      b_consensus_update(bb); if(tVERBOSE) printf("Consensused\n");
      b_fam_update(bb); if(tVERBOSE) printf("Fam Updated\n");
      b_p_update(bb);  if(tVERBOSE) printf("Pvaled\n");
    }
  
//    b_print(bb);
//    printf("Dumping\n");
//    b_dump(bb, "test.clus");
    
    for(i=0;i<4;i++) { for(j=0;j<4;j++) { prev_err[i][j] = err[i][j]; } } // copy err to prev
    b_update_err(bb, err);
    delta_err=0.0;
    for(i=0;i<4;i++) { for(j=0;j<4;j++) { delta_err += fabs(prev_err[i][j] - err[i][j]); } }
    printf("\nNew Errors (delta=%.2e):\n\n", delta_err); err_print(err);
  } // while(delta_err > DBL_MIN && erri<10)
  
  free(seq);
  free(seq1); free(seq2);
  uniques_free(uniques);
  return 0;
}
