#include "dada.h"
#include <Rcpp.h>
using namespace Rcpp;

int main() {
  int i;
  Uniques *uniques = uniques_from_file("test.uniques");
  char *seq = (char *) malloc(BUFFER_SIZE);
  
  for (i = 0; i < uniques_nseqs(uniques); i++) {
    uniques_sequence(uniques, i, seq);
    int2nt(seq, seq);
    printf("%d\t%s\n", uniques_reads(uniques, i), seq);
  }
  free(seq);

  /* char *s1 = (char *) malloc(BUFFER_SIZE);
  char *s2 = (char *) malloc(BUFFER_SIZE);
  float s[4][4] = { { 5, -4, -4, -4 }, { -4, 5, -4, -4 }, { -4, -4, 5, -4 }, { -4, -4, -4, 5 } };
  int gap_p = -8;
  uniques_sequence(uniques, 0, s1);
  uniques_sequence(uniques, 1, s2);
  char **al = nwalign_endsfree(s1, s2, s, gap_p);
  int2nt(al[0]);
  int2nt(al[1]);
  printf("%s\n%s\n",al[0],al[1]);
  */
  uniques_free(uniques);
  return 0;
}
