#include "dada.h"

int main() {
  int i, index, len, newi;
  double var;
  char seq1[BUFFER_SIZE];
  char seq2[BUFFER_SIZE];
  char **align;
  double SCORE[4][4] = {{5, -4, -4, -4}, {-4, 5, -4, -4}, {-4, -4, 5, -4}, {-4, -4, -4, 5}};
  double ERR[4][4] = {{0.991, 0.003, 0.003, 0.003}, {0.003, 0.991, 0.003, 0.003}, {0.003, 0.003, 0.991, 0.003}, {0.003, 0.003, 0.003, 0.991}};
  
//  strcpy(seq1, "GGTCGGCTAGTCTGGACTCTTTTGTACAGAGCGCGCCTCAAGTTTACGGACGTCCGTTGTAAGGGCATCCGTAACAAGGCCTGGTTAGATCCGGCGAAGC");
//  strcpy(seq2, "GGTTGGCTAGTCTGGACTCTTTTGTACAGAGCGCGCCTCAAGATTACGGACGTCCGTTGTAAGGGCATCCGTAACAAGGCCTGGTTAGACCCGGCGAAGC");
  strcpy(seq1, "GGTCGGCT"); nt2int(seq1, seq1);
  strcpy(seq2, "GGTTGGCT"); nt2int(seq2, seq2);
  
  align = nwalign_endsfree(seq1, seq2, SCORE, GAPPEN);
  align_print(align);
  align = nwalign_endsfree(seq1, seq1, SCORE, GAPPEN);
  align_print(align);
  
  return 0;
}
