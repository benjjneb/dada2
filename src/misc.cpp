#include "dada.h"
// [[Rcpp::interfaces(cpp)]]

void err_print(double err[4][4]) {
  for(int i=0;i<4;i++) {
    if(i==0)  { Rprintf("{"); }
    else      { Rprintf(" "); }
    Rprintf("{");
    for(int j=0;j<4;j++) {
      Rprintf("%.6f", err[i][j]);
      if(j<3) { Rprintf(", "); }
    }
    if(i<3) { Rprintf("},\n"); }
    else    { Rprintf("}}\n"); }
  }
}

void align_print(char **al) {
  char *al0, *al1;
  if(!al) { return; }
  al0 = al[0]; al1 = al[1];
  
  Rprintf("%s\n", ntstr(al0));
  for(int i=0;i<strlen(al0);i++) {
    if(al0[i] == al1[i]) { Rprintf("|"); }
    else { Rprintf(" "); }
  }
  Rprintf("\n%s\n", ntstr(al1));
  
}

/* nt2int takes an input string (iseq) and converts it to numeric indices,
  which are stored in the char output array (oseq).
 Need memory at oseq equal to strlen of iseq, which is not checked currently.
  nt2int(seq,seq) gives previous functionality.
*/

void nt2int(char *oseq, const char *iseq) {
  int i, len = strlen(iseq);

  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 'A':
      *oseq = 1;
      break;
    case 'C':
      *oseq = 2;
      break;
    case 'G':
      *oseq = 3;
      break;
    case 'T':
      *oseq = 4;
      break;
    case 'N':
      *oseq = 5;
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      Rprintf("invalid character in input:%c.\n",*iseq);
      *oseq = '\0';
    }
  }
  *oseq = '\0';
  return;
}

/* Converts seq in index form back to nts. */
void int2nt(char *oseq, const char *iseq) {
  int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 1:
      *oseq = 'A';
      break;
    case 2:
      *oseq = 'C';
      break;
    case 3:
      *oseq = 'G';
      break;
    case 4:
      *oseq = 'T';
      break;
    case 5:
      *oseq = 'N';
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      break;
    }
  }
  *oseq = '\0';
  return;
}

/* Convenience function for diagnostic output. */
void ntcpy(char *oseq, const char *iseq) {
  strcpy(oseq, iseq);
  int2nt(oseq, oseq);
}

/* Convenience function for diagnostic output. */
char *ntstr(const char *iseq) {
  char *oseq = (char *) malloc(strlen(iseq)+1); //E
  if (oseq == NULL)  Rcpp::stop("Memory allocation failed!\n");
  
  ntcpy(oseq, iseq);
  return oseq;
}

/* Convenience function for diagnostic input. */
char *intstr(const char *iseq) {
  char *oseq = (char *) malloc(strlen(iseq)+1); //E
  if (oseq == NULL)  Rcpp::stop("Memory allocation failed!\n");
  
  strcpy(oseq, iseq);
  nt2int(oseq, oseq);
  return oseq;
}
