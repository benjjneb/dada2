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

void b_dump(B *b, char *fn) {
  int i, f, r, index;
  FILE *fp;
  fp = fopen(fn, "w");
  if(fp == NULL) { Rprintf("Null file pointer!\n"); }
  
  for(index=0;index<b->nraw;index++) {
    for(i=0;i<b->nclust;i++) {
      for(f=0;f<b->bi[i]->nfam;f++) {
        for(r=0;r<b->bi[i]->fam[f]->nraw;r++) {
          if(index == b->bi[i]->fam[f]->raw[r]->index) { // The right raw
            fprintf(fp, "%i,%i,%i,%i\n", index, b->raw[index]->reads, i, f);
          }
        }
      }
    }
  }
  fclose(fp);
}

/* [OBJ]_print(OBJ *):
 These functions implement an eye-friendly printout of the cluster structure:
*/
void raw_print(Raw *raw) {
  Rprintf("\t\t\tRaw %i with %i reads\n", raw->index, raw->reads);
}

void fam_print(Fam *fam, int f) {
  Rprintf("\t\tFam %i with %i reads in %i raws.\n", f, fam->reads, fam->nraw);
  for(int r=0; r<fam->nraw; r++) {
    raw_print(fam->raw[r]);
  }
}

void bi_print(Bi *bi, int i) {
  Rprintf("\tCluster %i with %i reads in %i fams...\n", i, bi->reads, bi->nfam);
  for(int f=0; f<bi->nfam; f++) {
    fam_print(bi->fam[f], f);
  }
}

void b_print(B *b) {
  Rprintf("\nClustering of %i raws in %i clusters.........\n", b->nraw, b->nclust);
  for(int i=0; i<b->nclust; i++) {
    bi_print(b->bi[i], i);
  }
  Rprintf("\n");
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
      *oseq = 6;
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
    case 6:
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