#ifndef _DADA_H_
#define _DADA_H_

#define IMPLEMENTATION 'R'

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rcpp.h>
//#include <gsl/gsl_cdf.h>
#include <float.h>
#include "strmap.h" // an ANSI C hash table

// gcc dada.h test.c misc.c uniques.c cluster.c nwalign_endsfree.c strmap.c -lgsl

#define VERBOSE 0
#define tVERBOSE 0
#define BUFFER_SIZE 1250
#define GAPPEN -8
#define KEY_BUFSIZE 2000
#define SEQLEN 900 // Buffer size for DNA sequences read in from uniques files
#define HASHOCC 20
#define BAND 50 // Size of band in banded alignments. 0 means no banding.
#define USE_KMERS 1
#define KMER_SIZE 6
#define KDIST_CUTOFF 0.5

/* objects */
// typedef int bool;  // [C++]: has builtin bool type
#define TRUE  1
#define FALSE 0

/* Sub:
 A set of substitutions (position and identity) of one sequence
 in an alignment to another sequence.
 Note: positions will be 0-indexed in the alignment */
typedef struct {
  int nsubs;
  int *pos;
  char *nt0;
  char *nt1;
  char *key;
} Sub;

// Raw: This is where alignments are made/stored to each Bi
typedef struct { // Raw: The sub-oject of Fam for each unique sequence
  char *seq; // the sequence
  int *kmer; // the kmer vector of this sequence
  int reads; // number of reads of this unique sequence
  int index; // The index of this Raw in b->raw[index]
  double *lambda; // genotype error probability relative to each Bi
  double *e; // expected number of reads from each Bi
} Raw;

// Fam: This is where comparison to the current Bi happens
typedef struct { // Fam: The sub-object of Bi for each indel family
  char *seq; /* the indel family consensus sequence */
  Raw *center; /* the raw that sets the consensus sequence */
  int reads; /* number of reads in this indel family */
  double lambda; /* the genotype error probability relative to current Bi */
  double p; /* the abundance pval relative to the current Bi */
  Sub *sub; /* The Sub struct to current Bi */
  int nraw;
  int maxraw;
  Raw **raw; // An array of pointers to raws.
} Fam;

// Bi: This is one cluster or partition. Contains raws grouped in fams.
// Makes/stores cluster consensus sequence in seq.
// Flagged to recalculate various when changes made to its contents.
typedef struct { // Bi: The sub-object of B for each individual cluster
  char *seq; /* consensus sequence for this cluster */
  Raw *center;
  int reads; /* total number of reads in this cluster */
  int basecount[4]; /* number of ACGTs */
  double self; /* self-production genotype error probability */
  double pS;  /* singleton pval */
  double pSf; /* ?? */
  double pSval; /* ?? */
  bool update_lambda; /* set to true when consensus changes */
  bool update_fam; /* set to true when consensus changes and when raws are shuffled */
  int nfam;
  int nraw;
  int maxfam;
  Sub **sub; // Array of pointers to subs with all raws.
  StrMap *sm; // Hash table from sub->key to string 'fam index+1' '\0'
  double *lambda; // Array of lambdas with all raws.
  double *e;  // Array of expected read numbers with all raws.
  Fam **fam; // Array of pointers to child fams.
} Bi;

// This holds all the clusters. The full partition.
typedef struct { // This is the full cluster object, B.
  int nclust;
  int nraw;
  int reads;
  int maxclust;
  double err[4][4];
  double score[4][4];
  double gap_pen;
  Raw **raw;
  Bi **bi;
} B;

// These structures hold unique sequences and their associated abundances.
typedef struct {
  char *seq; /* A unique sequence, stored in 1-based index form. */
  int length; /* The length of the sequence. */
  int reads; /* The number of reads. */
} Unique;

typedef struct { /* Uniques object. */
  Unique *unique; /* Pointer to array of Unique objects. */
  int nseqs; /* Total number of Unique objects */
} Uniques;

/* methods implemented in uniques.c */

Uniques *uniques_from_vectors(std::vector< std::string > strings, std::vector< int > abundances);
Uniques *uniques_from_file(const char *f);
//uniques *uniques_from_fasta(char *f); TO BE WRITTEN
int uniques_nseqs(Uniques *uniques);
int uniques_reads(Uniques *uniques, int n);
int uniques_length(Uniques *uniques, int n);
void uniques_sequence(Uniques *uniques, int n, char *seq);
void uniques_free(Uniques *uniques);

/* methods implemented in cluster.c */
B *b_new(Uniques *uniques, double err[4][4], double score[4][4], double gap_pen);
void b_free(B *b);
void b_init(B *b);
void b_shuffle(B *b);
void b_reads_update(B *b);
void b_lambda_update(B *b, bool use_kmers, double kdist_cutoff);
void b_fam_update(B *b);
void b_consensus_update(B *b);
void b_e_update(B *b);
void b_p_update(B *b);
void b_update_err(B *b, double err[4][4]);
void b_get_trans_matrix(B *b, int32_t obs[4][4]);
int b_bud(B *b, double omegaA);
char **b_get_seqs(B *b);
int *b_get_abunds(B *b);

void bi_census(Bi *bi);
int b_add_bi(B *b, Bi *bi);
Bi *bi_new(int totraw);
Raw *bi_pop_raw(Bi *bi, int f, int r);
void bi_shove_raw(Bi *bi, Raw *raw);

/* methods implemented in error.c
T *t_new(B *b);
void t_free(T *t);
void t_update(T *t, B *b);
*/

/* methods implemented in misc.c misc */
void nt2int(char *oseq, char *iseq);
void int2nt(char *oseq, char *iseq);
void ntcpy(char *oseq, char *iseq);
char *ntstr(char *iseq);
void b_print(B *b);
void align_print(char **al);
void b_dump(B *b, char *fn);
void err_print(double err[4][4]);
void test_fun(int i);

/* nwalign.c */
char **nwalign_endsfree(char *s1, char *s2, double s[4][4], int gap_p, int band);
char **b_align(char *s1, char *s2, double s[4][4], int gap_p, bool use_kmer, double kdist_cutoff); // DEPRECATED
double compute_lambda(Sub *sub, double self, double t[4][4]);
double get_self(char *seq, double err[4][4]);
Sub *al2subs(char **al);
int *get_kmer(char *seq, int k);
char **raw_align(Raw *raw1, Raw *raw2, double score[4][4], int gap_p, bool use_kmer, double kdist_cutoff);

int timesTwo(int x);

#endif
