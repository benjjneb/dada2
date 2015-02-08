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

#define MAXMAXD 18
#define TARGET_RAW 0
#define ALIGN_SQUAWK 100000
#define TESTING 0
#define VERBOSE 0
#define tVERBOSE 1
#define BUFFER_SIZE 1250
#define KEY_BUFSIZE 2000
#define SEQLEN 900 // Buffer size for DNA sequences read in from uniques files
#define HASHOCC 20
// #define BAND 50 // Size of band in banded alignments. 0 means no banding.
#define KMER_SIZE 6
#define NERRS 12
#define TRUE  1
#define FALSE 0


/* -------------------------------------------
   -------- STRUCTS OBJECTS STRUCTS ----------
   ------------------------------------------- */

typedef std::pair<double, double> Prob;
// typedef int bool;  // [C++]: has builtin bool type

/* Sub:
 A set of substitutions (position and identity) of one sequence
 in an alignment to another sequence.
 Note: positions will be 0-indexed in the alignment */
typedef struct {
  int nsubs;   // number of substitions
  int *pos;    // sequence position of the substitition
  char *nt0;   // nt in reference seq
  char *nt1;   // different nt in aligned seq
  char *key;   // string of all subs: concatenation of "%c%d%c," % nt0,pos,nt1
} Sub;

// Raw: Container for each unique sequence/abundance
typedef struct {
  char *seq;   // the sequence, stored as C-string with A=1,C=2,G=3,T=4
  int *kmer;   // the kmer vector of this sequence
  int reads;   // number of reads of this unique sequence
  int index;   // The index of this Raw in b->raw[index]
} Raw;

// Fam: Child of Bi, contains raws in Bi with same substitution structure
// An "indel family": sequences are the "same" up to indels
typedef struct {
  char *seq;   // representative sequence
  Raw *center; // representative raw (corresponds to seq)
  int reads;   // number of reads in this fam
  Sub *sub;    // struct of substitutions relative to the Bi
  double lambda; // probability of this sequence being produced from the Bi seq
  double p;    // abundance pval relative to the current Bi
  double pS;   // singleton pval relative to the current Bi
  Raw **raw;   // array of pointers to contained raws
  int nraw;    // number of contained raws
  int maxraw;  // number of raws currently allocated for in **raw
} Fam;

// Bi: This is one cluster or partition. Contains raws grouped in fams.
// Bi stores the sub/lambda/e to from the cluster seq/reads to every raw in B.
// Flagged to recalculate various when changes made to its contents.
typedef struct {
  char *seq;   // representative sequence for the cluster
  Raw *center; // representative raw for the cluster (corresponds to seq)
  int nraw;    // number of raws in Bi
  int reads;   // number of reads in this cluster
  Fam **fam;   // Array of pointers to child fams.
  int nfam;    // number of fams in Bi
  int maxfam;  // number of fams currently allocated for in **fam
  bool update_lambda; // set to true when consensus changes
  bool update_fam; // set to true when consensus changes and when raws are shuffled
  double self; // self-production genotype error probability
  StrMap *sm;  // Hash table from sub->key to string 'fam index+1' '\0'
  Sub **sub;   // Array of pointers to subs with all raws.
  double *lambda; // Array of lambdas with all raws.
  double *e;   // Array of expected read numbers with all raws.
} Bi;

// B: holds all the clusters. The full clustering (or partition).
typedef struct {
  int nclust;
  int nraw;
  int reads;
  int maxclust;
  double err[4][4];
  double score[4][4];
  double gap_pen;
  double omegaA;
  bool use_singletons;
  double omegaS;
  double *lams;
  double *cdf;
  size_t nlam;
  Raw **raw;
  Bi **bi;
} B;

// Unique: holds a sequence and associated abundance
typedef struct {
  char *seq;   // A unique sequence, stored in 1-based index form.
  int length;  // The length of the sequence.
  int reads;   // The number of reads.
} Unique;

// Uniques: contains a collection of uniques.
typedef struct {
  Unique *unique; // Pointer to array of Unique objects
  int nseqs;      // Total number of Unique objects
} Uniques;

/* -------------------------------------------
   -------- METHODS METHODS METHODS ----------
   ------------------------------------------- */

// methods implemented in uniques.c
Uniques *uniques_from_vectors(std::vector< std::string > strings, std::vector< int > abundances);
Uniques *uniques_from_file(const char *f);
int uniques_nseqs(Uniques *uniques);
int uniques_reads(Uniques *uniques, int n);
int uniques_length(Uniques *uniques, int n);
void uniques_sequence(Uniques *uniques, int n, char *seq);
void uniques_free(Uniques *uniques);

// methods implemented in cluster.c
B *b_new(Uniques *uniques, double err[4][4], double score[4][4], double gap_pen, double omegaA, bool use_singletons, double omegaS);
void b_free(B *b);
void b_init(B *b);
void b_shuffle(B *b);
void b_lambda_update(B *b, bool use_kmers, double kdist_cutoff, int band_size);
void b_fam_update(B *b);
void b_consensus_update(B *b);
void b_e_update(B *b);
void b_p_update(B *b);
void b_update_err(B *b, double err[4][4]);
void b_get_trans_matrix(B *b, int32_t obs[4][4]);
int b_bud(B *b);
char **b_get_seqs(B *b);
int *b_get_abunds(B *b);

// methods implemented in misc.c
void nt2int(char *oseq, const char *iseq);
void int2nt(char *oseq, const char *iseq);
void ntcpy(char *oseq, const char *iseq);
char *ntstr(const char *iseq);
char *intstr(const char *iseq);
void b_print(B *b);
void align_print(char **al);
void b_dump(B *b, char *fn);
void err_print(double err[4][4]);
void test_fun(int i);

// method implemented in nwalign_endsfree.c
char **nwalign_endsfree(char *s1, char *s2, double s[4][4], double gap_p, int band);
char **raw_align(Raw *raw1, Raw *raw2, double score[4][4], double gap_p, bool use_kmer, double kdist_cutoff, int band);
int *get_kmer(char *seq, int k);
double kmer_dist(int *kv1, int len1, int *kv2, int len2, int k);
Sub *al2subs(char **al);
double compute_lambda(Sub *sub, double self, double t[4][4]);
double get_self(char *seq, double err[4][4]);

// methods implemented in pval.cpp
void getCDF(std::vector<double>& ps, std::vector<double>& cdf, double err[4][4], int nnt[4], int maxD);

#endif
