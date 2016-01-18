#ifndef _DADA_H_
#define _DADA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rcpp.h>
//#include <gsl/gsl_cdf.h>
#include "strmap.h" // an ANSI C hash table

#define TRACKING 0
#define MAXMAXD 18
#define ALIGN_SQUAWK 100000
#define TESTING 0
#define VERBOSE 0
#define SEQLEN 1000 // Buffer size for DNA sequences read in from uniques files
// SEQLEN MAY NOT BE INCREASED BEYOND 1000 WITHOUT REVISITING AL2SUBS
#define MIN_BUCKETS 10
#define BUCKET_SCALE 0.5
#define TAIL_APPROX_CUTOFF 1e-7 // Should test to find optimal
#define DBL_PRECISION 1e-15 // precision of doubles
#define KMER_SIZE 5
#define NERRS 12
#define TRUE  1
#define FALSE 0
#define MAX_SHUFFLE 10
#define QMIN 0
//#define QMAX 40       // Now a dynamic variable
#define QSTEP 1
#define GAP_GLYPH 9999
#define ALL_SHUFFLE_STEP 50


/* -------------------------------------------
   -------- STRUCTS OBJECTS STRUCTS ----------
   ------------------------------------------- */

/* Comparison:
 A brief summary of the comparison between a cluster and a raw */
typedef struct {
  unsigned int i;
  unsigned int index;
  double lambda;
  unsigned int hamming;
} Comparison;

/* Sub:
 A set of substitutions (position and identity) of one sequence
 in an alignment to another sequence.
 Note: positions will be 0-indexed in the alignment */
typedef struct {
  unsigned int nsubs;   // number of substitions
  unsigned int len0;    // The length of the ref seq
  uint16_t *map;    // map of the sequence position in the ref seq to that in the aligned seq
  uint16_t *pos;    // sequence position of the substitition: index in the reference seq
  char *nt0;   // nt in reference seq
  char *nt1;   // different nt in aligned seq
  double *q0;  // quality in reference seq
  double *q1;  // quality in aligned seq
  char *key;   // string of all subs: concatenation of "%c%d%c," % nt0,pos,nt1
} Sub;

// Raw: Container for each unique sequence/abundance
typedef struct {
  char *seq;   // the sequence, stored as C-string with A=1,C=2,G=3,T=4
  float *qual; // the average qualities at each position for this unique
  uint16_t *kmer;   // the kmer vector of this sequence
  unsigned int length;  // the length of the sequence
  unsigned int reads;   // number of reads of this unique sequence
  unsigned int index;   // The index of this Raw in b->raw[index]
  double p;    // abundance pval relative to the current Bi
  double E_minmax;
  Comparison comp;
} Raw;

// Bi: This is one cluster or partition. Contains raws grouped in fams.
// Bi stores the sub/lambda/e to from the cluster seq/reads to every raw in B.
// Flagged to recalculate various when changes made to its contents.
typedef struct {
  char seq[SEQLEN]; // representative sequence for the cluster
  Raw *center; // representative raw for the cluster (corresponds to seq)
  unsigned int nraw;    // number of raws in Bi
  unsigned int reads;   // number of reads in this cluster
  unsigned int i;       // the cluster number in the total clustering
  Raw **raw;   // Array of pointers to child fams.
  unsigned int maxraw;  // number of fams currently allocated for in **fam
  bool update_lambda; // set to true when consensus changes
  bool update_e; // set to true when consensus changes and when raws are shuffled
  bool shuffle; // set to true when e-values are updated
  double self; // self-production genotype error probability
  unsigned int totraw; // number of total raws in the clustering
  char birth_type[2]; // encoding of how this Bi was created: "I": Initial cluster, "A": Abundance pval, "S": Singleton pval
  double birth_pval; // the Bonferonni-corrected pval that led to this cluster being initialized
  double birth_fold; // the multiple of expectations at birth
  double birth_e; // the expected number of reads at birth
  Comparison birth_comp; // the Comparison object at birth
  std::vector<Comparison> comp;
  std::map<unsigned int, unsigned int> comp_index;
} Bi;

// B: holds all the clusters. The full clustering (or partition).
typedef struct {
  unsigned int nclust;
  unsigned int nraw;
  unsigned int reads;
  unsigned int maxclust;
  int band_size;
  unsigned int nalign;
  unsigned int nshroud;
  int score[4][4];
  int gap_pen;
  bool vectorized_alignment;
  double omegaA;
  bool use_quals;
  double *lams;
  double *cdf;
  size_t nlam;
  Raw **raw;
  Bi **bi;
} B;

/* -------------------------------------------
   -------- METHODS METHODS METHODS ----------
   ------------------------------------------- */

// methods implemented in cluster.c
B *b_new(Raw **raws, unsigned int nraw, int score[4][4], int gap_pen, double omegaA, int band_size, bool vectorized_alignment, bool use_quals);
Raw *raw_new(char *seq, double *qual, unsigned int reads);
void raw_free(Raw *raw);
void b_free(B *b);
void b_init(B *b);
bool b_shuffle2(B *b);
void b_compare(B *b, unsigned int i, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat, bool verbose, unsigned int qmax);
void b_consensus_update(B *b);
//void b_e_update(B *b);
void b_p_update(B *b);
int b_bud(B *b, double min_fold, int min_hamming, bool verbose);
char **b_get_seqs(B *b);
int *b_get_abunds(B *b);
//void b_make_consensus(B *b);

// methods implemented in misc.c
void nt2int(char *oseq, const char *iseq);
void int2nt(char *oseq, const char *iseq);
void ntcpy(char *oseq, const char *iseq);
char *ntstr(const char *iseq);
char *intstr(const char *iseq);
void align_print(char **al);
void err_print(double err[4][4]);
void test_fun(int i);

// method implemented in nwalign_endsfree.c
char **nwalign(char *s1, char *s2, int score[4][4], int gap_p, int band);
char **nwalign_endsfree(char *s1, char *s2, int score[4][4], int gap_p, int band);
char **nwalign_endsfree_vectorized(char *s1, char *s2, int16_t match, int16_t mismatch, int16_t gap_p, int band);
char **raw_align(Raw *raw1, Raw *raw2, int score[4][4], int gap_p, bool use_kmer, double kdist_cutoff, int band, bool vectorized_alignment);
uint16_t *get_kmer(char *seq, int k);
double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k);
Sub *al2subs(char **al);
Sub *sub_new(Raw *raw0, Raw *raw1, int score[4][4], int gap_p, bool use_kmers, double kdist_cutoff, int band, bool vectorized_alignment);
Sub *sub_copy(Sub *sub);
void sub_free(Sub *sub);

// methods implemented in pval.cpp
double calc_pA(int reads, double E_reads);
double get_pA(Raw *raw, Bi *bi);
double compute_lambda(Raw *raw, Sub *sub, Rcpp::NumericMatrix errMat, bool use_quals, unsigned int qmax);
double get_self(char *seq, double err[4][4]);

// methods implemented in error.cpp
Rcpp::DataFrame b_make_clustering_df(B *b, Sub **subs, Sub **birth_subs, bool has_quals);
Rcpp::IntegerMatrix b_make_transition_by_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int qmax);
Rcpp::NumericMatrix b_make_cluster_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int seqlen);
Rcpp::DataFrame b_make_positional_substitution_df(B *b, Sub **subs, unsigned int seqlen, Rcpp::NumericMatrix errMat, bool use_quals);
Rcpp::DataFrame b_make_birth_subs_df(B *b, Sub **birth_subs, bool has_quals);

// methods implemented in taxify.cpp
Rcpp::List C_assign_taxonomy(std::vector<std::string> seqs, std::vector<std::string> refs, std::vector<int> ref_to_genus, Rcpp::IntegerMatrix genusmat);

#endif
