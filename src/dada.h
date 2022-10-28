#ifndef _DADA_H_
#define _DADA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_map>
#include <unordered_set>
//#include <gsl/gsl_cdf.h>

// Define a variable that is true of SSE2 supported
// Assuming X64 has SSE2, effectively dropping support for pre-SSE X64 arch
#ifdef __x86_64
#define X64 1
#include "emmintrin.h"
#else
#define X64 0
#endif

#define VERBOSE 0
#define SEQLEN 9999 // Buffer size for DNA sequences read in from uniques files
#define TAIL_APPROX_CUTOFF 1e-7 // Should test to find optimal
#define DBL_PRECISION 1e-15 // precision of doubles
#define KMER_SIZE 5
#define TRUE  1
#define FALSE 0
#define MAX_SHUFFLE 10
#define GAP_GLYPH 9999
#define GRAIN_SIZE 10
#define CACHE_STRIDE 64 // up to 64kb of distance between compared kmer8 vectors


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
  uint8_t *q0;  // quality score in reference seq
  uint8_t *q1;  // quality score in aligned seq
} Sub;

// Raw: Container for each unique sequence/abundance
typedef struct {
  char *seq;   // the sequence, stored as C-string with A=1,C=2,G=3,T=4
  uint8_t *qual; // the rounded average qualities at each position for this unique
  bool prior;  // there are (not) prior reasons to expect this sequence to exist
  uint16_t *kmer;   // the kmer vector of this sequence
  uint8_t *kmer8;   // the kmer vector of this sequence
  uint16_t *kord;   // the kmers in order of this sequence
  unsigned int length;  // the length of the sequence
  unsigned int reads;   // number of reads of this unique sequence
  unsigned int index;   // The index of this Raw in b->raw[index]
  double p;    // abundance pval relative to the current Bi
  double E_minmax;
  Comparison comp;  // the comparison between this Raw and...
  bool lock;   // "Locks" the raw to its current Bi in greedy-mode
  bool correct; // this Raw will be corrected to the center of its partition
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
  bool update_e; // set to true when consensus changes and when raws are shuffled
  bool shuffle; // set to true when e-values are updated
  bool check_locks; // set to true when raws should be checked for locking to that Bi
  double self; // self-production genotype error probability
  unsigned int totraw; // number of total raws in the clustering
  char birth_type[2]; // encoding of how this Bi was created: "I": Initial cluster, "A": Abundance pval, "S": Singleton pval
  unsigned int birth_from; // the cluster from which this new cluster was born
  double birth_pval; // the Bonferonni-corrected pval that led to this cluster being initialized
  double birth_fold; // the multiple of expectations at birth
  double birth_e; // the expected number of reads at birth
  Comparison birth_comp; // the Comparison object at birth
  std::vector<Comparison> comp; // all Comparisons with raws that could potentially join this Bi
} Bi;

// B: holds all the clusters. The full clustering (or partition).
typedef struct {
  unsigned int nclust;
  unsigned int nraw;
  unsigned int reads;
  unsigned int maxclust;
  unsigned int nalign;
  unsigned int nshroud;
  double omegaA;
  double omegaP;
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

// methods implemented in containers.cpp
B *b_new(Raw **raws, unsigned int nraw, double omegaA, double omegaP, bool use_quals);
void b_free(B *b);
void b_init(B *b);
Raw *raw_new(char *seq, double *qual, unsigned int reads, bool prior);
void raw_free(Raw *raw);
Bi *bi_new(unsigned int totraw);
void bi_free(Bi *bi);
unsigned int b_add_bi(B *b, Bi *bi);
Raw *bi_pop_raw(Bi *bi, unsigned int r);
unsigned int bi_add_raw(Bi *bi, Raw *raw);

// methods implemented in cluster.cpp
void b_compare(B *b, unsigned int i, Rcpp::NumericMatrix errMat, int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, bool vectorized_alignment, int SSE, bool gapless, bool greedy, bool verbose);
void b_compare_parallel(B *b, unsigned int i, Rcpp::NumericMatrix errMat, int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, bool vectorized_alignment, int SSE, bool gapless, bool greedy, bool verbose);
bool b_shuffle2(B *b);
// void b_p_update_parallel(B *b);
int b_bud(B *b, double min_fold, int min_hamming, int min_abund, bool verbose);
void bi_census(Bi *bi);
void bi_assign_center(Bi *bi);

// methods implemented in misc.cpp
void nt2int(char *oseq, const char *iseq);
void int2nt(char *oseq, const char *iseq);
void ntcpy(char *oseq, const char *iseq);
char *ntstr(const char *iseq);
char *intstr(const char *iseq);
void align_print(char **al);
void err_print(double err[4][4]);
void test_fun(int i);

// method implemented in nwalign_endsfree.cpp
char **nwalign(const char *s1, size_t len1, const char *s2, size_t len2, int score[4][4], int gap_p, int band);
char **nwalign_endsfree(const char *s1, size_t len1, const char *s2, size_t len2, int score[4][4], int gap_p, int band);
char **nwalign_endsfree_homo(const char *s1, size_t len1, const char *s2, size_t len2, int score[4][4], int gap_p, int gap_homo_p, int band);
char **nwalign_vectorized2(const char *s1, size_t len1, const char *s2, size_t len2, int16_t match, int16_t mismatch, int16_t gap_p, int16_t end_gap_p, int band);
char **nwalign_gapless(const char *s1, size_t len1, const char *s2, size_t len2);
char **raw_align(Raw *raw1, Raw *raw2, int match, int mismatch, int gap_p, int homo_gap_p, bool use_kmers, 
                 double kdist_cutoff, int band, bool vectorized_alignment, int SSE, bool gapless);
Sub *sub_new(Raw *raw0, Raw *raw1, int match, int mismatch, int gap_p, int homo_gap_p, bool use_kmers, double kdist_cutoff, 
             int band, bool vectorized_alignment, int SSE, bool gapless);
Sub *al2subs(char **al);
Sub *sub_copy(Sub *sub);
void sub_free(Sub *sub);

// methods implemented in kmers.cpp
void assign_kmer(uint16_t *kvec, const char *seq, int k);
void assign_kmer_order(uint16_t *kord, char *seq, int k);
void assign_kmer8(uint8_t *kvec8, const char *seq, int k);
double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k);
double kmer_dist_SSEi(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k);
double kmer_dist_SSEi_8(uint8_t *kv1, int len1, uint8_t *kv2, int len2, int k);
double kord_dist(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k);
double kord_dist_SSEi(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k);
///TEST uint16_t kmer_dist2(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k);

// methods implemented in pval.cpp
void b_p_update(B *b, bool greedy, bool detect_singletons);
double calc_pA(int reads, double E_reads, bool prior);
double get_pA(Raw *raw, Bi *bi, bool detect_singletons);
double compute_lambda(Raw *raw, Sub *sub, Rcpp::NumericMatrix errMat, bool use_quals, unsigned int ncol);
double compute_lambda_ts(Raw *raw, Sub *sub, unsigned int ncol, double *err_mat, bool use_quals);
double get_self(char *seq, double err[4][4]);

// methods implemented in error.cpp
Rcpp::DataFrame b_make_clustering_df(B *b, Sub **subs, Sub **birth_subs, bool has_quals);
Rcpp::IntegerMatrix b_make_transition_by_quality_matrix(B *b, Sub **subs, bool has_quals, int ncol);
Rcpp::NumericMatrix b_make_cluster_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int seqlen);
Rcpp::DataFrame b_make_positional_substitution_df(B *b, Sub **subs, unsigned int seqlen, Rcpp::NumericMatrix errMat, bool use_quals);
Rcpp::DataFrame b_make_birth_subs_df(B *b, Sub **birth_subs, bool has_quals);

#endif
