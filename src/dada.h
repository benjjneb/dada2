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
#define TARGET_RAW -1
#define ALIGN_SQUAWK 100000
#define TESTING 0
#define VERBOSE 0
#define tVERBOSE 0
#define SEQLEN 999 // Buffer size for DNA sequences read in from uniques files
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
#define QMAX 40
#define QSTEP 1
#define GAP_GLYPH 9999


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
  int length;  // the length of the sequence
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
  int i;       // the cluster number in the total clustering
  Fam **fam;   // Array of pointers to child fams.
  int nfam;    // number of fams in Bi
  int maxfam;  // number of fams currently allocated for in **fam
  bool update_lambda; // set to true when consensus changes
  bool update_fam; // set to true when consensus changes and when raws are shuffled
  bool update_e; // set to true when consensus changes and when raws are shuffled
  bool shuffle; // set to true when e-values are updated
  double self; // self-production genotype error probability
  StrMap *sm;  // Hash table from sub->key to string 'fam index+1' '\0'
  Sub **sub;   // Array of pointers to subs with all raws.
  double *lambda; // Array of lambdas with all raws.
  double *e;   // Array of expected read numbers with all raws.
  size_t totraw; // number of total raws in the clustering
  char birth_type[2]; // encoding of how this Bi was created: "I": Initial cluster, "A": Abundance pval, "S": Singleton pval
  double birth_pval; // the Bonferonni-corrected pval that led to this cluster being initialized
  double birth_fold; // the multiple of expectations at birth
  double birth_e; // the expected number of reads at birth
  Sub *birth_sub; // the Sub object at birth
} Bi;

// B: holds all the clusters. The full clustering (or partition).
typedef struct {
  int nclust;
  int nraw;
  int reads;
  int maxclust;
  int band_size;
  int nalign;
  int nshroud;
//  double err[4][4];
  int score[4][4];
  int gap_pen;
  double omegaA;
  bool use_singletons;
  double omegaS;
  bool use_quals;
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
B *b_new(Raw **raws, int nraw, int score[4][4], int gap_pen, double omegaA, bool use_singletons, double omegaS, int band_size, bool use_quals);
Raw *raw_new(char *seq, int reads);
Raw *raw_qual_new(char *seq, double *qual, int reads);
void raw_free(Raw *raw);
void b_free(B *b);
void b_init(B *b);
bool b_shuffle(B *b);
void b_lambda_update(B *b, bool use_kmers, double kdist_cutoff, Rcpp::NumericMatrix errMat);
void b_fam_update(B *b);
void b_consensus_update(B *b);
void b_e_update(B *b);
void b_p_update(B *b);
int b_bud(B *b, double min_fold, int min_hamming);
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
char **nwalign_endsfree(char *s1, char *s2, int score[4][4], int gap_p, int band);
char **raw_align(Raw *raw1, Raw *raw2, int score[4][4], int gap_p, bool use_kmer, double kdist_cutoff, int band);
uint16_t *get_kmer(char *seq, int k);
double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k);
Sub *al2subs(char **al);
Sub *sub_new(Raw *raw0, Raw *raw1, int score[4][4], int gap_p, bool use_kmers, double kdist_cutoff, int band);
void sub_free(Sub *sub);

// methods implemented in pval.cpp
void b_make_pS_lookup(B *b);
void getCDF(std::vector<double>& ps, std::vector<double>& cdf, double err[4][4], int nnt[4], int maxD);
double calc_pA(int reads, double E_reads);
double get_pA(Fam *fam, Bi *bi);
double get_pS(Fam *fam, Bi *bi, B *b);
double compute_lambda(Sub *sub, double self, double t[4][4], bool use_quals);
double compute_lambda3(Raw *raw, Sub *sub, Rcpp::NumericMatrix errMat, bool use_quals);
double get_self(char *seq, double err[4][4]);

// methods implemented in error.cpp
void b_get_trans_matrix(B *b, int32_t obs[4][4]);
Rcpp::DataFrame b_get_positional_subs(B *b);
Rcpp::DataFrame b_get_quality_subs(B *b);
Rcpp::IntegerMatrix b_get_quality_subs2(B *b, bool has_quals, int qmin, int qmax);
Rcpp::DataFrame get_sublong(B *b, bool has_quals);

#endif
