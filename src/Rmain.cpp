#include "dada.h"
#include <Rcpp.h>
#ifdef _WIN32 //  Windows
#define cpuid(info, x)    __cpuidex(info, x, 0)
#else //  GCC Intrinsics
#include <cpuid.h>
void cpuid(int info[4], int InfoType){
  __cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
}
#endif

using namespace Rcpp;
//' @useDynLib dada2
//' @importFrom Rcpp evalCpp

B *run_dada(Raw **raws, int nraw, Rcpp::NumericMatrix errMat, int score[4][4], int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, int max_clust, double min_fold, int min_hamming, bool use_quals, bool final_consensus, bool vectorized_alignment, bool multithread, bool verbose, bool testing);

//------------------------------------------------------------------
// C interface to run DADA on the provided unique sequences/abundance pairs. 
// 
// [[Rcpp::export]]
Rcpp::List dada_uniques(std::vector< std::string > seqs, std::vector<int> abundances,
                        Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix quals,
                        Rcpp::NumericMatrix score, int gap,
                        bool use_kmers, double kdist_cutoff,
                        int band_size,
                        double omegaA, 
                        int max_clust,
                        double min_fold, int min_hamming,
                        bool use_quals,
                        bool final_consensus,
                        bool vectorized_alignment,
                        int homo_gap,
                        bool multithread,
                        bool verbose,
                        bool testing) {

  unsigned int i, j, r, index, pos, nraw, maxlen, minlen;
  
  /********** INPUT VALIDATION *********/
  // Check lengths of seqs and abundances vectors
  if(seqs.size() != abundances.size()) {
    Rcpp::stop("Sequence and abundance vectors had different lengths.");
  }
  nraw = seqs.size();
  if(nraw == 0) {
    Rcpp::stop("Zero input sequences.");
  }
  maxlen=0;
  minlen=SEQLEN;
  for(index=0;index<nraw;index++) {
    if(seqs[index].length() > maxlen) { maxlen = seqs[index].length(); }
    if(seqs[index].length() < minlen) { minlen = seqs[index].length(); }
  }
  if(maxlen >= SEQLEN) { Rcpp::stop("Input sequences exceed the maximum allowed string length."); }
  if(minlen <= KMER_SIZE) { Rcpp::stop("Input sequences must all be longer than the kmer-size (%i).", KMER_SIZE); }
  
  // Check for presence of quality scores and their lengths
  bool has_quals = false;
  if(quals.nrow() > 0) { // Each sequence is a COLUMN, each row is a POSITION
    has_quals = true;
    if(quals.nrow() != maxlen) {
      Rcpp::stop("Sequence must have associated qualities for each nucleotide position.");
    }
  }
  // Copy score matrix into a C style array
  if(score.nrow() != 4 || score.ncol() != 4) {
    Rcpp::stop("Score matrix must be 4x4.");
  }
  int c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = (int) score(i,j);
    }
  }
  // Check error matrix
  if(err.nrow() != 16) {
    Rcpp::stop("Error matrix must have 16 rows.");
  }
  // Check for SSE2+ support
  // Code adapted from Stack Overflow (Mysticial) 
  // https://stackoverflow.com/questions/6121792/how-to-check-if-a-cpu-supports-the-sse3-instruction-set
  bool HW_SSE = false;
  bool HW_SSE2 = false;
  bool HW_SSE3 = false;
  int info[4];
  cpuid(info, 0);
  int nIds = info[0];

  //  Detect Features
  if (nIds >= 0x00000001){
    cpuid(info,0x00000001);
    HW_SSE    = (info[3] & ((int)1 << 25)) != 0;
    HW_SSE2   = (info[3] & ((int)1 << 26)) != 0;
    HW_SSE3   = (info[2] & ((int)1 <<  0)) != 0;
  }
  if(!(HW_SSE && HW_SSE2)) { testing = false; }
  
  /********** CONSTRUCT RAWS *********/
  char seq[SEQLEN];
  double qual[SEQLEN];
  Raw **raws = (Raw **) malloc(nraw * sizeof(Raw *)); //E
  if (raws == NULL)  Rcpp::stop("Memory allocation failed.");
  // Construct a raw for each input sequence, store in raws[index]
  for (index = 0; index < nraw; index++) {
    strcpy(seq, seqs[index].c_str());
    nt2int(seq, seq);
    if(has_quals) {
      for(pos=0;pos<seqs[index].length();pos++) {
        qual[pos] = quals(pos, index);
      }
      raws[index] = raw_new(seq, qual, abundances[index]);
    } else {
      raws[index] = raw_new(seq, NULL, abundances[index]);
    }
    raws[index]->index = index;
  }

  /********** RUN DADA *********/
  B *bb = run_dada(raws, nraw, err, c_score, gap, homo_gap, use_kmers, kdist_cutoff, band_size, omegaA, max_clust, min_fold, min_hamming, use_quals, final_consensus, vectorized_alignment, multithread, verbose, testing);

  /********** MAKE OUTPUT *********/
  Raw *raw;
  
  // Create subs for all the relevant alignments
  Sub **subs = (Sub **) malloc(bb->nraw * sizeof(Sub *)); //E
  Sub **birth_subs = (Sub **) malloc(bb->nclust * sizeof(Sub *)); //E
  if(!subs || !birth_subs) Rcpp::stop("Memory allocation failed.");
  for(i=0;i<bb->nclust;i++) {
    // Make subs for members of that cluster
    for(r=0;r<bb->bi[i]->nraw;r++) {
      raw = bb->bi[i]->raw[r];
      subs[raw->index] = sub_new(bb->bi[i]->center, raw, c_score, gap, homo_gap, false, 1.0, band_size, vectorized_alignment, testing);
    }
    // Make birth sub for that cluster
    if(i==0) { birth_subs[i] = NULL; }
    else {
      birth_subs[i] = sub_new(bb->bi[bb->bi[i]->birth_comp.i]->center, bb->bi[i]->center, c_score, gap, homo_gap, false, 1.0, band_size, vectorized_alignment, testing);
    }
  }
  Rcpp::DataFrame df_clustering = b_make_clustering_df(bb, subs, birth_subs, has_quals);
  Rcpp::IntegerMatrix mat_trans = b_make_transition_by_quality_matrix(bb, subs, has_quals, err.ncol());
  Rcpp::NumericMatrix mat_quals = b_make_cluster_quality_matrix(bb, subs, has_quals, maxlen);
  //  Rcpp::DataFrame df_expected = b_make_positional_substitution_df(bb, subs, seqlen, err, use_quals);
  Rcpp::DataFrame df_birth_subs = b_make_birth_subs_df(bb, birth_subs, has_quals);

  // Free the created subs
  for(index=0;index<bb->nraw;index++) {
    sub_free(subs[index]);
  }
  for(i=0;i<bb->nclust;i++) {
    sub_free(birth_subs[i]);
  }
  
  // Make map from uniques to cluster
  Rcpp::IntegerVector Rmap(nraw);
  for(i=0;i<bb->nclust;i++) {
    for(r=0;r<bb->bi[i]->nraw;r++) {
      Rmap(bb->bi[i]->raw[r]->index) = i+1; // +1 for R 1-indexing
    }
  }

  // Free memory
  b_free(bb);
  for(index=0;index<nraw;index++) {
    raw_free(raws[index]);
  }
  free(raws);
  
  // Organize return List  
  return Rcpp::List::create(_["clustering"] = df_clustering, _["birth_subs"] = df_birth_subs, _["subqual"] = mat_trans, _["clusterquals"] = mat_quals, _["map"] = Rmap);
}

B *run_dada(Raw **raws, int nraw, Rcpp::NumericMatrix errMat, int score[4][4], int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, int max_clust, double min_fold, int min_hamming, bool use_quals, bool final_consensus, bool vectorized_alignment, bool multithread, bool verbose, bool testing) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;

  B *bb;
  bb = b_new(raws, nraw, score, gap_pen, homo_gap_pen, omegaA, band_size, vectorized_alignment, use_quals); // New cluster with all sequences in 1 bi
  // Everyone gets aligned within the initial cluster, no KMER screen
  if(multithread) { b_compare_parallel(bb, 0, FALSE, 1.0, errMat, verbose, testing); }
  else { b_compare(bb, 0, FALSE, 1.0, errMat, verbose, testing); }
  b_p_update(bb);       // Calculates abundance p-value for each raw in its cluster (consensuses)
  
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, verbose)) ) {
    if(verbose) Rprintf("----------- New Cluster C%i -----------\n", newi);
    if(multithread) { b_compare_parallel(bb, newi, use_kmers, kdist_cutoff, errMat, verbose, testing); }
    else { b_compare(bb, newi, use_kmers, kdist_cutoff, errMat, verbose, testing); }
    // Keep shuffling and updating until no more shuffles
    nshuffle = 0;
    do {
      shuffled = b_shuffle2(bb);
      if(verbose) { Rprintf("S"); }
    } while(shuffled && ++nshuffle < MAX_SHUFFLE);
    if(verbose && nshuffle >= MAX_SHUFFLE) { Rprintf("Warning: Reached maximum (%i) shuffles.\n", MAX_SHUFFLE); }

    b_p_update(bb);
    Rcpp::checkUserInterrupt();
  } // while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, verbose)) )
  
  if(verbose) Rprintf("\nALIGN: %i aligns, %i shrouded (%i raw).\n", bb->nalign, bb->nshroud, bb->nraw);
  
  return bb;
}

