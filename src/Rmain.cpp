#include "dada.h"
#include <Rcpp.h>

using namespace Rcpp;
//' @useDynLib dada2
//' @importFrom Rcpp evalCpp

B *run_dada(Raw **raws, int nraw, Rcpp::NumericMatrix errMat, int score[4][4], int gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust, double min_fold, int min_hamming, bool use_quals, int qmax, bool final_consensus, bool vectorized_alignment, bool verbose);

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
                        bool use_singletons, double omegaS,
                        int max_clust,
                        double min_fold, int min_hamming,
                        bool use_quals,
                        int qmax,
                        bool final_consensus,
                        bool vectorized_alignment,
                        bool verbose) {

  unsigned int i, j, index, pos, seqlen, nraw;
  
  /********** INPUT VALIDATION *********/
  // Check lengths of seqs and abundances vectors
  if(seqs.size() != abundances.size()) {
    Rcpp::stop("Sequence and abundance vectors had different lengths.");
  }
  nraw = seqs.size();
  if(nraw == 0) {
    Rcpp::stop("Zero input sequences.");
  }
  // Check sequence lengths
  seqlen = seqs[0].length();
  for(index=1;index<nraw;index++) {
    if(seqs[index].length() != seqlen) {
      Rcpp::stop("All input sequences must be the same length.");
    }
  }
  if(seqlen >= SEQLEN) {   // Need one extra byte for the string termination character
    Rcpp::stop("Input sequences exceed the maximum allowed string length.");
  }
  // Check for presence of quality scores and their lengths
  bool has_quals = false;
  if(quals.nrow() > 0) { // Each sequence is a COLUMN, each row is a POSITION
    has_quals = true;
    if(quals.nrow() != seqlen) {
      Rcpp::stop("Sequence lengths and quality lengths must be the same.");
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
  if(err.ncol() == 0) {
    Rcpp::stop("Error matrix must have >0 columns.");
  }

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
      for(pos=0;pos<seqlen;pos++) {
        qual[pos] = quals(pos, index);
      }
      raws[index] = raw_new(seq, qual, abundances[index]);
    } else {
      raws[index] = raw_new(seq, NULL, abundances[index]);
    }
    raws[index]->index = index;
  }

  /********** RUN DADA *********/
  B *bb = run_dada(raws, nraw, err, c_score, gap, use_kmers, kdist_cutoff, band_size, omegaA, use_singletons, omegaS, max_clust, min_fold, min_hamming, use_quals, qmax, final_consensus, vectorized_alignment, verbose);

  /********** MAKE OUTPUT *********/
  Rcpp::DataFrame df_clustering = b_make_clustering_df(bb, has_quals);
  Rcpp::IntegerMatrix mat_trans = b_make_transition_by_quality_matrix(bb, has_quals, qmax);
  Rcpp::NumericMatrix mat_quals = b_make_cluster_quality_matrix(bb, has_quals, seqlen);
  Rcpp::DataFrame df_birth_subs = b_make_birth_subs_df(bb, has_quals);
  Rcpp::DataFrame df_expected = b_make_positional_substitution_df(bb, seqlen, err);
  
  // Make map from uniques to cluster
  Rcpp::IntegerVector Rmap(nraw);
  unsigned int f, r;
  for(i=0;i<bb->nclust;i++) {
    for(f=0;f<bb->bi[i]->nfam;f++) {
      for(r=0;r<bb->bi[i]->fam[f]->nraw;r++) {
        Rmap(bb->bi[i]->fam[f]->raw[r]->index) = i+1; // +1 for R 1-indexing
      }
    }
  }

  // Free memory
  b_free(bb);
  for(index=0;index<nraw;index++) {
    raw_free(raws[index]);
  }
  free(raws);
  
  // Organize return List  
  return Rcpp::List::create(_["clustering"] = df_clustering, _["subpos"] = df_birth_subs, _["subqual"] = mat_trans, _["clusterquals"] = mat_quals, _["map"] = Rmap, _["exp"] = df_expected);
}

B *run_dada(Raw **raws, int nraw, Rcpp::NumericMatrix errMat, int score[4][4], int gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust, double min_fold, int min_hamming, bool use_quals, int qmax, bool final_consensus, bool vectorized_alignment, bool verbose) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;
  
  B *bb;
  bb = b_new(raws, nraw, score, gap_pen, omegaA, use_singletons, omegaS, band_size, vectorized_alignment, use_quals); // New cluster with all sequences in 1 bi and 1 fam
  b_lambda_update(bb, FALSE, 1.0, errMat, verbose); // Everyone gets aligned within the initial cluster, no KMER screen
  b_fam_update(bb, verbose);     // Organizes raws into fams, makes fam consensus/lambda
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
  
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, verbose)) ) {
    if(verbose) Rprintf("----------- New Cluster C%i -----------\n", newi);
    b_lambda_update(bb, use_kmers, kdist_cutoff, errMat, verbose);
    // Keep shuffling and updating until no more shuffles
    nshuffle = 0;
    do {
      shuffled = b_shuffle(bb);
      b_e_update(bb);
      if(verbose) { Rprintf("S"); }
    } while(shuffled && ++nshuffle < MAX_SHUFFLE);
    if(verbose && nshuffle >= MAX_SHUFFLE) { Rprintf("\nWarning: Reached maximum (%i) shuffles.\n", MAX_SHUFFLE); }
    
    b_fam_update(bb, verbose); // If centers can move, must have lambda_update before fam_update
    b_p_update(bb);
    bi_free_absent_subs(bb->bi[newi], bb->nraw); // Free subs in this cluster if they didn't join
  } // while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, verbose)) )
  
  if(final_consensus) { b_make_consensus(bb); }
  if(verbose) Rprintf("\nALIGN: %i aligns, %i shrouded (%i raw).\n", bb->nalign, bb->nshroud, bb->nraw);
  
  return bb;
}

