#include "dada.h"
#include <Rcpp.h>
///#ifdef _WIN32 //  Windows
/*  //do arm stuff
#endif
#if defined(_WIN32) && !defined(__MINGW32__) //  Windows and not MINGW
#define cpuid(info, x)    __cpuidex(info, x, 0)
#else //  GCC Intrinsics
#include <cpuid.h>
void cpuid(int info[4], int InfoType){
  __cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
}
#endif */

using namespace Rcpp;
//' @useDynLib dada2
//' @importFrom Rcpp evalCpp

B *run_dada(Raw **raws, int nraw, Rcpp::NumericMatrix errMat, 
            int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, 
            double omegaA, double omegaP, bool detect_singletons, 
            int max_clust, double min_fold, int min_hamming, int min_abund, 
            bool use_quals, bool final_consensus, bool vectorized_alignment, bool multithread, bool verbose, 
            int SSE, bool gapless, bool greedy);

//------------------------------------------------------------------
// C interface to run DADA on the provided unique sequences/abundance pairs. 
// 
// [[Rcpp::export]]
Rcpp::List dada_uniques(std::vector< std::string > seqs, std::vector<int> abundances, std::vector<bool> priors,
                        Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix quals,
                        int match, int mismatch, int gap,
                        bool use_kmers, double kdist_cutoff,
                        int band_size,
                        double omegaA, double omegaP, double omegaC, bool detect_singletons,
                        int max_clust,
                        double min_fold, int min_hamming, int min_abund,
                        bool use_quals,
                        bool final_consensus,
                        bool vectorized_alignment,
                        int homo_gap,
                        bool multithread,
                        bool verbose,
                        int SSE,
                        bool gapless,
                        bool greedy) {

  unsigned int i, r, index, pos, nraw, maxlen, minlen;
  Raw *raw;
  
  /********** INPUT VALIDATION *********/
  // Check lengths of seqs and abundances vectors
  nraw = seqs.size();
  if(nraw == 0) { Rcpp::stop("Zero input sequences."); }
  if(nraw != abundances.size()) { Rcpp::stop("Sequence and abundance vectors had different lengths."); }
  if(nraw != priors.size()) { Rcpp::stop("Sequence and priors vectors had different lengths."); }
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
  // Check error matrix
  if(err.nrow() != 16) {
    Rcpp::stop("Error matrix must have 16 rows.");
  }
  // Check for SSE2+ support
  // Code adapted from Stack Overflow (Mysticial) 
  // https://stackoverflow.com/questions/6121792/how-to-check-if-a-cpu-supports-the-sse3-instruction-set
/*  bool HW_SSE = false;
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
  if(!(HW_SSE && HW_SSE2)) { 
    Rprintf("SSE2 not supported.\n");
    SSE = 0;
  } */
  if(!X64) { SSE=0; }
  
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
      raws[index] = raw_new(seq, qual, abundances[index], priors[index]);
    } else {
      raws[index] = raw_new(seq, NULL, abundances[index], priors[index]);
    }
    raws[index]->index = index;
  }
  
  uint8_t *k8;
  uint16_t *k16;
  uint16_t *kord;
  if(use_kmers) {
    /// Add uint8_t kmer index in contiguous memory block
    size_t n_kmer = 1 << (2*KMER_SIZE);
    k8 = (uint8_t *) malloc(nraw * n_kmer * sizeof(uint8_t)); //E
    if (k8 == NULL)  Rcpp::stop("Memory allocation failed.");
    // Construct a raw for each input sequence, store in raws[index]
    for(index=0;index<nraw;index++) {
      raw = raws[index];
      raw->kmer8 = &k8[index*n_kmer];
      assign_kmer8(raw->kmer8, raw->seq, KMER_SIZE);
    }
    
    /// Add uint16_t kmer index in contiguous memory block
    k16 = (uint16_t *) malloc(nraw * n_kmer * sizeof(uint16_t)); //E
    if (k16 == NULL)  Rcpp::stop("Memory allocation failed.");
    // Construct a raw for each input sequence, store in raws[index]
    for(index=0;index<nraw;index++) {
      raw = raws[index];
      raw->kmer = &k16[index*n_kmer];
      assign_kmer(raw->kmer, raw->seq, KMER_SIZE);
    }
    
    /// Add uint16_t ordered kmer record in contiguous memory block
    kord = (uint16_t *) malloc(nraw * maxlen * sizeof(uint16_t)); //E
    if (kord == NULL)  Rcpp::stop("Memory allocation failed.");
    // Construct a raw for each input sequence, store in raws[index]
    for(index=0;index<nraw;index++) {
      raw = raws[index];
      raw->kord = &kord[index*maxlen];
      assign_kmer_order(raw->kord, raw->seq, KMER_SIZE);
    }
  } else {
    for(index=0;index<nraw;index++) {
      raw = raws[index];
      raw->kmer8 = NULL;
      raw->kmer = NULL;
      raw->kord = NULL;
    }
  }
  
  /********** RUN DADA *********/
  B *bb = run_dada(raws, nraw, err, 
                   match, mismatch, gap, homo_gap, use_kmers, kdist_cutoff, band_size, 
                   omegaA, omegaP, detect_singletons,
                   max_clust, min_fold, min_hamming, min_abund, use_quals, final_consensus, 
                   vectorized_alignment, multithread, verbose, SSE, gapless, greedy);

  /********** MAKE OUTPUT *********/
  // Create subs for all the relevant alignments
  Sub **subs = (Sub **) malloc(bb->nraw * sizeof(Sub *)); //E
  Sub **birth_subs = (Sub **) malloc(bb->nclust * sizeof(Sub *)); //E
  if(!subs || !birth_subs) Rcpp::stop("Memory allocation failed.");
  
  // Parallelize over clusters
  struct FinalSubsParallel : public RcppParallel::Worker
  {
    // source data
    B *b;

    // destination Sub arrays
    Sub **subs;
    Sub **birth_subs;
    
    // parameters
    int match, mismatch, gap, homo_gap;
    int band_size;
    bool use_kmers, vectorized_alignment;
    int SSE;
    bool gapless;
    
    // initialize with source and destination
    FinalSubsParallel(B *b, Sub **subs, Sub **birth_subs,
                    int match, int mismatch, int gap, int homo_gap,
                    int band_size, bool use_kmers, bool vectorized_alignment, int SSE, bool gapless) 
      : b(b), subs(subs), birth_subs(birth_subs), match(match), mismatch(mismatch), gap(gap), homo_gap(homo_gap), 
        band_size(band_size), use_kmers(use_kmers), vectorized_alignment(vectorized_alignment), SSE(SSE), gapless(gapless) {}
    
    // Perform sequence comparison
    void operator()(std::size_t begin, std::size_t end) {
      Raw *raw;

      for(std::size_t i=begin;i<end;i++) {
        for(unsigned int r=0;r<b->bi[i]->nraw;r++) {
          raw = b->bi[i]->raw[r];
          subs[raw->index] = sub_new(b->bi[i]->center, raw, match, mismatch, gap, homo_gap, false, 1.0, band_size, vectorized_alignment, SSE, gapless);
        }
        // Make birth sub for that cluster
        if(i==0) { birth_subs[i] = NULL; }
        else {
          birth_subs[i] = sub_new(b->bi[b->bi[i]->birth_comp.i]->center, b->bi[i]->center, match, mismatch, gap, homo_gap, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless);
        }
      } // for(std::size_t i=begin;i<end;i++)
    }
  };
  
  if(multithread) {
    FinalSubsParallel finalSubsParallel(bb, subs, birth_subs, match, mismatch, gap, homo_gap, band_size, use_kmers, vectorized_alignment, SSE, gapless);
    RcppParallel::parallelFor(0, bb->nclust, finalSubsParallel, GRAIN_SIZE);
  } else { // Non-Parallel implementation
    for(i=0;i<bb->nclust;i++) {
      // Make subs for members of that cluster
      for(r=0;r<bb->bi[i]->nraw;r++) {
        raw = bb->bi[i]->raw[r];
        subs[raw->index] = sub_new(bb->bi[i]->center, raw, match, mismatch, gap, homo_gap, false, 1.0, band_size, vectorized_alignment, SSE, gapless);
      }
      // Make birth sub for that cluster
      if(i==0) { birth_subs[i] = NULL; }
      else {
        birth_subs[i] = sub_new(bb->bi[bb->bi[i]->birth_comp.i]->center, bb->bi[i]->center, match, mismatch, gap, homo_gap, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless);
      }
    }
  }

  // Assign final w/in cluster P
  Rcpp::NumericVector Rpraw(nraw); ///! temporary output for evaluation purposes
  for(i=0;i<bb->nclust;i++) {
    for(r=0;r<bb->bi[i]->nraw;r++) {
      raw = bb->bi[i]->raw[r];
      index = raw->index;
      if(bb->bi[i]->center == raw) { // centers always have correct=true
        raw->p = 1.0; // Due to fp precision, possible to get 1's from non-centers in below calc
      } else {
        raw->p = calc_pA(raw->reads, raw->comp.lambda * bb->bi[i]->reads, true);
        if(raw->p < omegaC) { raw->correct = false; } // Don't correct!
      }
      Rpraw[index] = raw->p; ///! temporary output for evaluation purposes
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
      raw = bb->bi[i]->raw[r];
      if(raw->correct) { 
        Rmap(raw->index) = i+1; // +1 for R 1-indexing; 
      } else {
        Rmap(raw->index) = NA_INTEGER;
      }
    }
  }
  
  // Free memory
  b_free(bb);
  for(index=0;index<nraw;index++) {
    raw_free(raws[index]);
  }
  free(raws);
  if(use_kmers) {
    free(k8);
    free(k16);
    free(kord);
  }
  
  // Organize return List  
  return Rcpp::List::create(_["clustering"] = df_clustering, _["birth_subs"] = df_birth_subs, _["subqual"] = mat_trans, _["clusterquals"] = mat_quals, _["map"] = Rmap, _["pval"] = Rpraw);
}

B *run_dada(Raw **raws, int nraw, Rcpp::NumericMatrix errMat, 
            int match, int mismatch, int gap_pen, int homo_gap_pen, bool use_kmers, double kdist_cutoff, int band_size, 
            double omegaA, double omegaP, bool detect_singletons,
            int max_clust, double min_fold, int min_hamming, int min_abund, 
            bool use_quals, bool final_consensus, bool vectorized_alignment, bool multithread, bool verbose, 
            int SSE, bool gapless, bool greedy) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;

  B *bb;
  bb = b_new(raws, nraw, omegaA, omegaP, use_quals); // New cluster with all sequences in 1 bi
  // Everyone gets aligned within the initial cluster, no KMER screen
  if(multithread) { b_compare_parallel(bb, 0, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless, greedy, verbose); }
  else { b_compare(bb, 0, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, 1.0, band_size, vectorized_alignment, SSE, gapless, greedy, verbose); }
//  if(multithread) { b_p_update_parallel(bb); }
  b_p_update(bb, greedy, detect_singletons); // Calculates abundance p-value for each raw in its cluster (consensuses)
  
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, min_abund, verbose)) ) {
    if(verbose) Rprintf("\nNew Cluster C%i:", newi);
    if(multithread) { b_compare_parallel(bb, newi, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless, greedy, verbose); }
    else { b_compare(bb, newi, errMat, match, mismatch, gap_pen, homo_gap_pen, use_kmers, kdist_cutoff, band_size, vectorized_alignment, SSE, gapless, greedy, verbose); }
    // Keep shuffling and updating until no more shuffles
    nshuffle = 0;
    do {
      shuffled = b_shuffle2(bb);
      if(verbose) { Rprintf("S"); }
    } while(shuffled && ++nshuffle < MAX_SHUFFLE);
    if(verbose && nshuffle >= MAX_SHUFFLE) { Rprintf("Warning: Reached maximum (%i) shuffles.\n", MAX_SHUFFLE); }

//    if(multithread) { b_p_update_parallel(bb); }
    b_p_update(bb, greedy, detect_singletons);
    Rcpp::checkUserInterrupt();
  } // while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming, min_abund, verbose)) )
  
  if(verbose) Rprintf("\nALIGN: %i aligns, %i shrouded (%i raw).\n", bb->nalign, bb->nshroud, bb->nraw);
  
  return bb;
}

