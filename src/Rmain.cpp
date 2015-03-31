#include "dada.h"
#include <Rcpp.h>

using namespace Rcpp;
//' @useDynLib dadac
//' @importFrom Rcpp evalCpp

B *run_dada(Raw **raws, int nraw, double score[4][4], double err[4][4], double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust, double min_fold, int min_hamming, bool use_quals, Rcpp::Function lamfun);

//------------------------------------------------------------------
//' Run DADA on the provided unique sequences/abundance pairs. 
//'
//' @param seqs (Required). Character.
//'  A vector containing all unique sequences in the data set.
//'  Only A/C/G/T/N/- allowed. Ungapped sequences recommended.
//' 
//' @param abundances (Required). Numeric.
//'  A vector of the number of reads of each unique seuqences.
//'  NAs not tolerated. Must be same length as the seqs vector.
//'
//' @param err (Required). Numeric matrix (4x4).
//'
//' @param score (Required). Numeric matrix (4x4).
//' The score matrix used during the alignment.
//'
//' @param gap (Required). A \code{numeric(1)} giving the gap penalty for alignment.
//'
//' @return List.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List dada_uniques(std::vector< std::string > seqs,  std::vector< int > abundances,
                        Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix quals,
                        Rcpp::NumericMatrix score, Rcpp::NumericVector gap,
                        Rcpp::LogicalVector use_kmers, Rcpp::NumericVector kdist_cutoff,
                        Rcpp::NumericVector band_size,
                        Rcpp::NumericVector omegaA, 
                        Rcpp::LogicalVector use_singletons, Rcpp::NumericVector omegaS,
                        Rcpp::NumericVector maxClust,
                        Rcpp::NumericVector minFold, Rcpp::NumericVector minHamming,
                        Rcpp::LogicalVector useQuals,
                        Rcpp::Function lamfun) {
  int i, j, len1, len2, nrow, ncol;
  int f, r, nzeros, nones;
  size_t index;
  double tote;
  
  len1 = seqs.size();
  len2 = abundances.size();
  if(len1 != len2) {
    Rcpp::Rcout << "C: Different input lengths:" << len1 << ", " << len2 << "\n";
    return R_NilValue;
  }
  int nraw = len1;
  
  // Turn quals matrix into a c-style 1D array c_quals
  // column-by-column storage
  double *c_quals = NULL;
  int qlen = 0;
  bool has_quals = false;
  if(quals.nrow() > 0) {
    has_quals = true;
    qlen = quals.nrow();
    std::vector<double> cpp_quals = Rcpp::as<std::vector<double> >(quals);
    c_quals = &cpp_quals[0];
  }

  // Make the raws (uniques)
  char *seq;
  Raw **raws = (Raw **) malloc(nraw * sizeof(Raw *));
  for (index = 0; index < nraw; index++) {
    seq = (char *) malloc(seqs[index].length() + 1);
    strcpy(seq, seqs[index].c_str());
    nt2int(seq, seq);
    if(has_quals) {
      raws[index] = raw_qual_new(seq, &c_quals[qlen*index], abundances[index]);
    } else {
      raws[index] = raw_new(seq, abundances[index]);
    }
    raws[index]->index = index;
    free(seq);
  }
  
  // Copy score into a C style array
  nrow = score.nrow();
  ncol = score.ncol();
  if(nrow != 4 || ncol != 4) {
    Rcpp::Rcout << "C: Score matrix malformed:" << nrow << ", " << ncol << "\n";
    return R_NilValue;
  }
  double c_score[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_score[i][j] = score(i,j);
    }
  }

  // Copy err into a C style array
  nrow = err.nrow();
  ncol = err.ncol();
  if(nrow != 4 || ncol != 4) {
    Rcpp::Rcout << "C: Error matrix malformed:" << nrow << ", " << ncol << "\n";
    return R_NilValue;
  }
  double c_err[4][4];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      c_err[i][j] = err(i,j);
    }
  }

  // Check rest of params for Length-One-ness and make into C versions
  if(gap.size() != 1) { Rcpp::Rcout << "C: Gap penalty not length 1\n"; return R_NilValue; }
  double c_gap = as<double>(gap);
  
  // Copy use kmers into a C++ bool
  if(use_kmers.size() != 1) { Rcpp::Rcout << "C: Use_kmers not length 1\n"; return R_NilValue; }
  bool c_use_kmers = as<bool>(use_kmers);

  // Copy kdist_cutoff into a C double
  if(kdist_cutoff.size() != 1) { Rcpp::Rcout << "C: Kdist cutoff not length 1\n"; return R_NilValue; }
  double c_kdist_cutoff = as<double>(kdist_cutoff);

  if(band_size.size() != 1) { Rcpp::Rcout << "C: Band_size not length 1\n"; return R_NilValue; }
  int c_band_size = as<int>(band_size);

  if(omegaA.size() != 1) { Rcpp::Rcout << "C: OmegaA not length 1\n"; return R_NilValue; }
  double c_omegaA = as<double>(omegaA);

  if(use_singletons.size() != 1) { Rcpp::Rcout << "C: use_singletons not length 1\n"; return R_NilValue; }
  bool c_use_singletons = as<bool>(use_singletons);

  if(omegaS.size() != 1) { Rcpp::Rcout << "C: OmegaS not length 1\n"; return R_NilValue; }
  double c_omegaS = as<double>(omegaS);

  if(maxClust.size() != 1) { Rcpp::Rcout << "C: MaxClust not length 1\n"; return R_NilValue; }
  int c_max_clust = as<int>(maxClust);

  if(minFold.size() != 1) { Rcpp::Rcout << "C: MinFold not length 1\n"; return R_NilValue; }
  double c_min_fold = as<double>(minFold);

  if(minHamming.size() != 1) { Rcpp::Rcout << "C: MinHamming not length 1\n"; return R_NilValue; }
  int c_min_hamming = as<int>(minHamming);

  if(useQuals.size() != 1) { Rcpp::Rcout << "C: UseQuals not length 1\n"; return R_NilValue; }
  bool c_use_quals = as<bool>(useQuals);


  // Run DADA
  B *bb = run_dada(raws, nraw, c_score, c_err, c_gap, c_use_kmers, c_kdist_cutoff, c_band_size, c_omegaA, c_use_singletons, c_omegaS, c_max_clust, c_min_fold, c_min_hamming, c_use_quals, lamfun);

  // Extract output from Bi objects
  char **oseqs = (char **) malloc(bb->nclust * sizeof(char *));
  for(i=0;i<bb->nclust;i++) {
    oseqs[i] = (char *) malloc((strlen(bb->bi[i]->seq)+1) * sizeof(char));
    ntcpy(oseqs[i], bb->bi[i]->seq);
  }

  // Convert to R objects and return
  Rcpp::CharacterVector Rseqs;
  Rcpp::NumericVector Rabunds(bb->nclust);
  Rcpp::NumericVector Rzeros(bb->nclust);
  Rcpp::NumericVector Rones(bb->nclust); 
  Rcpp::NumericVector Rraws(bb->nclust); 
  Rcpp::NumericVector Rfams(bb->nclust); 
  Rcpp::NumericVector Rbirth_pvals(bb->nclust);
  Rcpp::NumericVector Rbirth_folds(bb->nclust);
  Rcpp::NumericVector Rbirth_hams(bb->nclust);
  Rcpp::NumericVector Rbirth_es(bb->nclust);
  Rcpp::CharacterVector Rbirth_types;
  Rcpp::NumericVector Rbirth_qaves(bb->nclust);
  Rcpp::NumericVector Rpvals(bb->nclust);

  for(i=0;i<bb->nclust;i++) {
    // The basics
    Rseqs.push_back(std::string(oseqs[i]));
    Rabunds[i] = bb->bi[i]->reads;
    
    // Extra info on the clusters
    nzeros = 0; nones = 0;
    for(f=0;f<bb->bi[i]->nfam;f++) {
      if(bb->bi[i]->fam[f]->sub) {
        if(bb->bi[i]->fam[f]->sub->nsubs == 0) { nzeros += bb->bi[i]->fam[f]->reads; }
        if(bb->bi[i]->fam[f]->sub->nsubs == 1) { nones += bb->bi[i]->fam[f]->reads; }
      }
    }
    Rzeros[i] = nzeros;
    Rones[i] = nones;
    Rraws[i] = bb->bi[i]->nraw;
    Rfams[i] = bb->bi[i]->nfam;
    
    // Record information from the cluster's birth
    Rbirth_types.push_back(std::string(bb->bi[i]->birth_type));
    Rbirth_pvals[i] = bb->bi[i]->birth_pval;
    Rbirth_folds[i] = bb->bi[i]->birth_fold;
    if(bb->bi[i]->birth_sub) { Rbirth_hams[i] = bb->bi[i]->birth_sub->nsubs; }
    else { Rbirth_hams[i] = NA_REAL; }
    Rbirth_es[i] = bb->bi[i]->birth_e;

    // Make qave column
    Sub *sub;
    if(has_quals) {
      double qave = 0.0;
      sub = bb->bi[i]->birth_sub;
      if(sub) {
        for(int s=0;s<sub->nsubs;s++) {
          qave += (sub->q1[s]);
        }
        qave = qave/((double)sub->nsubs);
      }
      Rbirth_qaves[i] = qave;
    }
    
    // Calculate post-hoc pval
    tote = 0.0;
    for(j=0;j<bb->nclust;j++) {
      if(i != j) {
        tote += bb->bi[j]->e[bb->bi[i]->center->index];  // A bit concerned about wahts happening with C0 here (NaN in $pval)
      }
    }
    Rpvals[i] = calc_pA(1+bb->bi[i]->reads, tote); // better to have expected number of center sequence here?
  }
  Rcpp::DataFrame df_clustering = Rcpp::DataFrame::create(_["sequence"] = Rseqs, _["abundance"] = Rabunds, _["n0"] = Rzeros, _["n1"] = Rones, _["nunq"] = Rraws, _["nfam"] = Rfams, _["pval"] = Rpvals, _["birth_type"] = Rbirth_types, _["birth_pval"] = Rbirth_pvals, _["birth_fold"] = Rbirth_folds, _["birth_ham"] = Rbirth_hams, _["birth_qave"] = Rbirth_qaves);

  // Get error (or substitution) statistics
  int32_t otrans[4][4];
  b_get_trans_matrix(bb, otrans);

  Rcpp::IntegerMatrix Rtrans(4, 4);  // R INTS ARE SIGNED 32 BIT
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      Rtrans(i,j) = otrans[i][j];
    }
  }
  
  Rcpp::DataFrame df_subpos = b_get_positional_subs(bb);
  Rcpp::DataFrame df_subqual = b_get_quality_subs(bb);
  
  // Get output average qualities
  Rcpp::NumericMatrix Rquals(qlen, bb->nclust);
  double qq;
  int nreads, pos0, pos1;
  Sub *sub;
  Raw *raw;
  if(has_quals) {
    for(i=0;i<bb->nclust;i++) {
      len1 = strlen(bb->bi[i]->seq);
      if(len1 > qlen) {
        if(tVERBOSE) { printf("Warning: C%i sequence too long for output q matrix (%i vs. %i).\n", i, len1, qlen); }
        len1 = qlen;
      }
      for(pos0=0;pos0<len1;pos0++) {
        qq=0.0;
        nreads=0;
        for(f=0;f<bb->bi[i]->nfam;f++) {
          for(r=0;r<bb->bi[i]->fam[f]->nraw;r++) {
            raw = bb->bi[i]->fam[f]->raw[r];
            sub = bb->bi[i]->sub[raw->index];
            if(sub) {
              pos1 = sub->map[pos0];
              if(pos1 == -999) { // Gap
                continue;
              }
              nreads += raw->reads;
              qq += (raw->qual[pos1] * raw->reads);
            } else {
              printf("Warning: No sub for R%i in C%i\n", r, i);
            }
          }
        }
        Rquals(pos0,i) = qq/nreads;
      } // for(pos0=0;pos0<len1;pos0++)
    }
  }
  
  /*
  // Get output average qualities -- Long form data frame
  std::vector< std::string > long_nt0s;
  std::vector< std::string > long_nt1s;
  std::vector< double > long_qaves;
  std::vector< int > long_abunds;
  std::vector< int > long_subpos;
  char nt0, nt1;
  char buf[2];
  char map_nt[5] = {'\0', 'A', 'C', 'G', 'T'};
  for(i=0;i<bb->nclust;i++) {
    for(pos0=0;pos0<strlen(bb->bi[i]->seq);pos0++) {
      for(f=0;f<bb->bi[i]->nfam;f++) {
        for(r=0;r<bb->bi[i]->fam[f]->nraw;r++) {
          raw = bb->bi[i]->fam[f]->raw[r];
          sub = bb->bi[i]->sub[raw->index];
          if(sub) {
            pos1 = sub->map[pos0];
            if(pos1 == -999) { // Gap
              continue;
            }
          } else {
            printf("Warning: No sub for R%i in C%i\n", r, i);
          }
          nt0 = bb->bi[i]->seq[pos0];  // Remember these are 1=A, 2=C, 3=G, 4=T
          nt1 = raw->seq[pos1];
          sprintf(buf, "%c", map_nt[(int) nt0]);
          long_nt0s.push_back(std::string(buf));
          sprintf(buf, "%c", map_nt[(int) nt1]);
          long_nt1s.push_back(std::string(buf));
          if(has_quals) {
            long_qaves.push_back(raw->qual[pos1]);
          } else {
            long_qaves.push_back(0.0);
          }
          long_abunds.push_back(raw->reads);
          long_subpos.push_back(pos1);
        }
      }
    } // for(pos0=0;pos0<len1;pos0++)
  }
  Rcpp::DataFrame df_sublong = Rcpp::DataFrame::create(_["nt0"] = long_nt0s, _["nt1"] = long_nt1s, _["reads"] = long_abunds, _["qave"] = long_qaves, _["pos"] = long_subpos);
  */
  
  // Free memory
  for(i=0;i<bb->nclust;i++) {
    free(oseqs[i]);
  }
  free(oseqs);
  b_free(bb);
  for(index=0;index<nraw;index++) {
    raw_free(raws[index]);
  }
  free(raws);
  
  // Organize return List  
  return Rcpp::List::create(_["clustering"] = df_clustering, _["subpos"] = df_subpos, _["subqual"] = df_subqual, _["trans"] = Rtrans, _["clusterquals"] = Rquals);
}

B *run_dada(Raw **raws, int nraw, double score[4][4], double err[4][4], double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust, double min_fold, int min_hamming, bool use_quals, Rcpp::Function lamfun) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;
  double inflation = 1.0;
  size_t index;
  
  B *bb;
  bb = b_new(raws, nraw, err, score, gap_pen, omegaA, use_singletons, omegaS, band_size, use_quals); // New cluster with all sequences in 1 bi and 1 fam
//  b_lambda_update(bb, FALSE, 1.0); // Everyone gets aligned within the initial cluster, no KMER screen
//  b_e_update(bb);
  b_fam_update(bb);     // Organizes raws into fams, makes fam consensus/lambda
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
    
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb, min_fold, min_hamming)) ) {
    if(tVERBOSE) printf("----------- New Cluster C%i -----------\n", newi);
    b_consensus_update(bb);
    b_lambda_update(bb, use_kmers, kdist_cutoff);
    // SHOULD GET_SELF ALSO USE QUALS?? I KIND OF THINK NOT...
    
    // Temporarily inflate the E's for the new cluster based on the expected number of reads from its center
    if((int) (bb->bi[newi]->center->reads/bb->bi[newi]->self) > bb->bi[newi]->reads) {
      inflation = (bb->bi[newi]->center->reads/bb->bi[newi]->self)/bb->bi[newi]->reads;
      for(index=0;index<bb->nraw;index++) {
        bb->bi[newi]->e[index] = bb->bi[newi]->e[index] * inflation;
        bb->bi[newi]->update_lambda = TRUE;
      }
    }
    
    // Keep shuffling and updating until no more shuffles
    nshuffle = 0;
    do {
      shuffled = b_shuffle(bb);
      b_consensus_update(bb);
      b_lambda_update(bb, use_kmers, kdist_cutoff);
      if(tVERBOSE) { printf("S"); }
    } while(shuffled && ++nshuffle < MAX_SHUFFLE);
    
    if(tVERBOSE && nshuffle >= MAX_SHUFFLE) { printf("\nWarning: Reached maximum (%i) shuffles.\n", MAX_SHUFFLE); }
    
    b_fam_update(bb); // must have lambda_update before fam_update
    b_p_update(bb);
  }
  return bb;
}

