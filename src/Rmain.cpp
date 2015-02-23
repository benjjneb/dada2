#include "dada.h"
#include <Rcpp.h>

using namespace Rcpp;
//' @useDynLib dadac
//' @importFrom Rcpp evalCpp

B *run_dada(Uniques *uniques, double score[4][4], double err[4][4], double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust);

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
                        Rcpp::NumericMatrix score, Rcpp::NumericVector gap,
                        Rcpp::LogicalVector use_kmers, Rcpp::NumericVector kdist_cutoff,
                        Rcpp::NumericVector band_size,
                        Rcpp::NumericVector omegaA, 
                        Rcpp::LogicalVector use_singletons, Rcpp::NumericVector omegaS,
                        Rcpp::NumericVector maxClust) {
  int i, j, f, s, len1, len2, nrow, ncol;
  double tote;
  Fam *fam;
  
  // Load the seqs/abundances into a Uniques struct
  len1 = seqs.size();
  len2 = abundances.size();
  if(len1 != len2) {
    Rcpp::Rcout << "C: Different input lengths:" << len1 << ", " << len2 << "\n";
    return R_NilValue;
  }
  Uniques *uniques = uniques_from_vectors(seqs, abundances);


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
  
  // Copy gap into a C double
  len1 = gap.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Gap penalty not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_gap = as<double>(gap);
  
  // Copy use kmers into a C++ bool
  len1 = use_kmers.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Use_kmers not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  bool c_use_kmers = as<bool>(use_kmers);

  // Copy kdist_cutoff into a C double
  len1 = kdist_cutoff.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Kdist cutoff not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_kdist_cutoff = as<double>(kdist_cutoff);

  len1 = band_size.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: Band_size not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  int c_band_size = as<int>(band_size);

  len1 = omegaA.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: OmegaA not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_omegaA = as<double>(omegaA);

  len1 = use_singletons.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: use_singletons not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  bool c_use_singletons = as<bool>(use_singletons);

  len1 = omegaS.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: OmegaS not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  double c_omegaS = as<double>(omegaS);

  len1 = maxClust.size();
  if(len1 != 1) {
    Rcpp::Rcout << "C: MaxClust not length 1:" << len1 << "\n";
    return R_NilValue;
  }
  int c_max_clust = as<int>(maxClust);

  // Run DADA
  B *bb = run_dada(uniques, c_score, c_err, c_gap, c_use_kmers, c_kdist_cutoff, c_band_size, c_omegaA, c_use_singletons, c_omegaS, c_max_clust);
  uniques_free(uniques);

  // Extract output from Bi objects
  char **oseqs = (char **) malloc(bb->nclust * sizeof(char *));
  for(i=0;i<bb->nclust;i++) {
    oseqs[i] = (char *) malloc((strlen(bb->bi[i]->seq)+1) * sizeof(char));
    ntcpy(oseqs[i], bb->bi[i]->seq);
  }

  // Convert to R objects and return
  Rcpp::CharacterVector Rseqs;
  Rcpp::NumericVector Rabunds(bb->nclust);
  Rcpp::NumericVector Rbirth_pvals(bb->nclust);
  Rcpp::CharacterVector Rbirth_types;
  Rcpp::NumericVector Rpvals(bb->nclust);

  for(i=0;i<bb->nclust;i++) {
    Rseqs.push_back(std::string(oseqs[i]));
    Rabunds[i] = bb->bi[i]->reads;
    Rbirth_pvals[i] = bb->bi[i]->birth_pval;
    Rbirth_types.push_back(std::string(bb->bi[i]->birth_type));
    
    tote = 0.0;
    for(j=0;j<bb->nclust;j++) {
      if(i != j) {
        tote += bb->bi[j]->e[bb->bi[i]->center->index];
      }
    }
    Rpvals[i] = calc_pA(1+bb->bi[i]->reads, tote);
  }

  // Get error (or substitution) statistics
  int32_t otrans[4][4];
  b_get_trans_matrix(bb, otrans);

  Rcpp::IntegerMatrix Rtrans(4, 4);  // R INTS ARE SIGNED 32 BIT
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      Rtrans(i,j) = otrans[i][j];
    }
  }
  
  // Move this into a function
  int32_t nts_by_pos[1+4][SEQLEN] = {0};
  int32_t subs_by_pos[1+16][SEQLEN] = {0};
  int sub_ind;
  
  for(i=0;i<bb->nclust;i++) {
    for(j=0;j<strlen(bb->bi[i]->seq);j++) {
      nts_by_pos[0][j] += bb->bi[i]->reads;
      nts_by_pos[(int) bb->bi[i]->seq[j]][j] += bb->bi[i]->reads;
    }
    
    for(f=0;f<bb->bi[i]->nfam;f++) {
      fam = bb->bi[i]->fam[f];
      if(fam->sub) { // not a NULL sub
        for(s=0;s<fam->sub->nsubs;s++) {
          subs_by_pos[0][fam->sub->pos[s]] += fam->reads;
          sub_ind = 1 + 4* ((int) fam->sub->nt0[s] - 1) + ((int) fam->sub->nt1[s] - 1);
          subs_by_pos[sub_ind][fam->sub->pos[s]] += fam->reads;
        }
      } else {  // Fams should never have NULL subs
        printf("Warning: Output fam C%iF%i had a NULL sub.\n", i, f);
      }
    } // for(f=0;f<bb->bi[i]->nfam;f++)
  }
  
  Rcpp::NumericVector Rnts_by_pos;
  Rcpp::NumericVector RAs_by_pos;
  Rcpp::NumericVector RCs_by_pos;
  Rcpp::NumericVector RGs_by_pos;
  Rcpp::NumericVector RTs_by_pos;
  
  Rcpp::NumericVector Rsubs_by_pos;
  Rcpp::NumericVector RAAs_by_pos;
  Rcpp::NumericVector RACs_by_pos;
  Rcpp::NumericVector RAGs_by_pos;
  Rcpp::NumericVector RATs_by_pos;
  Rcpp::NumericVector RCAs_by_pos;
  Rcpp::NumericVector RCCs_by_pos;
  Rcpp::NumericVector RCGs_by_pos;
  Rcpp::NumericVector RCTs_by_pos;
  Rcpp::NumericVector RGAs_by_pos;
  Rcpp::NumericVector RGCs_by_pos;
  Rcpp::NumericVector RGGs_by_pos;
  Rcpp::NumericVector RGTs_by_pos;
  Rcpp::NumericVector RTAs_by_pos;
  Rcpp::NumericVector RTCs_by_pos;
  Rcpp::NumericVector RTGs_by_pos;
  Rcpp::NumericVector RTTs_by_pos;
  for(i=0;i<SEQLEN;i++) {
    if(nts_by_pos[i]==0) { break; }
    else {
      Rnts_by_pos.push_back(nts_by_pos[0][i]);
      RAs_by_pos.push_back(nts_by_pos[1][i]);
      RCs_by_pos.push_back(nts_by_pos[2][i]);
      RGs_by_pos.push_back(nts_by_pos[3][i]);
      RTs_by_pos.push_back(nts_by_pos[4][i]);
      
      Rsubs_by_pos.push_back(subs_by_pos[0][i]);
      RAAs_by_pos.push_back(nts_by_pos[1][i] - subs_by_pos[2][i] - subs_by_pos[3][i] - subs_by_pos[4][i]);
      RACs_by_pos.push_back(subs_by_pos[2][i]);
      RAGs_by_pos.push_back(subs_by_pos[3][i]);
      RATs_by_pos.push_back(subs_by_pos[4][i]);

      RCAs_by_pos.push_back(subs_by_pos[5][i]);
      RCCs_by_pos.push_back(nts_by_pos[2][i] - subs_by_pos[5][i] - subs_by_pos[7][i] - subs_by_pos[8][i]);
      RCGs_by_pos.push_back(subs_by_pos[7][i]);
      RCTs_by_pos.push_back(subs_by_pos[8][i]);

      RGAs_by_pos.push_back(subs_by_pos[9][i]);
      RGCs_by_pos.push_back(subs_by_pos[10][i]);
      RGGs_by_pos.push_back(nts_by_pos[3][i] - subs_by_pos[9][i] - subs_by_pos[10][i] - subs_by_pos[12][i]);
      RGTs_by_pos.push_back(subs_by_pos[12][i]);
      
      RTAs_by_pos.push_back(subs_by_pos[13][i]);
      RTCs_by_pos.push_back(subs_by_pos[14][i]);
      RTGs_by_pos.push_back(subs_by_pos[15][i]);
      RTTs_by_pos.push_back(nts_by_pos[4][i] - subs_by_pos[13][i] - subs_by_pos[14][i] - subs_by_pos[15][i]);
    }
  }
  
  // Free memory
  for(i=0;i<bb->nclust;i++) {
    free(oseqs[i]);
  }
  free(oseqs);
  b_free(bb);
  
  // Organize return List  
  Rcpp::DataFrame df_clustering = Rcpp::DataFrame::create(_["sequence"] = Rseqs, _["abundance"]  = Rabunds, _["pval"] = Rpvals, _["birth_type"] = Rbirth_types, _["birth_pval"] = Rbirth_pvals);
  Rcpp::DataFrame df_subs = Rcpp::DataFrame::create(_["nts"] = Rnts_by_pos, _["subs"] = Rsubs_by_pos, 
          _["A"] = RAs_by_pos, _["C"] = RCs_by_pos, _["G"] = RGs_by_pos, _["T"] = RTs_by_pos,
          _["A2C"] = RACs_by_pos, _["A2G"] = RAGs_by_pos, _["A2T"] = RATs_by_pos, 
          _["C2A"] = RCAs_by_pos, _["C2G"] = RCGs_by_pos, _["C2T"] = RCTs_by_pos,
          _["G2A"] = RGAs_by_pos, _["G2C"] = RGCs_by_pos, _["G2T"] = RGTs_by_pos,
          _["T2A"] = RTAs_by_pos, _["T2C"] = RTCs_by_pos, _["T2G"] = RTGs_by_pos);  // Max 20 cols in Rcpp::DataFrame::create

  return Rcpp::List::create(_["clustering"] = df_clustering, _["subpos"] = df_subs, _["trans"] = Rtrans);
}

B *run_dada(Uniques *uniques, double score[4][4], double err[4][4], double gap_pen, bool use_kmers, double kdist_cutoff, int band_size, double omegaA, bool use_singletons, double omegaS, int max_clust) {
  int newi=0, nshuffle = 0;
  bool shuffled = false;
  double inflation = 1.0;
  size_t index;
  
  B *bb;
  bb = b_new(uniques, err, score, gap_pen, omegaA, use_singletons, omegaS, band_size); // New cluster with all sequences in 1 bi and 1 fam
  b_fam_update(bb);     // Organizes raws into fams, makes fam consensus/lambda
  b_p_update(bb);       // Calculates abundance p-value for each fam in its cluster (consensuses)
//  newi = b_bud(bb);
  
  if(max_clust < 1) { max_clust = bb->nraw; }
  
  while( (bb->nclust < max_clust) && (newi = b_bud(bb)) ) {
    if(tVERBOSE) printf("----------- New Cluster C%i -----------\n", newi);
    b_consensus_update(bb);
    b_lambda_update(bb, use_kmers, kdist_cutoff);
    
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
//    newi = b_bud(bb);
  }
  return bb;
}

