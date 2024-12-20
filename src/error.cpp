#include <Rcpp.h>
#include "dada.h"
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// This function constructs the output "clustering" data.frame for the dada(...) function.
// This contains core and diagnostic information on each partition (or cluster, or Bi).
Rcpp::DataFrame b_make_clustering_df(B *b, Sub **subs, Sub **birth_subs, bool has_quals) {
  unsigned int i, j, r, s, cind, max_reads;
  Raw *max_raw, *raw;
  Sub *sub;
  double q_ave;
  size_t index;
  
  // Create output character-vector of representative sequences for each partition (Bi)
  Rcpp::CharacterVector Rseqs;
  char oseq[SEQLEN];
  for(i=0;i<b->nclust;i++) {
    max_reads=0;
    max_raw = NULL;
    for(r=0;r<b->bi[i]->nraw;r++) {
      if(b->bi[i]->raw[r]->reads > max_reads) {
        max_raw = b->bi[i]->raw[r];
        max_reads = max_raw->reads;
      }
    }
    if(!max_raw) {
      Rseqs.push_back(std::string(""));
    } else {
      ntcpy(oseq, max_raw->seq);
      Rseqs.push_back(std::string(oseq));
    }
  }
  
  // Create output vectors for other columns of clustering data.frame
  // Each has one entry for each partition (Bi)
  Rcpp::IntegerVector Rabunds(b->nclust);      // abundances
  Rcpp::IntegerVector Rzeros(b->nclust);       // n0
  Rcpp::IntegerVector Rones(b->nclust);        // n1
  Rcpp::IntegerVector Rraws(b->nclust);        // nraw
  Rcpp::NumericVector Rbirth_pvals(b->nclust); // pvalue at birth
  Rcpp::NumericVector Rbirth_folds(b->nclust); // fold over-abundance at birth
  Rcpp::IntegerVector Rbirth_hams(b->nclust);  // hamming distance at birth
  Rcpp::NumericVector Rbirth_es(b->nclust);    // expected number at birth
///  Rcpp::CharacterVector Rbirth_types;           // DEPRECATED
  Rcpp::IntegerVector Rbirth_froms(b->nclust); // cluster from which this new cluster was born
  Rcpp::NumericVector Rbirth_qaves(b->nclust); // average quality of substitutions that drove birth
  Rcpp::NumericVector Rpvals(b->nclust);       // post-hoc pvalue

  // Assign values to the output vectors
  for(i=0;i<b->nclust;i++) {
    // abundance, nraw, n0, n1
    Rzeros[i] = 0; Rones[i] = 0;
    for(r=0;r<b->bi[i]->nraw;r++) {
      raw = b->bi[i]->raw[r];
      if(raw->correct) { // Use only correct=true raws (can assume true of center, but 2x check)
        Rabunds[i] += raw->reads;
        Rraws[i]++;
        sub = subs[raw->index];
        if(sub) {
          if(sub->nsubs == 0) { Rzeros[i] += b->bi[i]->raw[r]->reads; }
          if(sub->nsubs == 1) { Rones[i] += b->bi[i]->raw[r]->reads; }
        }
      }
    }
    // Record information from the cluster's birth
///    Rbirth_types.push_back(std::string(b->bi[i]->birth_type));
    if(i==0) {  // 0-clust wasn't born normally
      Rbirth_pvals[i] = NA_REAL; 
      Rbirth_froms[i] = NA_INTEGER;
      Rbirth_folds[i] = NA_REAL; 
      Rbirth_hams[i] = NA_INTEGER; 
      Rbirth_es[i] = NA_REAL;
      Rbirth_qaves[i] = NA_REAL;
    } else { 
      Rbirth_froms[i] = b->bi[i]->birth_from + 1; // R 1-based indexing
      Rbirth_pvals[i] = b->bi[i]->birth_pval;
      Rbirth_folds[i] = b->bi[i]->birth_fold;
      Rbirth_hams[i] = b->bi[i]->birth_comp.hamming;
      Rbirth_es[i] = b->bi[i]->birth_e;
      // Calculate average quality of birth substitutions
      if(has_quals) {
        q_ave = 0.0;
        sub = birth_subs[i];
        if(sub && sub->q1) {
          for(s=0;s<sub->nsubs;s++) {
            q_ave += (sub->q1[s]); // Numeric vector for output
          }
          q_ave = q_ave/((double)sub->nsubs);
        }
        Rbirth_qaves[i] = q_ave;
      } else {
        Rbirth_qaves[i] = Rcpp::NumericVector::get_na();
      }
    }
  }
  
  // Calculate post-hoc pval
  // This is not as exhaustive anymore
  std::unordered_map<unsigned int, unsigned int> center_of;
  for(i=0;i<b->nclust;i++) {
    center_of[b->bi[i]->center->index] = i;
  }
  std::vector<double> tot_e(b->nclust);
  for(i=0;i<b->nclust;i++) {
    for(cind=0;cind<b->bi[i]->comp.size();cind++) {
      index = b->bi[i]->comp[cind].index;
      if(center_of.count(index)) { // This is the center of a Bi
        j = center_of[index];
        if(i != j) { // Only add from other Bis
          tot_e[j] += b->bi[i]->comp[cind].lambda * b->bi[i]->reads;
        }
      }
    }
  }
  for(i=0;i<b->nclust;i++) {
    Rpvals[i] = calc_pA(b->bi[i]->center->reads, tot_e[i], true); // prior=true to get non-conditional p-val
  }
  
  return(Rcpp::DataFrame::create(_["sequence"] = Rseqs, _["abundance"] = Rabunds, 
                                 _["n0"] = Rzeros, _["n1"] = Rones, _["nunq"] = Rraws, 
                                 _["pval"] = Rpvals, /// _["birth_type"] = Rbirth_types, 
                                 _["birth_from"] = Rbirth_froms,
                                 _["birth_pval"] = Rbirth_pvals, _["birth_fold"] = Rbirth_folds, 
                                 _["birth_ham"] = Rbirth_hams, _["birth_qave"] = Rbirth_qaves));
}

// Returns a 16xN matrix with the observed counts of each transition categorized by
// type (row) and quality (column). Assumes qualities start at 0.
Rcpp::IntegerMatrix b_make_transition_by_quality_matrix(B *b, Sub **subs, bool has_quals, int ncol) {
  unsigned int i, r, pos0, pos1, nti0, nti1, qual, t_ij;
  Sub *sub;
  Raw *raw, *center;
  
  if(!has_quals) { ncol = 1; }
  
  // Storage for counts for 0...(ncol-1) and each nti->ntj
  Rcpp::IntegerMatrix transMat(16, ncol);

  for(i=0;i<b->nclust;i++) {
    center = b->bi[i]->center;
    for(r=0;r<b->bi[i]->nraw;r++) {
      raw = b->bi[i]->raw[r];
      if(!raw->correct) { continue; } // Don't count if not being corrected... Need to make sure used w/ care R side
      sub = subs[raw->index]; // The sub object includes the map between the center and the raw positions
      if(!sub) {
        if(VERBOSE) { Rprintf("Warning: No sub for R%i in C%i.\n", r, i); }
        continue;
      }
      
      for(pos0=0;pos0<center->length;pos0++) {
        pos1 = sub->map[pos0];
        if(pos1 == GAP_GLYPH) { // A gap in the aligned seq
          continue; // Gaps excluded from the model
        }
        nti0 = (int) (center->seq[pos0] - 1);
        nti1 = (int) (raw->seq[pos1] - 1);
        qual = raw->qual[pos1]; // qual=unsigned int
        // And record these counts
        t_ij = (4*nti0)+nti1;
        if(has_quals) {
          transMat(t_ij, qual) += raw->reads;
        } else { 
          transMat(t_ij, 0) += raw->reads; 
        }
      } // for(pos0=0;pos0<center->length;pos0++)
    } // for(r=0;b->bi[i]->nraw)
  } // for(i=0;i<b->nclust;i++)
  
  return(transMat);
}

// Makes data.frame of number of substitutions by position on the sequence
// Also finds the expected number of substitutions at each position, based on quality scores
//    and the input error matrix
Rcpp::DataFrame b_make_positional_substitution_df(B *b, Sub **subs, unsigned int seqlen, Rcpp::NumericMatrix errMat, bool use_quals) {
  unsigned int i, pos, pos1, qind, j, r, s, nti0;
  Raw *raw;
  Sub *sub;
  Rcpp::IntegerVector nts_by_pos(seqlen);
  Rcpp::IntegerVector subs_by_pos(seqlen);
  Rcpp::NumericVector exp_by_pos(seqlen);

  for(i=0;i<b->nclust;i++) {
    // Iterate through raws
    for(r=0;r<b->bi[i]->nraw;r++) {
      raw = b->bi[i]->raw[r];
      sub = subs[raw->index];
      if(sub) { // not a NULL sub
        // Add to the subs count
        for(s=0;s<sub->nsubs;s++) {
          subs_by_pos(sub->pos[s]) += raw->reads;
        }
        
        for(pos=0;pos<b->bi[i]->center->length;pos++) {
          pos1 = sub->map[pos];
          if(pos1 == GAP_GLYPH) { // Gap
            continue;
          }
          // Add to the nts_by_pos count
          nts_by_pos(pos) += raw->reads;
          // Add expected error count
          if(use_quals) {
            // Turn quality into the index in the array
            qind = raw->qual[pos1];  // qind = unsigned int
          } else {
            qind = 0;
          }
          nti0 = (int) b->bi[i]->center->seq[pos] - 1;
          for(j=4*nti0;j<4*nti0+4;j++) {
            if(j%5 == 0) { continue; } // the same-nt transitions
            exp_by_pos(pos) += raw->reads*errMat(j, qind);
          }
        }
      } // if(sub)
    }
  }
  return(Rcpp::DataFrame::create(_["nts"] = nts_by_pos, _["subs"] = subs_by_pos, _["exp"] = exp_by_pos));
}


// Calculate the average positional qualities for each cluster/partition/Bi
// Return position (rows) by Bi (columns) matrix.
Rcpp::NumericMatrix b_make_cluster_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int maxlen) {
  unsigned int i, r, pos0, pos1, raw_reads, seqlen;
  std::vector<unsigned int> nreads(maxlen);
  Sub *sub;
  Raw *raw;
  Rcpp::NumericMatrix Rquals(maxlen, b->nclust);
  
  if(has_quals) {
    for(i=0;i<b->nclust;i++) {
      seqlen = b->bi[i]->center->length;
      for(pos0=0;pos0<seqlen;pos0++) { nreads[pos0] = 0; }
      for(r=0;r<b->bi[i]->nraw;r++) {
        raw = b->bi[i]->raw[r];
        if(!raw->correct) { continue; }
        raw_reads = raw->reads;
        sub = subs[raw->index];
        if(sub) {
          for(pos0=0;pos0<seqlen;pos0++) {
            pos1 = sub->map[pos0];
            if(pos1 == GAP_GLYPH) { // Gap
              continue;
            }
            nreads[pos0] += raw_reads;
            Rquals(pos0,i) += (raw->qual[pos1] * raw_reads);  // Output qual construction
          }
        }
      } // for(pos0=0;pos0<len1;pos0++)
      for(pos0=0;pos0<seqlen;pos0++) { Rquals(pos0,i) = Rquals(pos0,i)/nreads[pos0]; }
      for(pos0=seqlen;pos0<maxlen;pos0++) { Rquals(pos0,i) = NA_REAL; }
    }
  }
  
  return(Rquals);
}

// Make output data.frame of birth_subs
Rcpp::DataFrame b_make_birth_subs_df(B *b, Sub **birth_subs, bool has_quals) {
  Sub *sub;
  unsigned int i, j, s;
  unsigned int tot_subs=0;
  char buf[2] = {'\0','\0'};
  // Count total number of subs
  for(i=0;i<b->nclust;i++) {
    sub = birth_subs[i];
    if(sub) { tot_subs += sub->nsubs; }
  }
  // Initialize the columns
  Rcpp::IntegerVector bs_pos(tot_subs);
  std::vector<std::string> bs_nt0(tot_subs);
  std::vector<std::string> bs_nt1(tot_subs);
  Rcpp::NumericVector bs_qual(tot_subs);
  Rcpp::IntegerVector bs_clust(tot_subs);
  // Record data for each birth substitution
  for(i=0,j=0;i<b->nclust;i++) {
    sub = birth_subs[i];
    if(sub) {
      for(s=0;s<sub->nsubs;s++) {
        bs_pos[j] = sub->pos[s]+1; // R 1 indexing
        buf[0] = sub->nt0[s];
        int2nt(buf, buf);
        bs_nt0[j].assign(std::string(buf));
        buf[0] = sub->nt1[s];
        int2nt(buf, buf);
        bs_nt1[j].assign(std::string(buf));
        if(has_quals) {
          bs_qual[j] = sub->q1[s];  // NumericVector, making birth sub quals from sub object
        } else {
          bs_qual[j] = Rcpp::NumericVector::get_na();
        }
        bs_clust[j] = i+1; // R 1 indexing
        j++;
      }
    }
  }
  return(Rcpp::DataFrame::create(_["pos"] = bs_pos, _["ref"] = bs_nt0, _["sub"] = bs_nt1, _["qual"] = bs_qual, _["clust"] = bs_clust));
}
