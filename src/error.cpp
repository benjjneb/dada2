#include <Rcpp.h>
#include "dada.h"
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// This function constructs the output "clustering" data.frame for the dada(...) function.
// This contains core and diagnostic information on each partition (or cluster, or Bi).
Rcpp::DataFrame b_make_clustering_df(B *b, Sub **subs, Sub **birth_subs, bool has_quals) {
  unsigned int i, j, r, s, cind;
  Sub *sub;
  double q_ave, tot_e;
  
  // Create output character-vector of representative sequences for each partition (Bi)
  Rcpp::CharacterVector Rseqs;
  char oseq[SEQLEN];
  for(i=0;i<b->nclust;i++) {
    ntcpy(oseq, b->bi[i]->seq);
    Rseqs.push_back(std::string(oseq));
  }
  
  // Create output vectors for other columns of clustering data.frame
  // Each has one entry for each partition (Bi)
  Rcpp::IntegerVector Rabunds(b->nclust);      // abundances
  Rcpp::IntegerVector Rzeros(b->nclust);       // n0
  Rcpp::IntegerVector Rones(b->nclust);        // n1
  Rcpp::IntegerVector Rraws(b->nclust);        // nraw
  Rcpp::IntegerVector Rfams(b->nclust);        // nfam
  Rcpp::NumericVector Rbirth_pvals(b->nclust); // pvalue at birth
  Rcpp::NumericVector Rbirth_folds(b->nclust); // fold over-abundance at birth
  Rcpp::IntegerVector Rbirth_hams(b->nclust);  // hamming distance at birth
  Rcpp::NumericVector Rbirth_es(b->nclust);    // expected number at birth
  Rcpp::CharacterVector Rbirth_types;           // DEPRECATED
  Rcpp::NumericVector Rbirth_qaves(b->nclust); // average quality of substitutions that drove birth
  Rcpp::NumericVector Rpvals(b->nclust);       // post-hoc pvalue

  // Assign values to the output vectors
  for(i=0;i<b->nclust;i++) {
    Rabunds[i] = b->bi[i]->reads;
    Rraws[i] = b->bi[i]->nraw;
    Rfams[i] = 0;
    // n0 and n1
    Rzeros[i] = 0; Rones[i] = 0;
    for(r=0;r<b->bi[i]->nraw;r++) {
      sub = subs[b->bi[i]->raw[r]->index];
      if(sub) {
        if(sub->nsubs == 0) { Rzeros[i] += b->bi[i]->raw[r]->reads; }
        if(sub->nsubs == 1) { Rones[i] += b->bi[i]->raw[r]->reads; }
      }
    }
    // Record information from the cluster's birth
    Rbirth_types.push_back(std::string(b->bi[i]->birth_type));
    if(i==0) {  // 0-clust wasn't born normally
      Rbirth_pvals[i] = Rcpp::NumericVector::get_na(); 
      Rbirth_folds[i] = Rcpp::NumericVector::get_na(); 
      Rbirth_hams[i] = Rcpp::IntegerVector::get_na(); 
      Rbirth_es[i] = Rcpp::NumericVector::get_na();
      Rbirth_qaves[i] = Rcpp::NumericVector::get_na();
    } else { 
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
            q_ave += (sub->q1[s]);
          }
          q_ave = q_ave/((double)sub->nsubs);
        }
        Rbirth_qaves[i] = q_ave;
      } else {
        Rbirth_qaves[i] = Rcpp::NumericVector::get_na();
      }
    }
    
    // Calculate post-hoc pval
    // This is not as exhaustive anymore
    tot_e = 0.0;
    for(j=0;j<b->nclust;j++) {
      if(i != j && b->bi[j]->comp_index.count(b->bi[i]->center->index) > 0) {
        cind = b->bi[j]->comp_index[b->bi[i]->center->index];
        tot_e += b->bi[j]->comp[cind].lambda * b->bi[j]->reads;
//        b->bi[j]->e[b->bi[i]->center->index];
      }
    }
    Rpvals[i] = calc_pA(1+b->bi[i]->reads, tot_e); // Add 1 because calc_pA subtracts 1 (conditional p-val)
  }
  
  return(Rcpp::DataFrame::create(_["sequence"] = Rseqs, _["abundance"] = Rabunds, _["n0"] = Rzeros, _["n1"] = Rones, _["nunq"] = Rraws, _["nfam"] = Rfams, _["pval"] = Rpvals, _["birth_type"] = Rbirth_types, _["birth_pval"] = Rbirth_pvals, _["birth_fold"] = Rbirth_folds, _["birth_ham"] = Rbirth_hams, _["birth_qave"] = Rbirth_qaves));
}

// Returns a 16xN matrix with the observed counts of each transition categorized by
// type (row) and quality (column). Assumes qualities start at 0!!
Rcpp::IntegerMatrix b_make_transition_by_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int qmax) {
  unsigned int i, r, pos0, pos1, nti0, nti1, qual, t_ij;
  int ncol;
  Sub *sub;
  Raw *raw, *center;
  
  if(has_quals) { ncol = qmax+1; }
  else { ncol = 1; }
  
  if(ncol<=0) {
    Rcpp::stop("Invalid QMAX.");
  }
  
  // Storage for counts for each qual and each nti->ntj
  Rcpp::IntegerMatrix transMat(16, ncol);

  for(i=0;i<b->nclust;i++) {
    center = b->bi[i]->center;
    for(r=0;r<b->bi[i]->nraw;r++) {
      raw = b->bi[i]->raw[r];
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
        qual = round(raw->qual[pos1]);  // Need to change round here if qsteps implemented
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
Rcpp::DataFrame b_make_positional_substitution_df(B *b, Sub **subs, unsigned int seqlen, Rcpp::NumericMatrix errMat) {
  unsigned int i, pos, pos1, qind, j, r, s, nti0, ncol;
  Raw *raw;
  Sub *sub;
  Rcpp::IntegerVector nts_by_pos(seqlen);
  Rcpp::IntegerVector subs_by_pos(seqlen);
  Rcpp::NumericVector exp_by_pos(seqlen);

  ncol = errMat.ncol();
  float prefactor = ((float) (ncol-1))/((float) QMAX-QMIN);
  float fqmin = (float) QMIN;
  
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
          if(raw->qual) {
            // Turn quality into the index in the array
            qind = round(prefactor * (raw->qual[pos1] - fqmin));
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
Rcpp::NumericMatrix b_make_cluster_quality_matrix(B *b, Sub **subs, bool has_quals, unsigned int seqlen) {
  double q_sum;
  unsigned int i, r, nreads, pos0, pos1;
  Sub *sub;
  Raw *raw;
  Rcpp::NumericMatrix Rquals(seqlen, b->nclust);
  
  if(has_quals) {
    for(i=0;i<b->nclust;i++) {
      for(pos0=0;pos0<seqlen;pos0++) {
        q_sum=0.0;
        nreads=0;
        for(r=0;r<b->bi[i]->nraw;r++) {
          raw = b->bi[i]->raw[r];
          sub = subs[raw->index];
          if(sub) {
            pos1 = sub->map[pos0];
            if(pos1 == GAP_GLYPH) { // Gap
              continue;
            }
            nreads += raw->reads;
            q_sum += (raw->qual[pos1] * raw->reads);
          }
        }
        Rquals(pos0,i) = q_sum/nreads;
      } // for(pos0=0;pos0<len1;pos0++)
    }
  }
  return(Rquals);
}

// Make output data.frame of birth_subs
Rcpp::DataFrame b_make_birth_subs_df(B *b, Sub **birth_subs, bool has_quals) {
  Sub *sub;
  unsigned int i, s;
  char buf[2] = {'\0','\0'};
  // Initialize the columns
  Rcpp::IntegerVector bs_pos;
  Rcpp::CharacterVector bs_nt0;
  Rcpp::CharacterVector bs_nt1;
  Rcpp::NumericVector bs_qual;
  Rcpp::IntegerVector bs_clust;
  // Record data for each birth substitution
  for(i=0;i<b->nclust;i++) {
    sub = birth_subs[i];
    if(sub) {
      for(s=0;s<sub->nsubs;s++) {
        bs_pos.push_back(sub->pos[s]+1); // R 1 indexing
        buf[0] = sub->nt0[s];
        int2nt(buf, buf);
        bs_nt0.push_back(std::string(buf));
        buf[0] = sub->nt1[s];
        int2nt(buf, buf);
        bs_nt1.push_back(std::string(buf));
        if(has_quals) {
          bs_qual.push_back(sub->q1[s]);
        } else {
          bs_qual.push_back(Rcpp::NumericVector::get_na());
        }
        bs_clust.push_back(i+1); // R 1 indexing
      }
    }
  }
  return(Rcpp::DataFrame::create(_["pos"] = bs_pos, _["ref"] = bs_nt0, _["sub"] = bs_nt1, _["qual"] = bs_qual, _["clust"] = bs_clust));
}