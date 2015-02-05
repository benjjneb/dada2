#include <Rcpp.h>
#include "dada.h"

// [[Rcpp::interfaces(r, cpp)]]

int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

// [[Rcpp::export]]
Rcpp::DataFrame getProbs(Rcpp::NumericMatrix err, std::vector<int> nnt, int maxD) {
  // THIS DOESNT DO ANY SORTING YET, NOT DONE
  int i, j, k, d;
  size_t index;
  
  // Copy err "matrix" into errs vector
  k=0;
  double errs[NERRS];
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      if(i!=j) { errs[k++] = err(i,j); }
    }
  }

  std::vector<Prob> probs;

  // Get prob of error for each nt, and the self-trans prob
  double pa = errs[0]+errs[1]+errs[2];
  double pc = errs[3]+errs[4]+errs[5];
  double pg = errs[6]+errs[7]+errs[8];
  double pt = errs[9]+errs[10]+errs[11];
  double self = pow((1.-pa), nnt[0]) * pow((1.-pc), nnt[1]) * pow((1.-pg), nnt[2]) * pow((1.-pt), nnt[3]);

  // make errors relative to the non-err prob
  for(i=0;i<3;i++) { errs[i] = errs[i]/(1.0 - pa); }
  for(i=3;i<6;i++) { errs[i] = errs[i]/(1.0 - pc); }
  for(i=6;i<9;i++) { errs[i] = errs[i]/(1.0 - pg); }
  for(i=9;i<12;i++) { errs[i] = errs[i]/(1.0 - pt); }
  
  // Declare variables needed to iterate though all d-aways
  int nerr[NERRS];
  double nopen[4];
  int first, store;
  double p, n;

  for(d=0;d<=maxD;d++) {
    // init partition
    for(i=0;i<NERRS;i++) { nerr[i] = 0; }
    nerr[NERRS-1] = d;
    first = NERRS-1;
    if(VERBOSE) {
      Rcpp::Rcout << "---- D = " << d << " ----\n";
    }

    while(1) {
      // Calc and store p/n/pval contribution
      p = self;
      n = 1;
      for(i=0;i<4;i++) { nopen[i] = nnt[i]; }
      
      for(i=0;i<NERRS;i++) {
        if(nerr[i]>0) {
          for(j=0;j<nerr[i];j++) {
            p = p*errs[i];
            n = n*nopen[i/3]/(j+1.0);
            nopen[i/3]--; // one fewer of that base available for future errors
            // going below zero no prob, since the times zero makes everything zero
            // IS IT THE EXCHANGEABILITY OF THE BASE DESTINATION CLASSES?
          }
        }
      }
      probs.push_back(Prob(p,n));

      if(nerr[0] >= d) { break; } // Should break when all d in first
      
      // Advance to next partition
      if(first > 0) {
        nerr[first]--;
        first--;
        nerr[first]++;
      } else {
        store = nerr[0];
        nerr[0] = 0;
        for(first=1;first<NERRS;first++) {
          if(nerr[first]>0) { break; }
          if(first==(NERRS-1)) {  // Error checking
            printf("Bad news!\n");
            return R_NilValue;
          }
        }
        nerr[first]--;
        first--;
        nerr[first] += (store+1);
      }
    } // while(1)

  } // for(d=0;d<maxD;d++)

  // Sort probs in descending order and create pval vector
  std::sort(probs.rbegin(), probs.rend());   // note: reverse iterators -> decreasing sort
  // For some reason this (very slightly) changes the cumulative sum (~10^-15 difference)!?
  // Could be a precision issue? Yep, double's have 15-17 digit significand precision.

  std::vector<double> ps;
  std::vector<double> ns;
  std::vector<double> cdf;
  double cum = 0;
  for(index=0;index < probs.size();index++) {
    cum += (probs[index].first * probs[index].second);
    ps.push_back(probs[index].first);
    ns.push_back(probs[index].second);
    cdf.push_back(cum);
  }

  // Find largest p beyond maxD (= max_err**(maxD+1))
  double max_err = 0.0;
  for(i=0;i<12;i++) { if(errs[i] > max_err) max_err = errs[i]; }
  double min_p = pow(max_err, maxD+1);
  
  for(index=0;index < ps.size();index++) {
    if(ps[index] <= min_p) { break; }
  }
  ps.resize(index);
  ns.resize(index);
  cdf.resize(index);

  return Rcpp::DataFrame::create(Rcpp::_["p"]=ps, Rcpp::_["n"]=ns, Rcpp::_["cdf"]=cdf);
}


Rcpp::NumericVector sampleProbs(std::vector<double> errs, std::vector<int> nnt, int32_t nsam) {
  // binom draw of #A errors, #C errors, #G errors and #T errors
  int n, e;
  double self, p;
  Rcpp::RNGScope scope;
  
  // Get prob of error for each nt, and the self-trans prob
  double pa = errs[0]+errs[1]+errs[2];
  double pc = errs[3]+errs[4]+errs[5];
  double pg = errs[6]+errs[7]+errs[8];
  double pt = errs[9]+errs[10]+errs[11];
  double cp_ac = errs[0]/pa;  // conditional prob of A->C given error at an A
  double cp_acg = (errs[0]+errs[1])/pa;  // conditional prob of A->C or A->G given error at an A
  double cp_ca = errs[3]/pc;
  double cp_cag = (errs[3]+errs[4])/pc;
  double cp_ga = errs[6]/pg;
  double cp_gac = (errs[6]+errs[7])/pg;
  double cp_ta = errs[9]/pt;
  double cp_tac = (errs[9]+errs[10])/pt;
  self = pow((1.-pa), nnt[0]) * pow((1.-pc), nnt[1]) * pow((1.-pg), nnt[2]) * pow((1.-pt), nnt[3]);
  
  // sample the number of errors from each nucleotide
  Rcpp::NumericVector n_Aerr = Rcpp::rbinom(nsam, nnt[0], pa);
  Rcpp::NumericVector n_Cerr = Rcpp::rbinom(nsam, nnt[1], pc);
  Rcpp::NumericVector n_Gerr = Rcpp::rbinom(nsam, nnt[2], pg);
  Rcpp::NumericVector n_Terr = Rcpp::rbinom(nsam, nnt[3], pt);
  
  // DOUBLE CHECK THAT THESE ARE UNAMBIGUOUSLY COMPARED TO INTEGERS CORRECTLY
  // IE 2 < (RV=2.) WORKS CORRECTLY EVERY TIME

  // sample a uniform number for each error
  size_t sumErr = 0;
  for(n=0;n<nsam;n++) { sumErr = sumErr + n_Aerr[n] + n_Cerr[n] + n_Gerr[n] + n_Terr[n]; }
  Rcpp::NumericVector rvs = Rcpp::runif(sumErr);

  size_t ii = 0;
  for(n=0;n<nsam;n++) {
    p = self;
    for(e=0;e<n_Aerr[n];e++) {
      if(rvs[ii]<cp_ac) {
        p *= errs[0];  // !!!!!PROBLEM, NOT PROPERLY DIVIDING BY THE NON-ERROR PROB!!
      } else if (rvs[ii] < cp_acg) {
        p *= errs[1];
      } else if (rvs[ii] <= 1.) {
        p *= errs[2];
      } else {
        printf("BAD UNIFORM RANDOM NUMBER!\n");
      }
      ii++;
    }
    /// REPEAT REPEAT REPAAT
  }
}






