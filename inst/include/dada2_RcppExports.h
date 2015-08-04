// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef __dada2_RcppExports_h__
#define __dada2_RcppExports_h__

#include <Rcpp.h>

namespace dada2 {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("dada2", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("dada2", "dada2_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in dada2");
            }
        }
    }

    inline Rcpp::DataFrame getSingletonCDF(Rcpp::NumericMatrix err, std::vector<int> nnt, int maxD) {
        typedef SEXP(*Ptr_getSingletonCDF)(SEXP,SEXP,SEXP);
        static Ptr_getSingletonCDF p_getSingletonCDF = NULL;
        if (p_getSingletonCDF == NULL) {
            validateSignature("Rcpp::DataFrame(*getSingletonCDF)(Rcpp::NumericMatrix,std::vector<int>,int)");
            p_getSingletonCDF = (Ptr_getSingletonCDF)R_GetCCallable("dada2", "dada2_getSingletonCDF");
        }
        RObject __result;
        {
            RNGScope __rngScope;
            __result = p_getSingletonCDF(Rcpp::wrap(err), Rcpp::wrap(nnt), Rcpp::wrap(maxD));
        }
        if (__result.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (__result.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(__result).c_str());
        return Rcpp::as<Rcpp::DataFrame >(__result);
    }

}

#endif // __dada2_RcppExports_h__