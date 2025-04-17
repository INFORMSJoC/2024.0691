#ifndef __COSPA_UNIT__
#define __COSPA_UNIT__

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

double cospa_unit(double &d,
                  arma::vec &u,
                  arma::vec &v,
                  arma::uvec &aset_u,
                  arma::uvec &iset_u,
                  arma::uvec &aset_v,
                  arma::uvec &iset_v,
                  arma::mat &XX,
                  arma::mat &XY,
                  double tol = 1e-4,
                  unsigned max_iter = 200,
                  bool verbose = false);
#endif
