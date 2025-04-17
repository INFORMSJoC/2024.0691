#include <chrono>
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

inline arma::vec col_norm(arma::mat A)
{
  arma::vec out = arma::zeros<arma::vec>(A.n_cols);
  for (unsigned i = 0; i < A.n_cols; i++)
  {
    out.at(i) = arma::norm(A.col(i), 2);
  }
  return out;
}

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
                  bool verbose = false)
{
  double eps = INFINITY;
  unsigned iteration = 0;
  double loss = d * d / 2 - d * arma::dot(v.elem(aset_v), XY.submat(aset_u, aset_v).t() * u.elem(aset_u));
  auto start = chrono::steady_clock::now();
  while (eps > tol && iteration < max_iter)
  {
    double loss_old = loss;
    arma::uvec aset_u_old = aset_u, aset_v_old = aset_v;
    
    // update A_v and v
    arma::vec eval_v = col_norm(XY.rows(aset_u));
    arma::uvec ind_v = arma::sort_index(eval_v, "descend");
    aset_v = ind_v.subvec(0, aset_v.n_elem - 1);
    iset_v = ind_v.subvec(aset_v.n_elem, ind_v.n_elem - 1);
    arma::mat tmp_U, tmp_V;
    arma::vec tmp_d;
    svd_econ(tmp_U, tmp_d, tmp_V, XY.submat(aset_u, aset_v), "right");
    v.zeros();
    v.elem(aset_v) = tmp_V.col(0);
    
    // update A_u
    arma::vec eta_u = arma::zeros<arma::vec>(u.n_elem);
    arma::vec gamma_u = eta_u;
    if(arma::rank(XX.submat(aset_u, aset_u)) < aset_u.n_elem){
      loss = d * d / 2 - d * arma::dot(v.elem(aset_v), XY.submat(aset_u, aset_v).t() * u.elem(aset_u));
      break;
    }
    eta_u.elem(aset_u) = arma::solve(XX.submat(aset_u, aset_u), XY.submat(aset_u, aset_v) * v.elem(aset_v)); // arma::inv_sympd(XX.submat(aset_u, aset_u)) * XY.submat(aset_u, aset_v) * v.elem(aset_v);
    // update u and d
    arma::vec du = eta_u.elem(aset_u);
    u.zeros();
    d = sqrt(arma::dot(du, XX.submat(aset_u, aset_u) * du));
    u.elem(aset_u) = du / d; // normalize to sqrt n for Xu
    loss = d * d / 2 - d * arma::dot(v.elem(aset_v), XY.submat(aset_u, aset_v).t() * u.elem(aset_u));
    
    gamma_u(iset_u) = XY.submat(iset_u, aset_v) * v.elem(aset_v) - XX.submat(iset_u, aset_u) * du;
    arma::vec eval_u = arma::abs(eta_u + gamma_u);
    arma::uvec ind_u = arma::sort_index(eval_u, "descend");
    aset_u = ind_u.subvec(0, aset_u.n_elem - 1);
    iset_u = ind_u.subvec(aset_u.n_elem, ind_u.n_elem - 1);
    
    //
    aset_u = arma::sort(aset_u);
    aset_v = arma::sort(aset_v);
    bool flag = (arma::sum(aset_u - aset_u_old) == 0) && (arma::sum(aset_v - aset_v_old) == 0);
    if (abs(loss - loss_old) / XY.n_cols < tol || flag)
    {
      break;
    }
    iteration++;
  }
  auto end = chrono::steady_clock::now();
  if (verbose)
  {
    Rcpp::Rcout << "unit estimate: iteration = " << iteration << ", time = " << chrono::duration_cast<chrono::microseconds>(end - start).count() << " um" << endl;
  }
  return loss;
}


// [[Rcpp::export]]
Rcpp::List cospa_para(arma::mat &X,
                      arma::mat &Y,
                      arma::mat &C,
                      unsigned nrank,
                      const unsigned smin,
                      const unsigned smax,
                      double alpha = 1,
                      double tol = 1e-4,
                      unsigned max_iter = 200,
                      bool verbose = false)
{
    arma::mat XX = X.t() * X / X.n_rows;
    double penalty = log(X.n_cols * Y.n_cols) * log(log(X.n_rows * Y.n_cols)) / X.n_rows;
    arma::mat U, V;
    arma::vec d;
    arma::mat Yres;
    arma::svds(U, d, V, arma::sp_mat(X * C / sqrt(X.n_rows)), nrank);
    U.reset();
    U = C * V;
    for (unsigned i = 0; i < nrank; i++)
    {
        d.at(i) = norm(X * U.col(i), 2) / sqrt(X.n_rows);
        U.col(i) = U.col(i) / d.at(i);
    }
    Yres = Y - X * C;
    arma::vec time = arma::zeros<arma::vec>(nrank);

    for (unsigned k = 0; k < nrank; k++)
    {
        clock_t start = clock();
        if (verbose)
        {
            Rcpp::Rcout << "factor k = " << k << "......" << endl;
        }
        // inilization for the k-th layer
        arma::vec u = U.col(k);
        arma::vec v = V.col(k);
        double d_inner = d.at(k);
        arma::mat Y_inner = Yres + (X * (d_inner * u)) * v.t();
        arma::mat XY = X.t() * Y_inner / X.n_rows;
        const double loss_const = pow(arma::norm(Y_inner, "fro"), 2) / (2 * Y.n_rows);

        // storing path solution
        unsigned path_length = smax - smin + 1;
        arma::vec d_path = arma::zeros<arma::vec>(path_length);
        arma::mat U_path = arma::zeros<arma::mat>(X.n_cols, path_length);
        arma::mat V_path = arma::zeros<arma::mat>(Y.n_cols, path_length);
        arma::uvec aset_u_path[path_length];
        arma::uvec iset_u_path[path_length];
        arma::uvec aset_v_path[path_length];
        arma::uvec iset_v_path[path_length];
        unsigned df;

        // initilize the active sets
        arma::uvec ind_u = arma::sort_index(arma::abs(u), "descend");
        arma::uvec ind_v = arma::sort_index(arma::abs(v), "descend");
        arma::uvec aset_u = ind_u.subvec(0, smin - 1);
        arma::uvec iset_u = ind_u.subvec(smin, X.n_cols - 1);
        arma::uvec aset_v = ind_v.subvec(0, smax - 1);
        arma::uvec iset_v = ind_v.subvec(smax, Y.n_cols - 1);
        arma::uvec combine;

        // tuning su
        arma::vec gic = arma::zeros<arma::vec>(path_length);
        gic.fill(INFINITY);
        gic.at(0) = 0;
        for (unsigned i = 0; i < path_length; i++)
        {
            unsigned df = smin + i + smax;
            gic.at(i) = cospa_unit(d_inner, u, v, aset_u, iset_u, aset_v, iset_v, XX, XY, tol, max_iter, verbose) + alpha * df * penalty + loss_const;
            aset_u_path[i] = aset_u;
            iset_u_path[i] = iset_u;
            aset_v_path[i] = aset_v;
            iset_v_path[i] = iset_v;
            U_path.col(i) = u;
            V_path.col(i) = v;
            d_path.at(i) = d_inner;
            if (i != (path_length - 1))
            {
                combine = arma::join_cols(aset_u, iset_u);
                aset_u = combine.subvec(0, smin + i);
                iset_u = combine.subvec(smin + i + 1, X.n_cols - 1);
            }
        }

        arma::uword ind = gic.index_min();
        aset_u = aset_u_path[ind];
        iset_u = iset_u_path[ind];
        aset_v = aset_v_path[ind];
        iset_v = iset_v_path[ind];
        u = U_path.col(ind);
        v = V_path.col(ind);
        d_inner = d_path.at(ind);

        // tuning sv
        combine = arma::join_cols(aset_v, iset_v);
        aset_v = combine.subvec(0, smin - 1);
        iset_v = combine.subvec(smin, Y.n_cols - 1);
        gic.fill(INFINITY);
        gic.at(0) = 0;
        for (unsigned i = 0; i < path_length; i++)
        {
            df = smin + i + aset_u.n_elem;
            gic.at(i) = cospa_unit(d_inner, u, v, aset_u, iset_u, aset_v, iset_v, XX, XY, tol, max_iter, verbose) + alpha * df * penalty + loss_const;
            aset_u_path[i] = aset_u;
            iset_u_path[i] = iset_u;
            aset_v_path[i] = aset_v;
            iset_v_path[i] = iset_v;
            U_path.col(i) = u;
            V_path.col(i) = v;
            d_path.at(i) = d_inner;
            if (i != (path_length - 1))
            {
                combine = arma::join_cols(aset_v, iset_v);
                aset_v = combine.subvec(0, smin + i);
                iset_v = combine.subvec(smin + i + 1, Y.n_cols - 1);
            }
        }

        ind = gic.index_min();
        aset_u = aset_u_path[ind];
        iset_u = iset_u_path[ind];
        aset_v = aset_v_path[ind];
        iset_v = iset_v_path[ind];
        u = U_path.col(ind);
        v = V_path.col(ind);
        d_inner = d_path.at(ind);

        time.at(k) = (clock() - start) / (double)CLOCKS_PER_SEC;

        // check singular value to trucate rank
        if (gic.at(ind) > loss_const)
        {
            int k_trucated = k - 1;
            Rcpp::List out;
            if (k_trucated < 0)
            {
                out["d"] = 0;
                out["U"] = arma::sp_mat(u.zeros());
                out["V"] = arma::sp_mat(v.zeros());
                out["nrank"] = 0;
                out["time"] = 0;
            }
            else
            {
                out["d"] = d.subvec(0, k_trucated);
                out["U"] = arma::sp_mat(U.cols(0, k_trucated));
                out["V"] = arma::sp_mat(V.cols(0, k_trucated));
                out["time"] = time.subvec(0, k_trucated);
                out["nrank"] = k;
            }
            return out;
        }
        d.at(k) = d_inner;
        U.col(k) = u;
        V.col(k) = v;
    }

    Rcpp::List out;
    out["d"] = d;
    out["U"] = U;
    out["V"] = V;
    out["time"] = time;
    out["nrank"] = nrank;

    return out;
}
