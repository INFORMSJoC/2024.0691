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
