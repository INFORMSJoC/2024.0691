#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "cospa_unit.h"

using namespace std;

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
