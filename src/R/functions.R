#' initialize the data for simulation
#' @param n sample size of the observation
#' @param p dimension of predictors
#' @param q dimension of responses 
#' @param d singular values of the coefficient matrix
#' @param su sparsity level of left singular vectors 
#' @param sv sparsity level of right singular vectors
#' @param snr signal-to-noise rate for the coefficient matrix 
#' @param corrX covariance matrix type of the predictors
#' @param rhoX the parameter for the predictor covariance matrix 
#' @param corrE covariance matrix type of the noise vector
#' @param rhoE the parameter for the noise covariance matrix
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm 
#' @importFrom stats runif
#' @export 
init.data <- function(n = 100,
                      p = 200,
                      q = 100,
                      d = c(10, 15, 20),
                      su = 3,
                      sv = 4,
                      snr = 0.5,
                      corrX = c("AR", "CS"),
                      rhoX = 0.5,
                      corrE = c("AR", "I", "CS"),
                      rhoE = 0.3) {
  nrank <- length(d)
  if (corrX == "AR") {
    SigmaX <- rhoX ^ abs(outer(1:p, 1:p, "-"))
  }
  else {
    SigmaX <- matrix(rhoX, p, p) + diag(1 - rhoX, p)
  }

  if (corrE == "I") {
    E <- matrix(rnorm(n * q), nrow = n, ncol = q)
  }
  else {
    if (corrE == "AR") {
      SigmaE <- rhoE ^ abs(outer(1:q, 1:q, "-"))
    }
    else {
      SigmaE <- matrix(rhoE, q, q) + diag(1 - rhoE, q)
    }
    E <- mvrnorm(n, rep(0, q), SigmaE)
  }

  X <- mvrnorm(n, rep(0, p), SigmaX)

  U <- matrix(0, nrow = p, ncol = nrank)
  V <- matrix(0, nrow = q, ncol = nrank)
  for (i in 1:nrank) {
    U[, i] <-
      c(rep(0, (i - 1) * su), sample(c(1, -1), su, replace = TRUE), rep(0, p - i * su))
    V[, i] <-
      c(rep(0, (i - 1) * sv),
        sample(c(1, -1), sv, replace = TRUE) * runif(sv, 0.3, 1),
        rep(0, q - i * sv))
    U[, i] <- U[, i] / norm(U[, i], "2")
    V[, i] <- V[, i] / norm(V[, i], "2")
  }
  if (nrank == 1) {
    C <- U %*% t(d * V)
  } else {
    C <- U %*% diag(d) %*% t(V)
  }

  res.svd <- svd(X %*% C, nu = nrank, nv = nrank)
  V <- res.svd$v
  UD <- C %*% V
  for (i in 1:nrank) {
    d[i] <- norm(X %*% UD[, i] / sqrt(n), "2")
    U[, i] <- UD[, i] / d[i]
  }
  res <- sort(d, decreasing = TRUE, index.return = TRUE)
  U <- U[, res$ix]
  V <- V[, res$ix]
  d <- res$x

  sigma <- res.svd$d[nrank] / norm(E, "F") / snr
  E <- sigma * E

  Y <- X %*% C + E
  return(list(
    X = X,
    Y = Y,
    C = C,
    U = U,
    V = V,
    d = d,
    nrank = nrank
  ))
}

#' estimation function for the co-sparse matrix 
#' @param X the predictor matrix 
#' @param Y the response matrix 
#' @param smin the minimum sparsity level 
#' @param smax the maximum sparsity level 
#' @param nrank the rank of the coefficient matrix
#' @param alpha the paramter to tune the penalty level, default is 1
#' @param verbose indicator to print iteration procedure 
#' @return A list includes the estimation result
#' \itemize{
#'   \item d - the singular values 
#'   \item U - the left singular vectors 
#'   \item V - the right singular vectors 
#'   \item time - the computation time 
#'   \item nrank - the estimated rank 
#' }
#' @importFrom glmnet glmnet
#' @importFrom scalreg scalreg
#' @importFrom stats coef
#' @useDynLib cospa
#' @export 
cospa <- function(X, Y, smin, smax, nrank, alpha = 1, verbose = FALSE) {
  p <- ncol(X)
  q <- ncol(Y)
  t1 <- Sys.time()
  fit <- scalreg(X, Y[, 1])
  lambda <- max(abs(fit$residuals)) * 1e-2
  C.init <- matrix(0, p, q)
  for (i in c(1:q)) {
    fit <- glmnet(X, Y[, i], lambda = lambda, intercept = FALSE, standardize = FALSE)
    C.init[, i] <- as.vector(coef(fit))[-1]
  }
  delta_t <- difftime(Sys.time(), t1, units = "secs")
  out <- cospa_para(X, Y, C.init, nrank, smin, smax, alpha, verbose = verbose)
  d <- as.vector(out$d)
  time <- max(out$time) + delta_t
  U <- out$U
  V <- out$V
  nrank <- out$nrank
  if (verbose) {
    cat("Estimated rank:", out$nrank, "\n")
  }
  return(list(d = d, U = U, V = V, time = time, nrank = nrank))
}