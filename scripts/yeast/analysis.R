rm(list = ls())
cat('\f')
require(foreach)
require(doParallel)

path <- "~/"
path <- paste0(path, "2024.0691/scripts/")
setwd(path)

data.path <- "../data/"
result.path <- "../results/raw/"
dir.create(result.path, showWarnings = FALSE, recursive = TRUE)

# repetition and number of cores
repetition <- 200
cl.num <- 50
cl <- makeCluster(cl.num)
registerDoParallel(cl)

load(paste0(data.path, "yeast_preprocess_data.RData"))

pkg.list <- c('cospa', 'rrpack', 'secure')

table <- foreach(
  i = 1:repetition,
  .combine = 'rbind',
  .packages = pkg.list
) %dopar% {
  set.seed(i)
  
  # choose 20% as out of sample samples
  train_ind <-
    sample(1:nrow(X), floor(nrow(X) * 0.2), replace = FALSE)
  X_test <- X[train_ind, ]
  Y_test <- Y[train_ind, ]
  X_train <- X[-train_ind, ]
  Y_train <- Y[-train_ind, ]
  
  n_train <- nrow(X_train)
  p_train <- ncol(X_train)
  q_train <- ncol(Y_train)
  
  n_test <- nrow(X_test)
  p_test <- ncol(X_test)
  q_test <- ncol(Y_test)
  
  pre.rank <- 4
  
  cospa.fit <- secure.fit <- rssvd.fit <- data.frame(
    U20 = 0,
    V20 = 0,
    MSE = 0,
    Prediction = 0,
    Rank = 0,
    Time = 0
  )
  
  eval <- function(fit, coef) {
    c(
      sum(apply(fit$U, 1, function(a)
        sum(a ^ 2)) != 0),
      sum(apply(fit$V, 1, function(a)
        sum(a ^ 2)) != 0),
      (norm(Y_train - X_train %*% coef, "F")) ^ 2 / (n_train * q_train),
      (norm(Y_test - X_test %*% coef, "F")) ^ 2 / (n_test * q_test)
    )
  }
  
  # Fit secure
  start_time <- Sys.time()
  fit <-
    secure.path(Y_train,
                X_train,
                nrank = pre.rank,
                control = secure.control(spU = 6 * 1e-3))
  end_time <- Sys.time()
  coef <- fit$C.est
  secure.fit$Time <- difftime(end_time, start_time, units = "secs")
  secure.fit$Rank <- length(diag(fit$D))
  secure.fit[1:4] <- eval(fit, coef)
  
  # Fit IEEA
  start_time <- Sys.time()
  fit <- rssvd(Y_train, X_train, nrank = secure.fit$Rank)
  end_time <- Sys.time()
  coef <- fit$U %*% diag(fit$D, nrow = length(fit$D)) %*% t(fit$V)
  rssvd.fit$Time <- difftime(end_time, start_time, units = "secs")
  rssvd.fit$Rank <- secure.fit$Rank
  rssvd.fit[1:4] <- eval(fit, coef)
  
  # Fit CoSPA
  fit <-
    cospa(
      X = X_train,
      Y = Y_train,
      smin = 3,
      smax = 15,
      nrank = pre.rank
    )
  cospa.fit$Time <- fit$time
  cospa.fit$Rank <- fit$nrank 
  if (length(fit$d) > 1) {
    coef <- fit$U %*% diag(fit$d) %*% t(fit$V)
  } else {
    coef <- fit$d * fit$U %*% t(fit$V)
  }
  cospa.fit[1:4] <- eval(fit, coef)
  
  res <- rbind(secure.fit, rssvd.fit, cospa.fit)
  tmp <- data.frame(method = c("SeCURE", "IEEA", "CoSPA"))
  res <- cbind(res, tmp)
  
  res
}

# end parallel computing
stopImplicitCluster()
stopCluster(cl)

filename <- paste0(result.path, "table-real-screening.RData")
save(table, file = filename)
