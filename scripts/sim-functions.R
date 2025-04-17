#### estimation error ####
est.err <- function(C, Chat) {
  err <- norm(C - Chat, "F") ^ 2 / (norm(C, "F") ^ 2)
  return(err)
}

#### prediction error ####
pred.err <- function(X, C, Chat) {
  n <- nrow(X)
  err <- norm(X %*% (C - Chat), "F") ^ 2 / (norm(X %*% C, "F") ^ 2)
  return(err)
}

#### FPR and FNR ####
frate <- function(U, V, Uhat, Vhat) {
  f <- function(x) {
    norm(x, "2")
  }
  u_rows <- apply(U, 1, f)
  v_rows <- apply(V, 1, f)
  v_hat_rows <- apply(Vhat, 1, f)
  u_hat_rows <- apply(Uhat, 1, f)
  fp <-
    sum(u_rows == 0 &
          u_hat_rows != 0) + sum(v_rows == 0 & v_hat_rows != 0)
  tp <-
    sum(u_rows != 0 &
          u_hat_rows != 0) + sum(v_rows != 0 & v_hat_rows != 0)
  fn <-
    sum(u_rows != 0 &
          u_hat_rows == 0) + sum(v_rows != 0 & v_hat_rows == 0)
  tn <-
    sum(u_rows == 0 &
          u_hat_rows == 0) + sum(v_rows == 0 & v_hat_rows == 0)
  
  return(list(
    fprate = fp / (fp + tn) * 100,
    fnrate = fn / (fn + tp) * 100
  ))
}

# main function for simulations
sim.table <-
  function(n,
           p,
           q,
           snr,
           corrE,
           d = c(10, 15, 20),
           only.cospa = FALSE) {
    data <-
      init.data(
        n = n,
        p = p,
        q = q,
        d = d,
        snr = snr,
        corrX = "AR",
        corrE = corrE
      )
    Y <- data$Y
    X <- data$X
    C <- data$C
    U <- data$U
    V <- data$V
    d <- data$d
    nrank <- data$nrank
    nlambda <- 200
    rank.max <- 10
    s.max <- 20
    s.min <- 1
    
    rrr.fit <-
      cospa.fit <- secure.fit <- rssvd.fit <- srrr.fit <- data.frame(
        ErC = 0,
        ErXC = 0,
        FPR = 0,
        FNR = 0,
        Rank = 0,
        Time = 0
      )
    
    # cospa
    fit <- cospa(X,
                 Y,
                 nrank = rank.max,
                 smin = s.min,
                 smax = s.max)
    Uhat <- as.matrix(fit$U)
    Vhat <- as.matrix(fit$V)
    Dhat <- fit$d
    if (length(Dhat) > 1) {
      Chat <- Uhat %*% diag(Dhat) %*% t(Vhat)
    } else {
      Chat <- Dhat[1] * Uhat %*% t(Vhat)
    }
    cospa.fit$ErC <- est.err(C, Chat)
    cospa.fit$ErXC <- pred.err(X, C, Chat)
    cospa.fit$Rank <- fit$nrank
    rate <- frate(U, V, Uhat, Vhat)
    cospa.fit$FPR <- rate$fprate
    cospa.fit$FNR <- rate$fnrate
    cospa.fit$Time <- max(fit$time)
    
    
    if (!only.cospa) {
      # update erank for orther methods
      erank <- fit$nrank
      
      # rrr
      start_time <- Sys.time()
      Chat <- ginv(t(X) %*% X) %*% t(X) %*% Y
      fit <- svd(X %*% Chat, nu = nrank, nv = nrank)
      Vhat <- fit$v
      Dhat <- fit$d[1:nrank]
      Uhat <- Chat %*% Vhat
      Chat <- Uhat %*% t(Vhat)
      for (i in 1:nrank) {
        Dhat[i] <- norm(X %*% Uhat[, i] / sqrt(nrow(X)), "2")
        Uhat[, i] <- Uhat[, i] / Dhat[i]
      }
      end_time <- Sys.time()
      rrr.fit$ErC <- est.err(C, Chat)
      rrr.fit$ErXC <- pred.err(X, C, Chat)
      rate <- frate(U, V, Uhat, Vhat)
      rrr.fit$FPR <- rate$fprate
      rrr.fit$FNR <- rate$fnrate
      rrr.fit$Rank <- ncol(Vhat)
      rrr.fit$Time <- difftime(end_time, start_time, units = "secs")
      
      # secure(D is diagonal matrix decrease)
      start_time <- Sys.time()
      fit <- secure.path(Y, X, nrank = erank, nlambda = nlambda)
      end_time <- Sys.time()
      Chat <- fit$C.est
      Uhat <- as.matrix(fit$U)
      Vhat <- as.matrix(fit$V)
      Dhat <- diag(fit$D)
      secure.fit$ErC <- est.err(C, Chat)
      secure.fit$ErXC <- pred.err(X, C, Chat)
      rate <- frate(U, V, Uhat, Vhat)
      secure.fit$FPR <- rate$fprate
      secure.fit$FNR <- rate$fnrate
      secure.fit$Rank <- ncol(Vhat)
      secure.fit$Time <-
        difftime(end_time, start_time, units = "secs")
      
      # rssvd (d is vector, decrease)
      start_time <- Sys.time()
      fit <-
        rssvd(Y,
              X,
              nrank = erank,
              control = list(nlambda = nlambda))
      end_time <- Sys.time()
      Chat <- fit$U %*% diag(fit$D, nrow = length(fit$D)) %*% t(fit$V)
      Uhat <- as.matrix(fit$U)
      Vhat <- as.matrix(fit$V)
      Dhat <- fit$D
      rssvd.fit$ErC <- est.err(C, Chat)
      rssvd.fit$ErXC <- pred.err(X, C, Chat)
      rate <- frate(U, V, Uhat, Vhat)
      rssvd.fit$FPR <- rate$fprate
      rssvd.fit$FNR <- rate$fnrate
      rssvd.fit$Rank <- ncol(Vhat)
      rssvd.fit$Time <- difftime(end_time, start_time, units = "secs")
      
      # srrr (D diagonal matrix, may not increase)
      start_time <- Sys.time()
      fit <-
        srrr(
          Y,
          X,
          nrank = erank,
          method = "adglasso",
          modstr = list(nlam = nlambda)
        )
      end_time <- Sys.time()
      Chat <- fit$U %*% fit$D %*% t(fit$V)
      Uhat <- as.matrix(fit$U)
      Vhat <- as.matrix(fit$V)
      Dhat <- diag(fit$D)
      srrr.fit$ErC <- est.err(C, Chat)
      srrr.fit$ErXC <- pred.err(X, C, Chat)
      rate <- frate(U, V, Uhat, Vhat)
      srrr.fit$FPR <- rate$fprate
      srrr.fit$FNR <- rate$fnrate
      srrr.fit$Rank <- ncol(Vhat)
      srrr.fit$Time <- difftime(end_time, start_time, units = "secs")
      
      table <-
        rbind(rrr.fit, srrr.fit, secure.fit, rssvd.fit, cospa.fit)
      table <-
        cbind(table, data.frame(method = c("RRR", "SRRR", "SeCURE", "IEEA", "CoSPA")))
    } else {
      table <- cospa.fit
    }
    
    
    return(table)
  }