require(foreach)
require(doParallel)

path <- "~/"
path <- paste0(path, "2024.0691/scripts/")
setwd(path)

result.path <- "../results/raw/sim-time/"
dir.create(path, showWarnings = FALSE, recursive = TRUE)
dir.create(result.path, showWarnings = FALSE, recursive = TRUE)
source(paste0(path, "sim-functions.R"))

# simulation settings
d <- c(10, 15, 20)

# functions list for simulations
func.list <- c('est.err', 'frate', 'pred.err', 'sim.table')
pkg.list <- c('cospa')

# repetition and number of cores
repetition <- 200
cl.num <- 50
cl <- makeCluster(cl.num)
registerDoParallel(cl)

# vary n in c(100,200,300,400)
table <- data.frame()
for (n in c(100, 300, 500, 700, 900)) {
  temp.table <-
    foreach(
      i = 1:repetition,
      .combine = 'rbind',
      .packages = pkg.list,
      .export = func.list
    ) %dopar% {
      set.seed(i)
      time.result <-
        sim.table(
          n = n,
          p = 100,
          q = 100,
          snr = 0.5,
          corrE = "AR",
          d = d,
          only.cospa = TRUE
        )
      time.result$Time
    }
  temp <-
    data.frame(Time = temp.table, n = rep(n, nrow(temp.table)))
  rownames(temp) <- NULL
  table <- rbind(table, temp)
}
filename <- paste0(result.path, "time-n", ".RData")
save(table, file = filename)

# vary p in c(100,200,300,400)
table <- data.frame()
for (p in c(100, 300, 500, 700, 900)) {
  temp.table <-
    foreach(
      i = 1:repetition,
      .combine = 'rbind',
      .packages = pkg.list,
      .export = func.list
    ) %dopar% {
      set.seed(i)
      time.result <-
        sim.table(
          n = 100,
          p = p,
          q = 100,
          snr = 0.5,
          corrE = "AR",
          d = d,
          only.cospa = TRUE
        )
      time.result$Time
    }
  temp <-
    data.frame(Time = temp.table, p = rep(p, nrow(temp.table)))
  rownames(temp) <- NULL
  table <- rbind(table, temp)
}
filename <- paste0(result.path, "time-p", ".RData")
save(table, file = filename)

# vary q in c(100,200,300,400)
table <- data.frame()
for (q in c(100, 300, 500, 700, 900)) {
  temp.table <-
    foreach(
      i = 1:repetition,
      .combine = 'rbind',
      .packages = pkg.list,
      .export = func.list
    ) %dopar% {
      set.seed(i)
      time.result <-
        sim.table(
          n = 100,
          p = 100,
          q = q,
          snr = 0.5,
          corrE = "AR",
          d = d,
          only.cospa = TRUE
        )
      time.result$Time
    }
  temp <-
    data.frame(Time = temp.table, q = rep(q, nrow(temp.table)))
  rownames(temp) <- NULL
  table <- rbind(table, temp)
}
filename <- paste0(result.path, "time-q", ".RData")
save(table, file = filename)

# end parallel computing
stopImplicitCluster()
stopCluster(cl)
