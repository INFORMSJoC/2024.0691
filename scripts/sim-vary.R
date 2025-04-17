require(foreach)
require(doParallel)

path <- "~/"
path <- paste0(path, "2024.0691/scripts/")
setwd(path)

result.path <- "../results/raw/sim-vary/"
dir.create(path, showWarnings = FALSE, recursive = TRUE)
dir.create(result.path, showWarnings = FALSE, recursive = TRUE)
source(paste0(path, "sim-functions.R"))

z1 <- file(paste0(result.path, "message.Rout"), open = "wt")
sink(z1, type = "message")
z2 <- file(paste0(result.path, "output.Rout"), open = "wt")
sink(z2, type = "output", split = TRUE)

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

# vary p in c(100,200,300,400)
table <- data.frame()
for (p in c(100, 200, 300, 400)) {
  temp.table <-
    foreach(
      i = 1:repetition,
      .combine = 'rbind',
      .packages = pkg.list,
      .export = func.list
    ) %dopar% {
      set.seed(i)
      sim.table(
        n = 100,
        p = p,
        q = 100,
        snr = 0.5,
        corrE = "AR",
        d = d,
        only.cospa = TRUE
      )
    }
  temp <-
    cbind(temp.table, data.frame(p = rep(p, nrow(temp.table))))
  table <- rbind(table, temp)
}
filename <- paste0(result.path, "vary-p", ".RData")
save(table, file = filename)

# vary q in c(100,200,300,400)
table <- data.frame()
for (q in c(100, 200, 300, 400)) {
  temp.table <-
    foreach(
      i = 1:repetition,
      .combine = 'rbind',
      .packages = pkg.list,
      .export = func.list
    ) %dopar% {
      set.seed(i)
      sim.table(
        n = 100,
        p = 100,
        q = q,
        snr = 0.5,
        corrE = "AR",
        d = d,
        only.cospa = TRUE
      )
    }
  temp <-
    cbind(temp.table, data.frame(q = rep(q, nrow(temp.table))))
  table <- rbind(table, temp)
}
filename <- paste0(result.path, "vary-q", ".RData")
save(table, file = filename)

# vary snr in c(0.4, 0.5, 0.6, 0.7)
table <- data.frame()
for (snr in c(0.4, 0.5, 0.6, 0.7)) {
  temp.table <-
    foreach(
      i = 1:repetition,
      .combine = 'rbind',
      .packages = pkg.list,
      .export = func.list
    ) %dopar% {
      set.seed(i)
      sim.table(
        n = 100,
        p = 100,
        q = 100,
        snr = snr,
        corrE = "AR",
        d = d,
        only.cospa = TRUE
      )
    }
  temp <-
    cbind(temp.table, data.frame(snr = rep(snr, nrow(temp.table))))
  table <- rbind(table, temp)
}
filename <- paste0(result.path, "vary-snr", ".RData")
save(table, file = filename)

sink(type = "message")
close(z1)
sink(type = "output")
close(z2)

# end parallel computing
stopImplicitCluster()
stopCluster(cl)
