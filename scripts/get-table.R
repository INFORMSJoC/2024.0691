path <- "~/"
path <- paste0(path, "2024.0691/scripts/")
setwd(path)

data.path <- "../results/raw/sim-table/"
result.path <- "../results/tables/"
dir.create(result.path, showWarnings = FALSE, recursive = TRUE)

for (file in c("table-AR", "table-SC")) {
  load(paste0(data.path, file, ".RData"))
  sink(paste0(result.path, file, ".tex"))
  
  table <- table[table$snr == 0.5, ]
  table$Rank <- NULL
  table$Time <- as.double(table$Time)
  
  # rescale ErC and ErXC
  er.scale <- 100
  table$ErC <- table$ErC * er.scale
  table$ErXC <- table$ErXC * er.scale
  
  performance <- colnames(table)[1:5]
  
  cat("\\begin{tabular}{l", rep("c", length(performance) + 1), "} \n", sep = "")
  cat("\\toprule \n")
  
  cat(
    "& Method & Er$(\\bC)\\times ",
    er.scale,
    "$ & Er$(\\bX\\bC)\\times ",
    er.scale,
    "$ & FPR (\\%) & FNR (\\%) & Time (s) \\\\ \n",
    sep = ""
  )
  
  for (snr in unique(table$snr)) {
    cat(
      "\\midrule \n & \\multicolumn{",
      length(performance) + 1,
      "}{c}{",
      "$\\text{SNR}=",
      snr,
      "$} \\\\ \n",
      sep = ""
    )
    for (p in unique(table$p)) {
      cat("\\multirow{", length(unique(table$method)), "}{*}{$p=", p, "$}", sep = "")
      ind <- (table$snr == snr) & (table$p == p)
      for (method in unique(table$method)) {
        cat(" & ", method, " & ", sep = "")
        sub.ind <- (table$method == method) & ind
        
        ave <- std <- data.frame()
        
        temp <- table[sub.ind, 1:length(performance)]
        ave <- apply(temp, 2, mean)
        std <- apply(temp, 2, sd)
        
        ave <- round(ave, 2)
        std <- round(std, 2)
        
        for (i in 1:length(ave)) {
          cat(ave[i], sep = "")
          if (i != length(ave)) {
            cat(" & ")
          }
          else {
            cat(" \\\\ \n")
          }
        }
      }
      if (p != unique(table$p)[length(unique(table$p))]) {
        cat("\\hline \n")
      }
    }
  }
  
  cat("\\bottomrule \n")
  cat("\\end{tabular} \n")
  
  sink()
}