rm(list = ls())
cat('\f')
require(ggplot2)
require(egg)

path <- "~/"
path <- paste0(path, "2024.0691/scripts/")
setwd(path)

data.path <- "../results/raw/sim-vary/"
result.path <- "../results/figures/"
dir.create(result.path, showWarnings = FALSE, recursive = TRUE)

er.scale <- 1e3
# vary p
file <- "vary-p"
load(paste0(data.path, file, ".RData"))
table$ErC <- table$ErC * er.scale
table$ErXC <- table$ErXC * er.scale
table$p <- as.character(table$p)

plot.ErC <- ggplot(table, aes(x = p, y = ErC)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab(expression(paste("Er(C)"%*%"10"^"3"))) +
  ylim(c(0, 25)) +
  theme(axis.title.x = element_blank())
plot.ErXC <- ggplot(table, aes(x = p, y = ErXC)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab(expression(paste("Er(XC)"%*%"10"^"3"))) +
  ylim(c(0, 25)) +
  theme(axis.title.x = element_blank())
plot.fpr <- ggplot(table, aes(x = p, y = FPR)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("FPR (%)") +
  ylim(c(0, 8)) +
  theme(axis.title.x = element_blank())
plot.fnr <- ggplot(table, aes(x = p, y = FNR)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("FNR (%)") +
  ylim(c(0, 3)) +
  theme(axis.title.x = element_blank())
plot.rank <- ggplot(table, aes(x = p, y = Rank)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("Rank") +
  ylim(c(0, 6)) +
  theme(axis.title.x = element_blank())

pdf(
  file = paste0(result.path, file, ".pdf"),
  width = 10,
  height = 10 / 5.291,
  onefile = FALSE
)
ggarrange(plot.ErC,
          plot.ErXC,
          plot.fpr,
          plot.fnr,
          plot.rank,
          nrow = 1,
          ncol = 5)
dev.off()

# vary q
file <- "vary-q"
load(paste0(data.path, file, ".RData"))
table$ErC <- table$ErC * er.scale
table$ErXC <- table$ErXC * er.scale
table$q <- as.character(table$q)

plot.ErC <- ggplot(table, aes(x = q, y = ErC)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab(expression(paste("Er(C)"%*%"10"^"3"))) +
  ylim(c(0, 25)) +
  theme(axis.title.x = element_blank())
plot.ErXC <- ggplot(table, aes(x = q, y = ErXC)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab(expression(paste("Er(XC)"%*%"10"^"3"))) +
  ylim(c(0, 25)) +
  theme(axis.title.x = element_blank())
plot.fpr <- ggplot(table, aes(x = q, y = FPR)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("FPR (%)") +
  ylim(c(0, 8)) +
  theme(axis.title.x = element_blank())
plot.fnr <- ggplot(table, aes(x = q, y = FNR)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("FNR (%)") +
  ylim(c(0, 3)) +
  theme(axis.title.x = element_blank())
plot.rank <- ggplot(table, aes(x = q, y = Rank)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("Rank") +
  ylim(c(0, 6)) +
  theme(axis.title.x = element_blank())

pdf(
  file = paste0(result.path, file, ".pdf"),
  width = 10,
  height = 10 / 5.291,
  onefile = FALSE
)
ggarrange(plot.ErC,
          plot.ErXC,
          plot.fpr,
          plot.fnr,
          plot.rank,
          nrow = 1,
          ncol = 5)
dev.off()

# vary snr
file <- "vary-snr"
load(paste0(data.path, file, ".RData"))
table$ErC <- table$ErC * er.scale
table$ErXC <- table$ErXC * er.scale
table$snr <- as.character(table$snr)

plot.ErC <- ggplot(table, aes(x = snr, y = ErC)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab(expression(paste("Er(C)"%*%"10"^"3"))) +
  ylim(c(0, 25)) +
  theme(axis.title.x = element_blank())
plot.ErXC <- ggplot(table, aes(x = snr, y = ErXC)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab(expression(paste("Er(XC)"%*%"10"^"3"))) +
  ylim(c(0, 25)) +
  theme(axis.title.x = element_blank())
plot.fpr <- ggplot(table, aes(x = snr, y = FPR)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("FPR (%)") +
  ylim(c(0, 8)) +
  theme(axis.title.x = element_blank())
plot.fnr <- ggplot(table, aes(x = snr, y = FNR)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("FNR (%)") +
  ylim(c(0, 3)) +
  theme(axis.title.x = element_blank())
plot.rank <- ggplot(table, aes(x = snr, y = Rank)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_grey(start = 0.4, end = 0.8) +
  ylab("Rank") +
  ylim(c(0, 6)) +
  theme(axis.title.x = element_blank())

pdf(
  file = paste0(result.path, file, ".pdf"),
  width = 10,
  height = 10 / 5.291,
  onefile = FALSE
)
ggarrange(plot.ErC,
          plot.ErXC,
          plot.fpr,
          plot.fnr,
          plot.rank,
          nrow = 1,
          ncol = 5)
dev.off()