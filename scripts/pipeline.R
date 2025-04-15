path <- "~/2024.0691/scripts/"

# # run all numerical simulations
# source(paste0(path, "sim-table.R"))
# source(paste0(path, "sim-time.R"))
# source(paste0(path, "sim-vary.R"))

# generate tables, line and box plots of numerical simulations
source(paste0(path, "get-table.R"))
source(paste0(path, "get-box.R"))
source(paste0(path, "get-line.R"))

# run the real-world data analysis 
source(paste0(path, "yeast/data_clean.R"))
source(paste0(path, "yeast/analysis.R"))

# generate tables of real data analysis 
source(paste0(path, "yeast/get-table-real.R"))