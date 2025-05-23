[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

<!-- # 2024.0691 -->

# Fast Association Recovery in High Dimensions by Parallel Learning

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported in the paper [Fast Association Recovery in High Dimensions by Parallel Learning](href) by Ruipeng Dong, and Canhong Wen.

**Important: This code is only be used to reproduce the results of the paper, and we will not maintain this project continuously. 
If you find any bugs, please email Ruipeng Dong (<drp@ustc.edu.cn> or <drpoct@163.com>).**

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/xxx/xxx

https://doi.org/xxx/xxx.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{cospa2025,
  author    = {Ruipeng Dong and Canhong Wen},
  publisher = {INFORMS Journal on Computing},
  title     = {Fast Association Recovery in High Dimensions by Parallel Learning},
  year      = {2025},
  doi       = {xxxx},
  note      = {Available for download at: https://github.com/INFORMSJoC/xxxx},
}  
```

## Required Packages and Tools
We use R language for the numerical simulations and real-world data analysis. To run this project, make sure you have the following R packages installed. You can install them using:

```R
install.packages("devtools")
install.packages("glmnet")
install.packages("MASS")
install.packages("scalreg")
install.packages("rrpack")
install.packages("secure")
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("foreach")
install.packages("doParallel")
install.packages("egg")
install.packages("ggplot2")
```
> **Note**: Besides the above packages, we also need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed to compile our R package <code>cospa</code>.

## Folders Organization

- folder <code>data</code>: an expression quantitative trait loci data for the real-world data analysis
- folder <code>results</code>:
  - subfolder <code>figures</code>: including all figures of the paper
  - subfolder <code>raw</code>: including all raw results of simulations and real-world data analysis that can be converted to the figures and tables in our paper
  - subfolder <code>tables</code>: including all tables of the paper
- folder <code>scripts</code>: including all scripts for running the simulation and data analysis in this paper
  - subfolder <code>yeast</code>: including all scripts of the real-world data analysis 
- folder <code>src</code>: the source codes of our R package named by <code>cospa</code>

## Running the Project

> [!CAUTION]
> We use parallel computing to speed up simulations. The number of cores is set by the variable <code>cl.num</code> in each script, with a default value of 50. Please set the number of cores suitable for your computing platform.

### Step 1: Installing R package 
To perform all experiments, please install our R package <code>cospa</code> first:

```R
devtools::install_github("INFORMSJoC/2024.0691/src")
```

### Step 2: Running Scripts

> * All numerical results can be generated by the "pipeline.R" in the scripts folder. 
> * All results are stored in the results folder. 
> * For each script, the variable <code>path</code> corresponds to the absolute path of the folder 2024.0691.
> * For details, we introduce the usage of each scripts in the scripts folder as follows.

#### (1) Numerical simulations
> 
> (a) To generate all raw simulation results, run the following scripts from the scripts folder:
> 
> ```R
> source("scripts/sim-table.R") # generate all table results
> source("scripts/sim-time.R") # generate the time comparison result
> source("scripts/sim-vary.R") # generate the boxplot results
> ```
> 
> (b) After obtaining the raw simulation results, run the following scripts from the scripts folder to generate tables, line and box plots:
> 
> ```R
> source("scripts/get-table.R") # generate latex files summarizing error tables 
> source("scripts/get-line.R") # generate the time performance figures
> source("scripts/get-box.R") # generate the box plots
> ```

#### (2) Real data analysis
To obtain the real data result, run the following scripts from the subfolder "/yeast" in the scripts folder:

> (a) First, run the script "data_clean.R" to obtain the screened data that translate "yeast.rda" into "/data/yeast_preprocess_data.RData".
> ```R
> source("scripts/yeast/data_clean.R") # clean the data and save the screened data as "yeast_preprocess_data.RData" into the data folder
> ```
> 
> (b) After obtaining the screened data, run the script "/yeast/analysis.R" that estimate models by "yeast_preprocess_data.RData".
> ```R
> source("scripts/yeast/analysis.R") # compare the estimations of different methods, and save the result as table-real-screening.RData
> ```
> 
> (c) Finally, to obtain the summary table in our paper, run the script "/yeast/get-table-real.R" that outputs a latex file including the summary table.
> ```R
> source("scripts/yeast/get-table-real.R")
> ```
