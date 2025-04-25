## Description 

The real data is a yeast eQTL mapping data that is archived in the R package [trigger](https://rdrr.io/bioc/trigger/man/yeast.html). 
For reproducing the result of data analysis, we first need to install the R package as follows 
```
install.packages("BiocManager")
BiocManager::install("trigger")
```
Run the following commands and then the yeast data can be loaded automatically.
```
require(trigger)
data(yeast)
```
Before the data analysis, we first need to the yeast data. It is implemented by running the script <code>data_clean.R</code> in folder <code>scripts/yeast</code>, 
which makes a screening for responses and covariates of yeast data, and then normalizes the responses and covariates. Correspondingly, the script <code>data_clean.R</code> will save the cleaned data in the <code>data</code> folder, which is named <code>yeast_preprocess_data.RData</code>. After obtaining <code>yeast_preprocess_data.RData</code>, running other scripts in the folder <code>scripts/yeast</code> can reproduce the results of data analysis (README.md in the parent folder introduces the usage of each script). 