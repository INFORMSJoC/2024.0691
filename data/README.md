## Description 

The real data is a yeast eOTL mapping data that is archived in the R package [trigger](https://rdrr.io/bioc/trigger/man/yeast.html). 
For reproducing the result of data analysis, we first need to install the R package as follows 
```
install.packages("BiocManager")
BiocManager::install("trigger")
```
Run the following commands and then the <code>yeast.rds</code> is loaded automatically.
```
require(trigger)
data(yeast)
```
After the above implementations, we run the script <code>data_clean.R</code> in folder <code>scripts/yeast</code>, 
which makes a screening for variables of <code>yeast.rds</code> and normalizes the responses and covariates. Correspondingly, the cleaned data will be saved in the <code>data</code> folder, which is named <code>yeast_preprocess_data.RData</code>. All results of data analysis are based on the cleaned data, and corresponding scripts are attached in <code>scripts/yeast</code> folder. 