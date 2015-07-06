# DMRcate-devel

Installation instructions:

```R
install.packages(c("devtools","git2r"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("ks", "limma", "minfi", "DMRcatedata", "impute", "qvalue", "rtracklayer", "matrixStats"))
install_github("ttriche/DMRcate-devel")
```

Numerous patches (private to Ramsingh lab / USC) to propagate back into BioC version of DMRcate. 

Backoff, bigWig generation, collapsing over DMRs/VMRs for prediction, etc. 

See also Ozymandias (https://github.com/RamsinghLab/ozymandias) for an omnibus integrative QC/analysis package. 
