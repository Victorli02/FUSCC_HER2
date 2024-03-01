

# ssGSEA ------------------------------------------------------------------

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) 
library(GSEABase)


HER2_ssGSEA_RNA <- gsva(as.matrix(RNA_all),
                            geneset,
                            kcdf= "Gaussian",
                            method="ssgsea",
                            min.sz = 10,
                            parallel.sz=1)






































