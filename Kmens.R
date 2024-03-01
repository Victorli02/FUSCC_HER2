

library(NbClust)
library(ConsensusClusterPlus)


# SD_15%RNA --------------------------------------------------------------

## Select genes for clustering (sd top 15%)
filter <- round(15*length(SDs)/100)

Clustering_Genes <- names(SDs[1:filter])

## Expression dataset for clustering
Kmeans_RNA <- Tes_Expr[Clustering_Genes, ]
dim(Kmeans_RNA)



# 1.ConsensusCluster ------------------------------------------------------

rcc = ConsensusClusterPlus(as.matrix(Kmeans_RNA),maxK=10,reps=1000,pItem=0.8,pFeature=1,title="consensus",distance="euclidean",clusterAlg="km",plot="pdf", verbose = T)





# 2. Clustering -----------------------------------------------------------

rm(list=ls())
SD_num <- 15


filter <- round(SD_num*length(SDs)/100)

Clustering_Genes <- names(SDs[1:filter])

## Expression dataset for clustering
Kmeans_RNA <- Tes_Expr[Clustering_Genes, ]



set.seed(123)
km1          <- kmeans(t(Kmeans_RNA),4,iter.max = 1000,nstart = 100)   
table(km1$cluster)
Tes_Clusters <- km1$cluster

Tes_Clin$Tes_Clusters <- Tes_Clusters[row.names(Tes_Clin)]
xtabs(~ Tes_Clin$NMF + Tes_Clin$Tes_Clusters)
