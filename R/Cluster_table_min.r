library(tidyr)
#library(dplyr)
library(dbscan)

df <- read.table("Matrix_disimilarity.txt", header=TRUE)


for (x in 1:length(df$X1)){ 
  for (y in x:length(df$X1)){
    if (y>x){
      df[y ,x] <- df[x ,y]
    }
  }
}


## Cluster with the chosen parameters DBSCAN
res <- hdbscan(df ,minPts = 2,gen_hdbscan_tree = FALSE, gen_simplified_tree = FALSE)


sink('C_IDs.txt')
for (id in 1:length(res$cluster)){ cat( res$cluster[id], "\t")
}
sink()

