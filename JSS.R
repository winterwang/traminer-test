library(TraMineR)

data("mvad")

mvad.seq <- seqdef(mvad, 17:86, xtstep = 6)
mvad.om <- seqdist(mvad.seq, method = "OM", indel = 1, sm = "TRATE")


library(cluster)

clusterward <- agnes(mvad.om, diss = TRUE, method = "ward")
mvad.cl4 <- cutree(clusterward, k = 4)
cl4.lab <- factor(mvad.cl4, labels = paste("Cluster", 1:4))

seqdplot(mvad.seq, group = cl4.lab)
