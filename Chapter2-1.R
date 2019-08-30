
##%######################################################%##
#                                                          #
####                   State sequence analysis          ####
#                                                          #
##%######################################################%##




library(TraMineR)

data("mvad")

mvad.labels <- c("employment", "further education", "higher education", 
                 "joblessness", "school", "training")

mvad.scode <- c("EM", "FE", "HE", "JL", "SC", "TR")


mvad.seq <- seqdef(mvad, 17:86, states = mvad.scode, 
                   labels = mvad.labels, xtstep = 6)


par(mfrow = c(2,2))

seqiplot(mvad.seq, withlegend = F, title = "Index plot (10 first sequences)", border = NA)

seqfplot(mvad.seq, with.legend = F, border = NA, main = "Sequence frequency plot")

seqdplot(mvad.seq, with.legend = F, border = NA, main = "State distribution plot")

seqlegend(mvad.seq, cex = 0.7)


par(mfrow = c(1,2))

seqHtplot(mvad.seq, main = "Entropy index")


Turbulences <- seqST(mvad.seq)

summary(Turbulences)
hist(Turbulences, col = "cyan", main = "Sequence turbulence")



submat <- seqsubm(mvad.seq, method = "TRATE")

dist.om1 <- seqdist(mvad.seq, method = "OM", indel = 1, sm = submat)



library(cluster)

clusterward1 <- agnes(dist.om1, diss = TRUE, method = "ward")
plot(clusterward1)

cl1.4 <- cutree(clusterward1, k = 4)
cl1.4fac <- factor(cl1.4, labels = paste("Type", 1:4))


seqdplot(mvad.seq, group = cl1.4fac, border = NA)


##%######################################################%##
#                                                          #
####              event sequence analysis               ####
#                                                          #
##%######################################################%##

mvad.seqe <- seqecreate(mvad.seq)


fsubseq <- seqefsub(mvad.seqe, pMinSupport = 0.05)


par(mfrow = c(1,1))


plot(fsubseq[1:15], col = "cyan")


discr <- seqecmpgroup(fsubseq, group = cl1.4fac)

plot(discr[1:6])
