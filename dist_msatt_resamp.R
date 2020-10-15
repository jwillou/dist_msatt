setwd("~/Documents/Projects/pupfish/dist_msatt//")
library(scales)

####setup####
data = read.table("msatt_genedist.csv", header=T, sep=",")

#does distance from closest gene influence observed heterozygosity?
mean(data$HoM, na.rm=T)
mean(data$HoS, na.rm=T) #mean across loci are very similar

#recode with 1 het per line and populaiton/locus as an indicator variable
pdata = data.frame(ho = c(data$HoM, data$HoS), dist = data$closegenedist, pop = c(rep("Malpais", nrow(data)), rep("SaltCreek", nrow(data))), locus = c(data$locus, data$locus))
pdata = pdata[complete.cases(pdata),]

####linear model with permutations - all data####
modelout = lm(ho~dist+pop, data=pdata)
summary(modelout)

#permute to estimate pvalues (since not normal)
pOUT = data.frame(rep=seq(1,1000), int=rep(NA, 1000), dist=rep(NA, 1000), pop=rep(NA, 1000))
for(i in 1:1000){ 
  pdata$newy = sample(pdata$ho, length(pdata$ho), replace=F)
  tempmodel  = lm(newy~dist+pop, data=pdata)
  pOUT$int[i]   = coef(tempmodel)[1]
  pOUT$dist[i]  = coef(tempmodel)[2]
  pOUT$pop[i]   = coef(tempmodel)[3]
}
p.dist=nrow(pOUT[abs(coef(modelout)[2])<abs(pOUT$dist),])/i
p.pop =nrow(pOUT[abs(coef(modelout)[3])<abs(pOUT$pop),])/i

#plot data nicely
pdf(file="dist_het.pdf", height = 4, width = 4)
plot(-100,-100, xlim=c(0,650), ylim=c(0,1), xlab="Kbp to nearest gene", ylab="Heterozygosity")
Mdata = pdata[pdata$pop=="Malpais",,drop=F]
Sdata = pdata[pdata$pop=="SaltCreek",,drop=F]
points(Mdata$dist, Mdata$ho, pch=21, col=alpha("dodgerblue3", 1), bg=alpha("dodgerblue3", 0.5))
points(Sdata$dist, Sdata$ho, pch=22, col=alpha("firebrick3", 1), bg=alpha("firebrick3", 0.5))
abline(modelout)
dev.off()

####linear model with permutations - less than 100 Kbp####
spdata = pdata[pdata$dist<=100,,drop=F]
modelout = lm(ho~dist+pop, data=spdata)
summary(modelout)

#permute to estimate pvalues (since not normal)
pOUT = data.frame(rep=seq(1,1000), int=rep(NA, 1000), dist=rep(NA, 1000), pop=rep(NA, 1000))
for(i in 1:1000){ 
  spdata$newy = sample(spdata$ho, length(spdata$ho), replace=F)
  tempmodel  = lm(newy~dist+pop, data=spdata)
  pOUT$int[i]   = coef(tempmodel)[1]
  pOUT$dist[i]  = coef(tempmodel)[2]
  pOUT$pop[i]   = coef(tempmodel)[3]
}
p.dist=nrow(pOUT[abs(coef(modelout)[2])<abs(pOUT$dist),])/i
p.pop =nrow(pOUT[abs(coef(modelout)[3])<abs(pOUT$pop),])/i

#plot data nicely
pdf(file="dist_het_smalldist.pdf", height = 4, width = 4)
plot(-100,-100, xlim=c(0,50), ylim=c(0,1), xlab="Kbp to nearest gene", ylab="Heterozygosity")
Mdata = spdata[spdata$pop=="Malpais",,drop=F]
Sdata = spdata[spdata$pop=="SaltCreek",,drop=F]
points(Mdata$dist, Mdata$ho, pch=21, col=alpha("dodgerblue3", 1), bg=alpha("dodgerblue3", 0.5))
points(Sdata$dist, Sdata$ho, pch=22, col=alpha("firebrick3", 1), bg=alpha("firebrick3", 0.5))
abline(modelout)
dev.off()
