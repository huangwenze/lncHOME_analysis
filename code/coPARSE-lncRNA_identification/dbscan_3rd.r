#!/Share/home/zhangqf/usr/R-3.1.3/bin/Rscript


suppressPackageStartupMessages(library(dbscan))
options(show.error.messages = F)

Args <- commandArgs()
path <- Args[6] 
fi<- Args[7]
fo<- Args[8]

inf  <- file.path(path, fi)
ouf  <- file.path(path, fo)

a=read.table(inf, header = F, sep ="\t")
a=as.data.frame(cbind(a$V3))

N=nrow(a)
d=1
MinPts=5*(d+2)
Vd=(max(a$V1)-min(a$V1))
finner=function(t) {t/exp(t)}
fouter=integrate(finner, lower = 0, upper = Inf)
Ep=(Vd*MinPts*as.numeric(fouter[1])/(N*(pi^(d/2))))^(1/d)

res <- dbscan(a, eps = Ep,  minPts = MinPts)
p=cbind(a,res$cluster)
q=p[order(p$`res$cluster`,p$V1),]
#plot(p$V2~p$V1, col = p$`res$cluster`)

write.table(q, file=ouf, append = F, col.names = F, row.names = F, sep="\t")
