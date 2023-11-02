#!/Share/home/zhangqf/usr/R-3.1.3/bin/Rscript

#suppressPackageStartupMessages(library(edgeR))
options(show.error.messages = F)

calpval<-function(la,lb,lout,indata){
  p = 0
  p25=sum(indata[,lb])/sum(indata[,la])
  rt = (1-p25)/p25
  for( i in 1:length(indata[,1])){
    p = 0
    r = indata[i,lb]
    x = indata[i,la]
    if(r == 0){
      indata[i,lout] = 1
    } 
    else if(r == 1){
      for(j in 1:x){
        p = p + p25*(1-p25)**(j-1)
      }
      indata[i,lout] = p
    }else{
      p = p + p25**(r)
      if(r < x){
        for(j in r+1:x){
          #p = p + (prod(r:(j-1))/prod(1:(j-r)))*p25**(r)*(1-p25)**(j-r)
          #p = p + (prod((j-r+1):(j-1))/prod(1:(r-1)))*p25**(r)*(1-p25)**(j-r)
          p = p + (10**( sum(log10((j-r+1):(j-1)))-sum(log10(1:(r-1))) + log10(p25**(r)) + log10((1-p25)**(j-r)) ))
#print( sum(log10((j-r+1):(j-1)))-sum(log10(1:(r-1))) + log10(p25**(r)) + log10((1-p25)**(j-r)) )
#print( sum(log(((j-r+1):(j-1)), rt*p25)) + j - r )
#print( sum(log(((j-r+1):(j-1)), rt*p25)) - sum(log((1:(r-1)), p25)) +j )
#          p = p + (rt ** (sum(log(((j-r+1):(j-1)), rt*p25)) + j - r)) * (p25 ** (sum(log(((j-r+1):(j-1)), rt*p25)) - sum(log((1:(r-1)), p25)) +j))
        } 
      }
      indata[i,lout] = p
    }     
  }     
  return(indata)
}  


Args <- commandArgs()
path <- Args[6] 
fi<- Args[7]
fo<- Args[8]

inf  <- file.path(path, fi)
ouf  <- file.path(path, fo)

if (FALSE){
rows=read.table(inf, header = F, sep ="\t")
t<-cbind(rows$V2, rows$V5)
g <- c(0,1)

libSizes <- as.vector(colSums(t))
##d <- DGEList(counts=t,group=g)
d <- DGEList(counts=t,group=g,lib.size=libSizes)
d <- calcNormFactors(d)
d$common.dispersion <- 1e-5
##a <- exactTest(d, pair=g)
design <- model.matrix(~g)
fit <- glmFit(d, design)
tr <- glmLRT(fit,coef=2)
a<-data.frame(tr$table)
##b<-a[(a$PValue<0.1 & a$logFC>0),]
b<-a[(a$PValue<1 & a$logFC>0),]
#b<-a
b<-cbind(rows$V1[as.numeric(rownames(b))],rows$V3[as.numeric(rownames(b))],b)
write.table(b,file=ouf, append = T, col.names = F, row.names = F, sep="\t")
}

if (TRUE){
dtest<-read.table(file=inf, header = F, sep ="\t")

dtest<-calpval(2,5,6,dtest)
dtest<-calpval(4,5,7,dtest)
for(i in 1:length(dtest[,1])){
  dtest[i,8] = dtest[i,6]*dtest[i,7]
}
b<-dtest
b<-dtest[(dtest$V8<=0.05 & dtest$V5>=3),]
#b<-cbind( dtest$V1[as.numeric(rownames(b))], dtest$V3[as.numeric(rownames(b))], dtest$V5[as.numeric(rownames(b))], dtest$V6[as.numeric(rownames(b))], dtest$V7[as.numeric(rownames(b))], dtest$V8[as.numeric(rownames(b))] )
b<-cbind( b$V1, b$V2, b$V3, b$V4, b$V5, b$V6, b$V7, b$V8 )
write.table(file=ouf, b, append = T, col.names = F, row.names = F, sep="\t")
}

