input = commandArgs()[8]
country = commandArgs()[9]
res = commandArgs()[10]

#setwd("/Users/tdo/Google_Drive/Papers/Active/1_Paper_VAR/Data_Results/tSNE/tmp")
#input="Matrix.CIDRg.1000.4bp.328.txt"
#res="test2"
d<-read.table(input,sep="\t", header=T)
dat=d[,-1]
c<-read.table(country,sep="\t", header=T)


library(RColorBrewer)
#install.packages(pkgs = 'Rtsne', dependencies = TRUE)

library(Rtsne)

tsne <- Rtsne(as.matrix(dat), dims = 2, check_duplicates = F, perplexity=30, verbose=TRUE, max_iter = 2500)

col=d[,1]
n=length(levels(unique(col)))
col_vector=palette()
#col_vector=rainbow(n);#brewer.pal(n,"rainbow") ;
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#col_vector=sample(color, n)

names(col_vector)=unique(col)


pdf(paste(res,".pdf",sep=""))
plot(tsne$Y, col=col_vector[d[,1]],pch=c[,2],main=res,xlab="", ylab="",cex=0.7)
#legend("topleft", inset=.03,  c(expression(paste("DBL", alpha, "0", sep="")),expression(paste("DBL", alpha, "1", sep="")),expression(paste("DBL", alpha, "2", sep=""))),col=unique(col_vector), pch=c(16,16,16,16,16,16,16,16,16), pt.cex=c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.4), box.lwd = 1 ,cex=0.6)
legend("topleft", inset=.03, pch=16, names(col_vector),col=col_vector, box.lwd = 1 ,cex=0.5)
legend("bottomleft", inset=.03, c("WAF","EAF","CAF","ESEA","WSEA","SAS","Reference"), pch=c(0,1,2,3,4,5,6), box.lwd = 1 ,cex=0.6)
dev.off()

