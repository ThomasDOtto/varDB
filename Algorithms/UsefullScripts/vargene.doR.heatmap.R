input = commandArgs()[8]
res = commandArgs()[9]

library(RColorBrewer)


d<-read.table(input,sep="\t", header=T)
nameRow=d[,1]
nameCol=names(d[,-1])
library(som)
toNorm=(d[,-1]);
toNorm <- toNorm[complete.cases(toNorm),]
nor <- normalize(toNorm, byrow=1)
nor<-nor[complete.cases(nor),]
dat=nor
dat=d[,-1]

type<-read.table(paste("Product.",res,".txt",sep=""),header=T,sep="\t");
head(type);
dat[dat>2]=2;

# change from true
library("gplots")
#col <- colorRampPalette(c("white","green","blue","yellow","red"))(n=100)
col <- colorRampPalette(c("red","blue","white"))(n=100)
#col <- colorRampPalette(c("white","blue"))(n=100)
pdf(paste(res,".DistMat.pdf",sep=""));
#png(paste(res,".DistMat.png",sep=""),width = 1480, height = 1480, units = "px", pointsize = 12)
Mat<-heatmap.2((as.matrix(dat)),distfun=dist, hclustfun=function(d) hclust(d, method="ward"),labCol=t(nameCol),labRow=type[,1],col=rev(col),trace="none",dendrogram="both",Colv=T,cexRow=0.3,cexCol=0.35,density.info= ('none'),main=paste(res," heatmap log "),symm = F)
dev.off();
write.table(paste("Order.",res,".txt",sep=""),nameRow[rev(Mat$rowInd)])            
           
# pch=ifelse(color[,2]!="black",16,3),cex=ifelse(color[,2]!="black",1,0.2)

