samples <- c(snakemake@input[["samples"]])
tissue <- snakemake@params[["tissue"]]

# snakemake param für cell-line (für allFreq)-> möglicherweise für ALLE cell-lines?
##################################

proteinAtlas <- read.table(snakemake@input[["proteinAtlas"]],header=T,as.is=T,sep="\t")
rownames(proteinAtlas) <- proteinAtlas$GeneID

ndata <- proteinAtlas[,-1]
ndata[ndata == 0] <- NA
# at least 3 tissues with non-zero values
ndata <- ndata[apply(ndata,1,FUN=function(x) { length(which(is.na(x))) < length(x)-2 } ),] 
ndata[is.na(ndata)] <- 0.04
logndata <- log2(ndata)
dim(logndata)

##################################

tLabels <- read.table(snakemake@input[["labels"]],header=T,as.is=T,sep="\t",quote="\"")

fftColumns <- 29:52 # 160-222
selFreq <- c("193","196","199")

library(gplots)

##################################
  pdf(snakemake@output[["allFreq"]],width=8,height=6)
  for (sample in samples)
  {
    fdata <- read.table(sample,as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]
    sample_name <- strsplit(tail(strsplit(sample, "/")[[1]],1),"-")[[1]][2]

    res <- cor(fdata[,fftColumns],logndata2,use="pairwise.complete.obs")

    matplot(as.numeric(sub("X","",names(fdata[,fftColumns]))),res,type="l",xlab="Frequency",ylab="Correlation",col="darkgrey",lwd=1,lty=1,main=sprintf("%s: Correlation of intensities across tissues",sample_name))
    lines(as.numeric(sub("X","",names(fdata[,fftColumns]))),cor(fdata[,fftColumns],logndata2[,tissue],method="pearson",use="pairwise.complete.obs"),col="black",lwd=2,type="b",pch=19,cex=0.6)
    legend("topright",tissue,fill="black")
  }
  dev.off()