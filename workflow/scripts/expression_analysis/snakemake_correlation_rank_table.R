samples <- snakemake@input[["samples"]]
refSample <- snakemake@params[["refSample"]]
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
#selFreq <- c("193","196","199")

library(gplots)

##################################
  # Replace BH01 with sample you are using as reference in the rank comparison
  refSample_name <- strsplit(tail(strsplit(refSample, "/")[[1]],1),"_")[[1]][2]
  fdata <- read.table(refSample,as.is=T,sep="\t",header=T,comment.char="~")
  colnames(fdata) <- sub("X","",colnames(fdata))
  rownames(fdata) <- fdata[,1]
  fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
  logndata2 <- logndata[fdata[,1],]

  selFreq <- c(colnames(fdata) > 189 & colnames(fdata) < 200)

  refCorrelation <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")

  pdf(sprintf(snakemake@output[["aveCorRank"]]),width=10,height=15)
  for (sample in samples)
  {
    fdata <- read.table(sample,as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]
    sample_name <- strsplit(tail(strsplit(sample, "/")[[1]],1),"-")[[1]][2]

    selFreq <- c(colnames(fdata) > 189 & colnames(fdata) < 200)

    res <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")
    res <- data.frame(category=tLabels$Category,description=tLabels$Type,tissue=colnames(res),correlation=as.numeric(res),rankDiff=rank(refCorrelation)-rank(res))
    par(mfrow=c(2,1))
    textplot(head(res[order(res$rankDiff,decreasing=T),],15))
    title(sprintf("By correlation rank difference: %s (vs. %s)",sample_name, refSample_name))
    textplot(head(res[order(res$correlation),],15))
    title(sprintf("By correlation: %s",sample_name))
  }
  dev.off()

