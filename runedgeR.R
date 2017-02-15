runedgeR=function(m,counts,group, geneLength){
  library(edgeR)
  library(DESeq)
  
  scale = c("none","tmm", "rle","mrn","tc","ps","uq","tbt")
  setCounts =c("qq","ruvg")
  if (m%in%(scale)){
    f=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[1]]
    y <- DGEList(counts=counts, group=group, norm.factors = f)
  }
  if (m%in%(setCounts)){
    df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
    y <- DGEList(counts=df, group=group, norm.factors = rep(1, ncol(df)))
  }
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  y.test = exactTest(y)
  y.out=topTags(y.test,n=nrow(counts),adjust.method="BH",sort.by="none")$table
  
  sum(y.out$FDR<0.05)
  
  y.out$id =rownames(y.out)
  colnames(y.out)[c(3,4)] =c("pval","padj")
  return(y.out)
}
# res= runedgeR("tmm",counts,group, geneLength)
runedgeRS = function(f,counts,group){
  y <- DGEList(counts=counts, group=group, norm.factors = f)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  y.test = exactTest(y)
  y.out=topTags(y.test,n=nrow(counts),adjust.method="BH",sort.by="none")$table
  
  sum(y.out$FDR<0.05)
  
  y.out$id =rownames(y.out)
  colnames(y.out)[c(3,4)] =c("pval","padj")
  return(y.out)
}

runedgeR1=function(m,counts,group, geneLength){
  library(edgeR)
  library(DESeq)
  
  scale = c("none","tmm", "rle","mrn","tc","ps","uq","tbt")
  setCounts =c("qq","ruvg")
  if (m%in%(scale)){
    f=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[1]]
    y <- DGEList(counts=counts, group=group, norm.factors = f)
  }
  if (m%in%(setCounts)){
    df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
    y <- DGEList(counts=df, group=group, norm.factors = rep(1, ncol(df)))
  }
  y <- estimateGLMCommonDisp(y)
  y <- estimateGLMTagwiseDisp(y)
  y.test = exactTest(y)
  
  
  #############rlm#################
  # design <- model.matrix(~group)
  # y <- estimateGLMCommonDisp(y,design, verbose = T)
  # y <- estimateGLMTagwiseDisp(y, design)
  # fit <- glmFit(y, design)
  # lrt <- glmLRT(fit)
  # # topTags(lrt)
  # 
  # y.test =lrt
  
  plotMDS(y, col=col[as.integer(group)])
  abline(h=0,v=0,lty=2,col="lightgray")
  clusterPlot(counts.trans=cpm(counts), group=group, outfile=F)
  clusterPlot(counts.trans=counts(y,normalized=T), group=group, outfile=F)
  
  return((y.test))
}

runedgeRS1 = function(counts,group, lib=lib.n){
  y <- DGEList(counts=counts, group=group, norm.factors = rep(1,ncol(counts)), lib.size=lib)  
  
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  y.test = exactTest(y)
  y.out=topTags(y.test,n=nrow(counts),adjust.method="BH",sort.by="none")$table
  
  sum(y.out$FDR<0.05)
  
  y.out$id =rownames(y.out)
  colnames(y.out)[c(3,4)] =c("pval","padj")
  return(y.out)
}
