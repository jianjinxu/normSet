runDeSeq=function(m,counts,group, geneLength){
  library(edgeR)
  library(DESeq)
  scale = c("none","tmm", "rle","mrn","tc","ps","uq","tbt")
  setCounts =c("qq","ruvg")
  if (m%in%(scale)){
    f=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[1]]
    y = newCountDataSet(counts,group)
    sizeFactors(y) =  f
  }
  if (m%in%(setCounts)){
    df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
    y <- newCountDataSet(df, group)
    sizeFactors(y) =  rep(1, ncol(df))
  }
  
  y = estimateDispersions(y,"per-condition","maximum","local")
  y.out = nbinomTest(y,"pre","post")
  # y.out = nbinomTest(y,levels(group))
  
}

# test =counts
runDEseq2.default=function(counts,group){
  library(DESeq2)
  colData = data.frame(colnames(counts), group)
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ group)
  
  dds <- DESeq(dds)
  deseq.test <- data.frame(results(dds))
  deseq.test$lfc = deseq.test$logFC
  return(deseq.test)
}

# head(deseq.test)
# sum(deseq.test$padj <0.05)


