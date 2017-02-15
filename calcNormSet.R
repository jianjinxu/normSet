

# calcNormSetScal(method ="none", counts= counts, group =group)
# f.count=calcNormSetCounts(method ="qq",counts =counts, group =group,geneLength =geneLength )
 # out=calcuNormSet(m="tpm",counts=counts, group=group,geneLength=geneLength)

calcuNormSet=function(m="tmm",counts=counts, group=group,geneLength=geneLength){
  if (m%in%c("tmm", "rle","mrn","ps","uq","tbt","none")) {
    f= calcNormSetScal(method=m, counts= counts, group =group)
    s1 = colSums(counts)*f/(exp(mean(log(colSums(counts)*f))))
   
   if (m!="none") f.counts = round(t(t(counts)/s1))
    if (m=="none") f.counts = counts
    f.lib = colSums(f.counts)
    return(list(f,f.lib, f.counts))
    }
  else if(m%in%c("tc","qq", "ruvg","rpkm","tpm")){
    f =rep(1,ncol(counts))
    f.counts=calcNormSetCounts(method = m,counts =counts, group =group,geneLength =geneLength)
    f.lib =colSums(f.counts)
    return(list(f,f.lib, f.counts))
        }else stop("No such method")
}
 


# t1=calcNormSetScal(method ="tbt", counts= counts, group =group)


calcNormSetScal=function(method ="tmm", counts= counts, group =group){
  
  f=switch(method,
            tmm =edgeR::calcNormFactors(counts, method="TMM"),
            rle = edgeR::calcNormFactors(counts, method="RLE"),
            mrn = mrnFactors(counts,conditions = as.integer(group))$sizef,
          
           
            
            ps = {
                  f=PS.Est.Depth(counts)/colSums(counts)
                  f/exp(mean(log(f)))
                  },
            uq = edgeR::calcNormFactors(counts, method="upperquartile"),
          
           tbt = {tcc <- new("TCC", counts, group)
                   tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edgeR", iteration = 1, samplesize = 100)
                  tcc$norm.factors},
            none = rep(1,ncol(counts))
  )
    f        
}


# t2=calcNormSetCounts("tpm",counts =counts, group =group,geneLength =geneLength)
calcNormSetCounts= function(method,counts =counts, group =group,geneLength =geneLength ){
  d=switch(method,
           tc <- {
             totalCounts <- colSums(counts)
             round(counts*(exp(mean(log(totalCounts)))/totalCounts))},
            qq = normalize.quantiles(as.matrix(counts)),
            ruvg = {
                      set = newSeqExpressionSet(as.matrix(counts),phenoData = data.frame(group, row.names=colnames(counts)))
                      empirical= cont1(counts=counts, group=group, n.control=100)
                      normCounts(RUVg(set, empirical, k=1))},
            rpkm =rpkm(counts, geneLength),
            tpm= do.call(cbind, lapply(1:ncol(counts), tpm, geneLength =geneLength))
  )
  rownames(d)=rownames(counts)
  colnames(d)=colnames(counts)
  d=as.data.frame(round(d))
  d
  
}

tpm =function(x,geneLength){
  rate = log(x) - log(geneLength)
  denom = log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

# a function used to generate the psedo-control genes
cont1=function(counts=counts, group=group, n.control=n.control){
  cri = 0.5
  keep = rowSums(cpm(counts)>= cri) >=2
  table(keep)
  counts = as.data.frame(counts[keep,])
  y <- DGEList(counts=counts, group=group)
  y <- estimateGLMCommonDisp(y)
  y <- estimateGLMTagwiseDisp(y)
  y.test = exactTest(y)
  y.pvalues = y.test$table$PValue
  y.test$table$adjP = p.adjust(y.pvalues, method = "BH")
  
  # top = topTags(y.test,(nrow(counts)-n.control),adjust.method="BH", sort.by="PValue", p.value=1)
  o=order(y.test$table$PValue,decreasing = T)[1:n.control]
  empirical <- rownames(counts)[o] 
  return(empirical)
}


