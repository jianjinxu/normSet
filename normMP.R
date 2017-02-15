

normMP=function(counts,n.p){
  f= NULL
  for (i in 1:n.p){
    f0= edgeR::calcNormFactors(counts[,c(i,n.p+i)], method="TMM")
  f=rbind(f, f0)
  }
  return(as.vector(f))
}

normMP.notrim=function(counts,n.p){
  f= NULL
  for (i in 1:n.p){
    f0= edgeR::calcNormFactors(counts[,c(i,n.p+i)], method="TMM",logratioTrim=0,sumTrim=0)
    f=rbind(f, f0)
  }
  return(as.vector(f))
}
norm.iter=function(counts,group,n.p, fdr=0.05){
  f=normMP(counts,n.p)
  f.0=runedgeRS(f,counts,group)
  is.de=(f.0$padj<fdr)
  return(normMP.notrim(counts[!is.de,],n.p))
}



