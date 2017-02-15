


# biasvar= do.call(rbind, lapply(norm, biVa,control.name))
biVa= function(m, control.name){
  df=calcuNormSet(m=m,counts=counts, group=group,geneLength=geneLength)[[3]]
  tmp.0= df[rownames(df)%in%control.name,]
  a = log2(tmp.0+1)
  b= rowMeans(a)
  
  logfc=a - rowMeans(a)
  
  bias = mean(apply(a, 1, sd))
  variance = mean(apply(logfc, 1, var))
  out=c(bias, variance)
  return(out)
}