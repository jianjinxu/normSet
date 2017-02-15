## common bar plot
commBarPlot=function(de.id,fn="barPlot_common_DEG",outfile=T){
  tmp= vennCounts(t(de.id))
  tmp=  tmp[tmp[,dim(de.id)[1]+1]>0,]
  out =matrix(0,ncol=length(norm.0),nrow=length(norm))
  colnames(out) =norm
  for (i in 1:length(norm.0)){
    d = tmp[tmp[,i]==1,]
    v = rowSums(d[,1:dim(de.id)[1]])
    v1 =d[,dim(de.id)[1]+1]
    x=data.frame(v,v1)
    x1=aggregate(v1 ~ v, x, sum)
    x1$v
    x1$v1
    names(f) =x1$v
    out[x1$v,i] =x1$v1
  }   
  
  out1= apply(out,2,rev)
  if (outfile) png(filename=paste0(figdir,fn,".png"),width=2400,height=1800,res=300) 
  barplot(out1,col=col[1:ncol(out)],main="Barplot of common DE genes")
  if (outfile) dev.off()
  
  return(out)
  
}