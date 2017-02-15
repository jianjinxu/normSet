# tmp =commonG(norm)

commonG =function(norm){
  alpha =0.05
  de.id = numeric()
  de.id.top =numeric()
  de.id.3c =numeric()
  aveAcc =numeric()
  topn= 300
  i=1
  for (i in 1:length(norm)){
    print(i)
    m=norm[i]
    res= runedgeR(m,counts,group, geneLength)
    
    ## common DE genes
    de.id.0= (res$pval<alpha)
    de.id=rbind(de.id, de.id.0)
    
    # res.sorted = res[order(res$FDR),]
    # id.selected = res.sorted[1:topn, "id"]
    # de.id.top = rbind(de.id.top, id.selected)
    
    
    ## 3 class classification
    tt= (res$pval<alpha)*sign(res$logFC)
    de.id.3c =rbind(de.id.3c ,tt)
  }
  ########## Table 2 discriminate analysis
  rownames(de.id.3c)= (norm)
  comm=matrix(nrow=length(norm),ncol=length(norm))
  # comm0=matrix(nrow=length(norm),ncol=length(norm))
  for(i in  1:length(norm)){
    for (j in  1:length(norm)){
      comm[i,j]=sum(de.id[i,] + de.id[j,] ==2)
      # comm0[i,j]=sum(de.id.top[i,]%in%de.id.top[j,])
    }
  }
  library(reshape2)
  t1=melt(comm)
  t1$ref= norm[t1$Var1]
  t1$alt= norm[t1$Var2]
  library(ggplot2)
  
  per = (round((rowSums(comm)/diag(comm) -1)/(ncol(comm)-1),4))
  dim(per) =c(length(norm),1)
  rownames(per) = norm
  colnames(per) ="Percentage of common genes"
  return(list(per,de.id,t1,de.id.3c))
  
}








## barplot

## classification of norm methods




