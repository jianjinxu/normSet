### some discriminate methods created on 12/8/2016

# 
# lapply(lib.list, install.packages, character.only = TRUE)
# install.packages("e1071")
# install.packages("party")
# norm=c("none","tmm", "rle","mrn","tc","ps","uq","tbt","qq","ruvg")
# res= runedgeR(m,counts,group, geneLength)
# aveAcc = disAnalysis(res,counts,group,topn=10)
discAnalysis=function(res,counts,group,topn){
  
  
  res.sorted = res[order(res$padj),]
  id.selected = res.sorted[1:topn, "id"]
  
  gene.selected= counts[rownames(counts)%in%id.selected,]
  
  data.0=data.frame(t(gene.selected),group)
  # install.packages("caret")
  library(caret)
  # load the iris dataset
  # define training control
  train_control <- trainControl(method="LOOCV")
  # train the model
  acc =numeric()
  f = c("nb","knn","svmLinear2","nnet","rf")
  
  # don't use nb methods
  for (i in 2:length(f)){
    print(i)
    fit.0= train(group~., data=data.0, trControl=train_control, method=f[i])
    # summarize results
    # print(fit.0)
    acc=c(acc,fit.0$results$Accuracy[1])
  }
  
  return(mean(acc))
}





