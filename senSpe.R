## roc comparison
## input res and real simulated information
# true.label=raw.0$DEGenes
# res =res1
#  senSpe(res1,true.label)
# binary classification DE or non-De
senSpe= function(res,true.label){
  alpha =0.05
  
  assign.label = (res$pval < alpha)
  assign.label[which(is.na(assign.label))]=0
  # crit = res$pvalue
  crit = res$pval
  roc.0 = roc(true.label, crit)
  auc= auc(roc.0)
    # plot(roc.0)
  
  (rate.sig = sum(assign.label!=0)/length(assign.label))
  tp = sum(assign.label==1 & true.label==1)
  tn = sum(assign.label==0 & true.label==0)
  p=sum(true.label)
  n=sum(true.label==0)
  sen =  tp/p
  spe = tn/n
  fp= sum(assign.label==1 & true.label==0)
  rate.fd = fp/sum(assign.label)
  acc =(tp+tn)/length(assign.label)
  f =c(auc, rate.sig,sen, spe, rate.fd, acc)
  names(f) =c("AUC","Significance rate","Sensitivity","Specificity", "false discovery rate","Accuracy")
  return(round(f,4))
  
}

senSpeN= function(res,true.label){
  alpha =0.05
  assign.label = (res$pval < alpha)
  crit = res$pval
  roc.0 = roc(true.label, crit)
  auc= auc(roc.0)
  # plot(roc.0)
  
  return(roc.0)
  
}

