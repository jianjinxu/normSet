# randomly select control genes from simulated data

# input raw.0 has all the gene information

# simCon(raw.0, size =10)
simCon=function(raw.0){
  size=10
  # rownames(raw.0)[sample(which(raw.0$DEGenes ==0),size =size)]
  rownames(raw.0)[which(raw.0$DEGenes ==0)]
  
}


