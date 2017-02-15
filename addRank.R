
# addRank(d=aveACC)
addRank=function(d){
  rank.bs=apply(-d,2,rank)
  out.rank= paste0(rank.bs, "(",round(d,3),")")
  out = matrix(out.rank,nrow = nrow(d),ncol= ncol(d), byrow = F)
  rownames(out) = rownames(d)
  colnames(out) = colnames(d)
  return(out)
}