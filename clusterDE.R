clusterDE = function(de.id.3c, outfile=TRUE, fn="cluster"){
  hc <- hclust(dist(de.id.3c), method="ward.D")
  if (outfile) png(filename=paste0(figdir,fn,".png"),width=1800,height=1800,res=300) 
  plot(hc, hang=-1, ylab=" ", las=2, xlab="Method: Euclidean distance - Ward criterion", main="Cluster dendrogram")
  if (outfile) dev.off()
}