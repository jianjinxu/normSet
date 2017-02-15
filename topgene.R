#' Most expressed top n sequences per sample
#'
#' Proportion of reads associated with lot the most expressed gene percentage
#'
#' @param counts \code{matrix} of counts
#' @param n number of most expressed sequences to return
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A \code{matrix} with the percentage of reads of the three most expressed sequences and a file named majSeq.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet modified by JX

topgene <- function(counts, n=3, group, geneName=geneName,col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
 
 if(length(geneName)==0) seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]}) else
   seqnames <- apply(counts, 2, function(x){o <- order(x, decreasing=TRUE)[1:n]; geneName[o]})
  
  p1 = apply(counts, 2, function(x){100*(sort(x, decreasing=TRUE)[1:n])/sum(x)})
  p <- apply(counts, 2, function(x){100*sum(sort(x, decreasing=TRUE)[1:n])/sum(x)})

  p <- round(p,digits=3)
  
  if (outfile) png(filename="figures/topngene.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
 
  
  x <- barplot(p, col=col[as.integer(group)], main=paste("Proportion of reads from top ",n," genes",sep=""),
               ylim=c(0, max(p)*1.2), las=2, ylab="Proportion of reads")
  legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  if (outfile) dev.off()
  
  return(invisible(cbind(seqnames,p1)))
}