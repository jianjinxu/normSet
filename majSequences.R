#' Most expressed sequences per sample
#'
#' Proportion of reads associated with the three most expressed sequences per sample, plot the most expressed gene percentage
#'
#' @param counts \code{matrix} of counts
#' @param n number of most expressed sequences to return
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A \code{matrix} with the percentage of reads of the three most expressed sequences and a file named majSeq.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet modified by JX

majSequences <- function(counts, n=1, group,col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
# select top 3 most expressed genes per sample
seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]}) 
# at least rank top3 in one sample, it's called most expressed genes
  seqnames <- unique(unlist(as.character(seqnames)))

  sum <- apply(counts,2,sum)
  ct <- counts[seqnames,]
  sum <- matrix(sum,nrow(ct),ncol(ct),byrow=TRUE)
  p <- round(100*ct/sum,digits=3)

  if (outfile) png(filename="figures/majSeq.png",width=min(3600,1800+800*ncol(ct)/10),height=1800,res=300)
    maj <- apply(p, 2, max)
    seqname <- rownames(p)[apply(p, 2, which.max)]
    x <- barplot(maj, col=col[as.integer(group)], main="Proportion of reads from most expressed genes (highest read count)",
	             ylim=c(0, max(maj)*1.2), las=2, ylab="Proportion of reads")
    legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
    for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex=0.8, srt=90, adj=0)
  if (outfile) dev.off()
  
  return(invisible(p))
}
