#' Descriptive visualization
#' 
#' Generate the descriptive plots for data and sample visualization, note the program will find figures folder in current directory, if not found, create figure folder
#' 
#' @param counts count data frame
#' @param group factor variable for condition
#' @param col colors used for plotting
#' @param  outfile (default=T), output figures to a figures folder

#' 
#' @return a count table with rows are genes and columns are samples
#' 
#' @examples 
#' setwd(outdir)
#' desVisu(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"),outfile=F)
#' 
#' @return topgenelist the gene list and the associated proportions of most abundant genes for each sample
#' 
#' @export
#' 

desVisu=function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"),outfile=T){
    # create the figures directory if does not exist
    if (!I("figures" %in% dir())) dir.create("figures", showWarnings=FALSE)
    
    # total number of reads per sample
    barplotTotal(counts=counts, group=group, col=col,outfile=outfile)
    
    # percentage of null counts per sample
    barplotNull(counts=counts, group=group, col=col,outfile=outfile)
    
    # distribution of counts per sample
    densityPlot(counts=counts, group=group, col=col,outfile=outfile)
    
    # features which catch the most important number of reads
    majSequences <- majSequences(counts=counts, group=group, geneName=geneName,col=col,outfile=outfile)
    # top n most abundant 
    
    topgenelist = topgene(counts=counts, group=group, geneName =rownames(counts),n=5,col=col,outfile=outfile)
    # SERE and pairwise scatter plots
    cat("Matrix of SERE statistics:\n")
    print(tabSERE(counts))
    if (ncol(counts)<=10){
      pairwiseScatterPlots(counts=counts, group=group, outfile=outfile)
    } else{
      warning("No pairwise scatter-plot produced because of a too high number of samples (>30).")
    }
    
    return(topgenelist)
  }
  
  
  
  
