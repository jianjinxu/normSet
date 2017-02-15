#' median ratio normalization method
#'
#' mrn normalization factor calculation
#'


#' @author Elie Maza et al.


#-------------------------------------------
# Median Ratio Normalization function (MRN)
#-------------------------------------------

mrnFactors=function(rawCounts,conditions) {
	rawCounts <- as.matrix(rawCounts)
	totalCounts <- colSums(rawCounts)
	normFactors <- totalCounts
	medianRatios <- rep(1,length(conditions))
	names(medianRatios) <- names(normFactors)
	if (sum(conditions==1)>1)
		meanA <- apply(rawCounts[,conditions==1]%*%diag(1/totalCounts[conditions==1]),1,mean)
	else
		meanA <- rawCounts[,conditions==1]/totalCounts[conditions==1]
	for (i in 2:max(conditions)) {
		if (sum(conditions==i)>1)
			meanB <- apply(rawCounts[,conditions==i]%*%diag(1/totalCounts[conditions==i]),1,mean)
		else
			meanB <- rawCounts[,conditions==i]/totalCounts[conditions==i]
		meanANot0 <- meanA[meanA>0&meanB>0]
		meanBNot0 <- meanB[meanA>0&meanB>0]
		ratios <- meanBNot0/meanANot0
		medianRatios[conditions==i] <- median(ratios)
		normFactors[conditions==i] <- medianRatios[conditions==i]*totalCounts[conditions==i]
	}
	medianRatios <- medianRatios/exp(mean(log(medianRatios)))
	
	# size factor used in edgeR
	f = normFactors/totalCounts
	f <- f/exp(mean(log(f)))
	# norm facor used in DEseq
	normFactors <- normFactors/exp(mean(log(normFactors)))
	## size factor used in edgeR

	return(list(medianRatios=medianRatios,normFactors=normFactors, sizef=f))
}

#---------
# The End
#---------
