

keepOnlySignificant <- function(matrixData, prot.names, outputFolderTemp, outputFileNameRmd, thresholdPVal){

	# call filter function
	flt.data <- pval.filterAndOrder(matrixData, prot.names, thresholdPVal)

    execLabel <- paste(
        c(format(Sys.time(), "%Y%m%d%H%M%S"), trunc( 
            runif(1) * 10000)), 
        collapse='')

    ## Write data in a file usable later when the report will be transformed 
    ## into HTML
    tempOutput <- paste(
        c(outputFolderTemp, '/flt_keep-sign_Rmd_data_', execLabel,'.txt'), 
        collapse='')
    write.table(flt.data, tempOutput, sep="\t") 
    
    ## prepare Rmd file
    cat(paste0('```{r results_corrected_',execLabel,', results="asis", echo=FALSE}'),
        paste0('flt.data <- ', 'read.table("', tempOutput, '", stringsAsFactors=FALSE)'),
        '',
        'kable(flt.data, digits=3)',
        '```',
        sep="\n", file=outputFileNameRmd, append=TRUE)
    
}

# function which filters and orders p-vals according to a given cutoff
pval.filterAndOrder <- function(matrixData, prot.names, thresholdPVal){

	# select significant on uncorrected results
	my.pval.sign <- subset(matrixData, p.values < thresholdPVal)
	my.pval.ord <- my.pval.sign[order(my.pval.sign$p.values),]
	my.pval.named <- my.pval.ord
	rownames(my.pval.named) <- featuredata$Majority.protein.IDs[as.integer(rownames(my.pval.ord))]

	# add simple fold.change (not log2)
	my.pval.named$fold.change.not.log2 <- 2^(-my.pval.named$fold.change)
	colnames(my.pval.named) <- c("p.values", "fold.change.log2", "fold.change")
	my.pval.named

}