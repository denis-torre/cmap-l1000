#################################################################
#################################################################
################ CMap L1000 Analysis R Scripts ##################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
source('/Users/denis/Documents/Projects/scripts/Support.R')

##### 2. Other libraries #####

#######################################################
#######################################################
########## S1. PCA Analysis
#######################################################
#######################################################

#############################################
########## 1. Run PCA
#############################################

run_pca <- function(infile, outfile) {
	
	# Read data
	expressionDataframe <- read.table(infile, sep='\t', header=TRUE, row.names='gene_id')

	# Get N
	nGenes <- 5000

	# Get top genes
	topGenes <- names(sort(apply(expressionDataframe, 1, var), decreasing=TRUE))[1:nGenes]

	# Run PCA
	pca <- runPCA(expressionDataframe[topGenes,])

	# Save results
	save(pca, file=outfile)
}

#############################################
########## 2. Celltype PCA 
#############################################

run_celltype_pca <- function(infiles, outfile) {

	# Load libraries
	require(data.table)

	# Read expression data
	expressionDataframe <- as.data.frame(fread(infiles[[1]]))
	rownames(expressionDataframe) <- expressionDataframe$gene_id
	expressionDataframe$gene_id <- NULL

	# Read annotation data
	sampleAnnotationDataframe <- read.table(infiles[[2]], sep='\t', header=TRUE, row.names='sample_id')

	# Get unique cell names
	cellTypes <- as.character(unique(sampleAnnotationDataframe$cell_id))

	# Split samples by cell line
	splitSamples <- sapply(cellTypes, function(x) rownames(sampleAnnotationDataframe[sampleAnnotationDataframe$cell_id==x,]))

	# Get N
	nGenes <- 5000
	    
	# Get top variable genes
	topGenes <- sapply(cellTypes, function(x) names(sort(apply(expressionDataframe[,splitSamples[[x]]], 1, var), decreasing=TRUE))[1:nGenes], simplify=FALSE)

	# Run PCA
	pcaResults <- sapply(cellTypes, function(x) runPCA(expressionDataframe[topGenes[[x]], splitSamples[[x]]]), simplify=FALSE)

	# Save
	save(pcaResults, file=outfile)
}

#######################################################
#######################################################
########## S3. Prepare Plotly Tables
#######################################################
#######################################################

#############################################
########## 1. Add sample annotations
#############################################

process_sample_annotations <- function (infile, outfile) {

	# Read sample annotation dataframe
	sampleAnnotationDataframe <- read.table(infile, sep='\t', header=TRUE)

	# Fix annotation dataframe
	sampleAnnotationDataframe <- sampleAnnotationDataframe[,c('sample_id', 'cell_id', 'pert_iname', 'pert_dose','pert_2_iname', 'pert_2_dose')]

	# Add plate names
	sampleAnnotationDataframe$plateName <- sapply(sampleAnnotationDataframe$sample_id, function(x) strsplit(gsub('.','_',x,fixed=TRUE), '_')[[1]][4])

	# Add perturbation 1
	sampleAnnotationDataframe$perturbation_1 <- paste0(sampleAnnotationDataframe$pert_iname, ' (', round(sampleAnnotationDataframe$pert_dose, digits=1), ' um)')
	sampleAnnotationDataframe$perturbation_1 <- gsub(' (-666 um)', '', sampleAnnotationDataframe$perturbation_1, fixed=TRUE)

	# Add perturbation 2
	sampleAnnotationDataframe$perturbation_2 <- paste0(sampleAnnotationDataframe$pert_2_iname, ' (', round(sampleAnnotationDataframe$pert_2_dose, digits=1), ' um)')
	sampleAnnotationDataframe$perturbation_2 <- gsub('-666 (-666 um)', '', sampleAnnotationDataframe$perturbation_2, fixed=TRUE)

	# Merge perturbations
	sampleAnnotationDataframe$perturbationLabel <- paste(sampleAnnotationDataframe$perturbation_1, sampleAnnotationDataframe$perturbation_2, sep='<br>')

	# Write
	write.table(sampleAnnotationDataframe, file=outfile, row.names=FALSE, quote=FALSE, sep='\t')
}

#############################################
########## 2. Prepare PCA dataframe
#############################################

prepare_pca_dataframe <- function(infiles, outfile) {

	# Load data
	load(infiles[[1]])
	sampleAnnotationDataframe <- read.table(infiles[[2]], sep='\t', header=TRUE, row.names='sample_id')

	# Fix rownames
	rownames(sampleAnnotationDataframe) <- gsub(':', '.', rownames(sampleAnnotationDataframe))

	# Get PCs
	PCs <- c('PC1', 'PC2', 'PC3')

	# Get PCA dataframe
	pcaDataframe <- as.data.frame(pca$x[,PCs])

	# Add variance labels
	pcaDataframe <- rbind(varLabels=pca$varLabels[colnames(pcaDataframe)], pcaDataframe)

	# Merge with sample annotations
	mergedPcaDataframe <- merge(pcaDataframe, sampleAnnotationDataframe, by='row.names', all.x=TRUE)[,c('Row.names','PC1','PC2','PC3','cell_id','perturbationLabel')]
	colnames(mergedPcaDataframe)[1] <- 'sample_id'

	# Write
	write.table(mergedPcaDataframe, file=outfile, row.names=FALSE, quote=FALSE, sep='\t')
}

#############################################
########## 3. Prepare celltype PCA dataframe
#############################################

prepare_celltype_pca_dataframe <- function(infiles, outfile) {

	# Load data
	load(infiles[[1]])
	sampleAnnotationDataframe <- read.table(infiles[[2]], sep='\t', header=TRUE, row.names='sample_id')

	# Get PCs
	PCs <- c('PC1', 'PC2', 'PC3')

	# Get cell types
	cellTypes <- names(pcaResults)

	# Prepare dataframes
	celltypePcaDataframes <- sapply(cellTypes, function(x){ resultDataframe <- as.data.frame(pcaResults[[x]]$x[,PCs]);
	                                                        resultDataframe <- rbind(varLabels=pcaResults[[x]]$varLabels[PCs], resultDataframe);
	                                                        resultDataframe$plate_name <- sapply(rownames(resultDataframe), function(x) strsplit(gsub('.','_',x,fixed=TRUE), '_')[[1]][4]);
	                                                        resultDataframe <- merge(resultDataframe, sampleAnnotationDataframe, by='row.names', all.x=TRUE)
	                                                        colnames(resultDataframe)[1] <- 'sample_id';
	                                                        resultDataframe$cell_id <- x;
	                                                        resultDataframe <- resultDataframe[,c('sample_id', 'PC1', 'PC2', 'PC3', 'cell_id', 'plate_name', 'perturbationLabel')]
	                                                        return(resultDataframe);}, simplify=FALSE)
	                                                       
	# Concatenate
	celltypePcaDataframe <- do.call('rbind', celltypePcaDataframes) 

	# Write
	write.table(celltypePcaDataframe, file=outfile, row.names=FALSE, quote=FALSE, sep='\t')
}



#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

