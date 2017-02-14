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
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

