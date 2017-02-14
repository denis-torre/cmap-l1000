#################################################################
#################################################################
################ CMap L1000 Analysis Pipeline ###################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, glob
import pandas as pd
import rpy2.robjects as robjects
import pandas.rpy.common as com

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineCmapL1000Analysis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
gctFiles = glob.glob('rawdata.dir/POL002_*x22268.gct')
gctGeneAnnotationCols = ['pr_analyte_id', 'pr_analyte_num', 'pr_bset_id', 'pr_lua_id', 'pr_pool_id', 'pr_gene_id', 'pr_gene_symbol', 'pr_gene_title', 'pr_is_bing', 'pr_is_inf', 'pr_is_lmark', 'pr_model_id']
gctSampleAnnotationRows = ['qc_slope', 'qc_f_logp', 'qc_iqr', 'bead_batch', 'bead_revision', 'bead_set', 'cell_id', 'det_mode', 'det_plate', 'det_well', 'pert_2_dose', 'pert_2_dose_unit', 'pert_2_id', 'pert_2_iname', 'pert_2_mfc_desc', 'pert_2_mfc_id', 'pert_2_type', 'pert_2_vehicle', 'pert_dose', 'pert_dose_unit', 'pert_id', 'pert_idose', 'pert_iname', 'pert_itime', 'pert_mfc_desc', 'pert_mfc_id', 'pert_time', 'pert_time_unit', 'pert_type', 'pert_vehicle', 'pool_id', 'rna_plate', 'rna_well', 'count_mean', 'count_cv', 'provenance_code', 'inf_model']

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-cmap-l1000-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Preprocess Data
#######################################################
#######################################################

#############################################
########## 1. Merge expression data
#############################################

@follows(mkdir('f1-processed_data.dir'))

@merge(gctFiles,
	   'f1-processed_data.dir/cmap_l1000-expression_data.txt')

def mergeExpressionData(infiles, outfile):

	# Read dataframes in list
	expressionDataframeList = [pd.read_table(infile, skiprows=2).set_index('id').drop(gctGeneAnnotationCols, axis=1).drop(gctSampleAnnotationRows, axis=0).reset_index() for infile in infiles]

	# Merge dataframes
	mergedExpressionDataframe = reduce(lambda left, right: pd.merge(left, right, on='id'), expressionDataframeList).set_index('id')

	# Save
	mergedExpressionDataframe.to_csv(outfile, sep='\t', index=True, index_label='gene_id')

#############################################
########## 2. Merge sample annotations
#############################################

@merge(gctFiles,
	   'f1-processed_data.dir/cmap_l1000-sample_annotations.txt')

def mergeSampleAnnotations(infiles, outfile):

	# Read dataframes in list
	sampleAnnotationDataframeList = [pd.read_table(infile, skiprows=2).set_index('id').loc[gctSampleAnnotationRows,:].drop(gctGeneAnnotationCols, axis=1).reset_index() for infile in infiles]

	# Merge dataframes
	mergedSampleAnnotationDataframe = reduce(lambda left, right: pd.merge(left, right, on='id'), sampleAnnotationDataframeList).set_index('id').T

	# Save
	mergedSampleAnnotationDataframe.to_csv(outfile, sep='\t', index=True, index_label='sample_id')

#############################################
########## 3. Merge gene annotations
#############################################

@merge(gctFiles,
	   'f1-processed_data.dir/cmap_l1000-gene_annotations.txt')

def mergeGeneAnnotations(infiles, outfile):

	# Read dataframes in list
	geneAnnotationDataframeList = [pd.read_table(infile, skiprows=2).set_index('id').loc[:,gctGeneAnnotationCols].drop(gctSampleAnnotationRows, 0).reset_index() for infile in infiles]

	# Concatenate dataframe
	concatenatedGeneAnnotationDataframe = pd.concat(geneAnnotationDataframeList).drop_duplicates().set_index('id')

	# Save
	concatenatedGeneAnnotationDataframe.to_csv(outfile, sep='\t', index=True, index_label='gene_id')

#######################################################
#######################################################
########## S2. PCA Analysis
#######################################################
#######################################################

#############################################
########## 1. Run PCA
#############################################

@follows(mkdir('f2-pca.dir'))

@transform(mergeExpressionData,
		   regex(r'.*/(.*)-expression_data.txt'),
		   r'f2-pca.dir/\1-pca.rda')

def runPca(infile, outfile):

	# Run function
	r.run_pca(infile, outfile)

#############################################
########## 2. Celltype PCA 
#############################################

@transform(mergeExpressionData,
		   regex(r'.*/(.*)-expression_data.txt'),
		   add_inputs(mergeSampleAnnotations),
		   r'f2-pca.dir/\1-celltype_pca.rda')

def runCelltypePca(infiles, outfile):

	# Run function
	r.run_celltype_pca(list(infiles), outfile)

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
