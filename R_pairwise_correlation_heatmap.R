# ****************************************************************
# USAGE: 
#Rscript scripts/R_pairwise_correlation_heatmap.R input/modelingIN/01_csf_filter/01_normalized_expFiltered/csf_all_samples_diagnosis_LS100K_age40_mirnaome_normalized.txt output/01_csf/03_correlation_heatmaps/csf_all_samples_diagnosis_LS100K_age40_mirnaome_normalized.pdf 11 2 initials p "initals,age,gender,detailed_diagnosis,library_size,diagnosis"

#Rscript scripts/R_pairwise_correlation_heatmap.R input/modelingIN/01_csf_filter/01_normalized_expFiltered/csf_all_samples_diagnosis_LS100K_age40_pirnaome_normalized.txt output/01_csf/03_correlation_heatmaps/csf_all_samples_diagnosis_LS100K_age40_pirnaome_normalized.pdf 11 2 initials p "initals,age,gender,detailed_diagnosis,library_size,diagnosis"

#Rscript scripts/R_pairwise_correlation_heatmap.R input/modelingIN/01_csf_filter/01_normalized_expFiltered/csf_AD_samples_diagnosis_LS100K_age40_pirnaome_normalized.txt output/01_csf/03_correlation_heatmaps/csf_AD_samples_diagnosis_LS100K_age40_pirnaome_normalized.png 11 2 initials p "initals,age,gender,detailed_diagnosis,library_size,diagnosis"

#Rscript scripts/R_pairwise_correlation_heatmap.R input/modelingIN/01_csf_filter/01_normalized_expFiltered/csf_AD_samples_diagnosis_LS100K_age40_mirnaome_normalized.txt output/01_csf/03_correlation_heatmaps/csf_AD_samples_diagnosis_LS100K_age40_mirnaome_normalized.png 11 2 initials p "initals,age,gender,detailed_diagnosis,library_size,diagnosis"

#Rscript scripts/R_pairwise_correlation_heatmap.R input/modelingIN/01_csf_filter/01_normalized_expFiltered/csf_AD_samples_diagnosis_LS100K_age40_pimirnaome_normalized.txt output/01_csf/03_correlation_heatmaps/csf_AD_samples_diagnosis_LS100K_age40_pimirnaome_normalized.png 11 2 initials p "initals,age,gender,detailed_diagnosis,library_size,diagnosis"



# ****************************************************************

################ load libraries ###############
suppressMessages( library(ggplot2))
suppressMessages( library(gplots))
suppressMessages( library(RColorBrewer))
suppressMessages( library(pheatmap))
suppressMessages( library(tools))        # for file path, basename_without_extension(file_path_sans_ext) and extensions(file_ext)
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("BHC"))
suppressPackageStartupMessages(library("dendsort"))
suppressPackageStartupMessages(library("coop"))
suppressPackageStartupMessages(library("effects"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("gtable"))
suppressPackageStartupMessages(library("stats"))
suppressPackageStartupMessages(library("grDevices"))
suppressPackageStartupMessages(library("graphics"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("reshape2"))

##### Parse command line arguments ############
warnings    <- warnings();
args        <- commandArgs(TRUE);
inputfile   <- args[1];
outputfile  <- args[2];
rnaStartCol <- as.numeric(args[3]);
rnaEndCol   <- as.numeric(args[4]);  # from the end
colName     <- args[5]               # header name of the column which will be used as names of rows and columns
cmethod     <- args[6]
annCols     <- as.vector(unlist(strsplit(args[7], ','))) # "initals,age,gender,detailed_diagnosis,library_size,diagnosis"
cohort      <- as.vector(unlist(strsplit(args[8], ','))) # "1 or 2 or 1,2 or 1,3"
diagnosis   <- as.vector(unlist(strsplit(args[9], ','))) # "0 or 1 or 0,1"

## Main logic
main <- function(){    
	## ############## DOWNLOAD AND LOAD THE DATA ###############

	## Get the input data
	cat("- Reading input file ...\n")
	allDataOrig <- data.frame(read.table(inputfile,header=TRUE,sep="\t", na.strings="NA"))

	# Subset the data for cohort and/or diagnosis
	# allData <- allDataOrig[allDataOrig$cohort%in%cohort & allDataOrig$diagnosis%in%diagnosis,]
	allData <- allDataOrig[allDataOrig$diagnosis%in%diagnosis,]

	# Get the number of colums
	ncols  <- length(names(allData))
	nrows  <- length(rownames(allData))
	cat("- Total samples : ", nrows,"\n- Total features: ", ncols, "\n")

	# Get the smallncRNA data
	cat("- rnaStartCol: ", rnaStartCol, "\n- rnaEndCol: ", rnaEndCol, "\n\n")
	Xrna <- allData[(rnaStartCol+1):(ncols-rnaEndCol)]

	# Get the rows and column headers
	column_header = names(Xrna) # ['hsa-let-7b-5p', 'hsa-let-7c-5p', 'hsa-let-7d-3p' ,...]
	row_header    = allData[[colName]]      # [1, ..., 72, 73, 74, 78, 80, 83, 84, 85, 86, 90, 91,...]

	# Set index as colName (patient information)
	rownames(Xrna) <- row_header

	# Get the annotation df needed to be ploted on the heatmap
	annDF <- subset(allData, select = names(allData) %in% annCols)
	annDF[is.na(annDF)] <- as.double("NA")
	rownames(annDF) <- row_header

	# Get the correlation matrix
	cat("- Computing pairwise correlation of columns, excluding NA/null values\n")
	if (cmethod == 'p'){
		corr_method <- 'pearson'
		cat ("\t- using pearson  : standard correlation coefficient\n")
	} else if (cmethod == 'k'){
		corr_method = 'kendall'
		cat("\t- using kendall  : Kendall Tau correlation coefficient\n")
	} else if (cmethod == 's'){
		corr_method = 'spearman'
		cat("\t- using spearman : Spearman rank correlation\n")
	}

	# Calculate the correlation matrix
	D            <- as.data.frame(coop::cosine(t(Xrna)))
	rownames(D)  <- row_header
	colnames(D)  <- row_header
	itemLabels   <- row_header
	percentiles  <- FindOptimalBinning(D, itemLabels, transposeData=TRUE, verbose=TRUE)
	discreteData <- DiscretiseData(D, percentiles=percentiles)
	hc1          <- bhc(t(discreteData), itemLabels, verbose=TRUE)
	bhc_order    <- order.dendrogram(hc1)
	bhc_ordered_discreteData <- D[bhc_order,bhc_order]

	# Plot the heatmap with the new cluster and order and save it to the output file
	# pngfile  <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(inputfile)),"_bhc.png",sep='');
	pdffile     <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_bhc.pdf",sep='');
	labelsfile  <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(inputfile)),"_bhc_cluster_labels.txt",sep='');
	psize=30

	# Flatten matrix to scale and then reshape
	matbhc              <- as.matrix(bhc_ordered_discreteData)
	scaled_ordered_data <- matrix(scale(c(matbhc)), nrow=nrows, byrow=TRUE)
	scaled_df           <- as.data.frame(scaled_ordered_data)
	rownames(scaled_df) <- rownames(bhc_ordered_discreteData)
	colnames(scaled_df) <- rownames(bhc_ordered_discreteData)
	colorP              <- colorRampPalette(c('dodgerblue3','deepskyblue2','lightblue','white','orange','red','darkred'))(256)
	#colorP             <- colorRampPalette(c('white','yellow','gold','orange','darkorange','red','darkred'))(256)
	res                 <- pheatmap(scaled_df , cluster_rows=F, cluster_cols=F, color = colorP, fontsize = psize/2, annotation_col=annDF, border_color = "grey50", main=file_path_sans_ext(basename(inputfile)), filename=pdffile, width=25, height=23)
	# res                 <- pheatmap(scaled_df , cluster_rows=F, cluster_cols=F, color = colorP, fontsize = psize/2, border_color = "grey50", main=file_path_sans_ext(basename(inputfile)), filename=pdffile, width=25, height=23)
	# Plot the dendrogram with the new cluster and order and save it to the output file
	# Create the png filename
	# pngfile  <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(inputfile)),"_dendrogram.png",sep='');
	# #Start PNG device driver to save output to figure.png
	# png(filename=pngfile, height=10,width=25, res=300, units="in");

	pdffile  <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_dendrogram.pdf",sep='');
	# Start PNG device driver to save output to figure.png
	pdf(pdffile, height=10,width=25);
	plot(hc1, xlab='', cex = 0.6, sub="")
	
	# Turn off device driver (to flush output to PNG file)
	dev.off()
	
	#	#Output cluster lables to a file
	#	labelsfile  <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(inputfile)),"_bhc_cluster_labels.txt",sep='');
	#	WriteClusterLabels(hc1, labelsfile, verbose=TRUE)

}

############ USER DEFINED FUNCTIONS ##########
##Function to write out the items labels for each cluster.
WriteClusterLabels <- function(dendro, outputFile="", verbose=FALSE){
  ##----------------------------------------------------------------------
  ## DEFINE SOME FUNCTIONS TO USE RECURSIVELY ON THE DENDROGRAM NODES ----
  ##----------------------------------------------------------------------
  ##for ease, we use discrete height labels here
  ##this hardwires the logEvidence threshold at zero, which is the point
  ##where merge and not-merge hypotheses are equally likely  
  WhereToCut <- function(n){
    attr(n,"height") <- 1##default
    if (!is.leaf(n)){
      attr(n,"height") <- 2   
      if (attr(n, "logEvidence")<0)
        attr(n,"height") <- 3
    }
    n
  }
  ##----------------------------------------------------------------------
  ## PROCESS THE DENDROGRAM NODES RECURSIVELY ----------------------------
  ##----------------------------------------------------------------------
  dendro <- dendrapply(dendro, WhereToCut);
  str(dendro)

  ##----------------------------------------------------------------------
  ## CUT THE DENDROGRAM AND PRINT THE LABELS IN EACH CLUSTER -------------
  ##----------------------------------------------------------------------
  cutDendro     <- cut(dendro, 2)
  nClusters     <- length(cutDendro$lower)
  nTotalLabels  <- length(labels(dendro))
  outputStrings <- rep("", nTotalLabels+nClusters)
  counter       <- 1
  print(outputStrings)
  quit()  
  for (i in 1:nClusters) {
	print(cutDendro$lower[[i]])
    ##extract the current dendrogram
    currentCluster <- cutDendro$lower[[i]]
    currentLabels  <- labels(currentCluster) 
    nLabels        <- length(currentLabels) 
    ##for each cluster, construct and store the labels
    outputStrings[counter] <- paste("---CLUSTER", i, "---")
    counter                <- counter + 1
    for (j in 1:nLabels){
      outputStrings[counter] <- currentLabels[j]
      counter                <- counter + 1
    }
  }
  ##----------------------------------------------------------------------
  ## IF REQUIRED, WRITE OUT THE CLUSTER LABELS TO A FILE -----------------
  ##----------------------------------------------------------------------
  if (outputFile!="") write.table(outputStrings, file=outputFile, quote=FALSE, row.names=FALSE)
  ##----------------------------------------------------------------------
  ## IF REQUIRED, PRINT THE CLUSTER LABELS OUT TO SCREEN -----------------
  ##----------------------------------------------------------------------
  if (verbose) for (i in 1:length(outputStrings)) print(outputStrings[i], quote=FALSE)
}

## Call the main function in the end
main()

