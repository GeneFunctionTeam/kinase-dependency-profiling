# ================================== #
# Analysis functions used in
# "Large scale profiling of kinase
# dependencies in cancer cell lines"
# jamesc@icr.ac.uk, 3rd March 2016
# ================================== #

# This file contains functions used for
# analysis. Please run the code included
# in the script: run_Intercell_analysis.R
# which should be in this folder.

require("preprocessCore")
require(gplots)
require(mixtools)

# Define colours used for plotting
# "#C9DD03" green
# "#FFD602" yellow
# "#F9A100" orange
# "#EE7EA6" pink
# "#A71930" red
# "#CCCCCC" grey
# "#726E20" olive
# "#6E273D" damson
# "#F00034" brightred
# "#ADAFAF" lightgrey
# "#003D4C" blue

icrpal <- palette(c(
	"#C9DD03",
	"#FFD602",
	"#F9A100",
	"#EE7EA6",
	"#A71930",
	"#CCCCCC",
	"#726E20",
	"#6E273D",
	"#F00034",
	"#ADAFAF",
	"#003D4C"))

legend_pretty_tissues = c(
	"Osteosarcoma",
	"Breast",
	"Lung",
	"Head & Neck",
	"Pancreas",
	"Cervical",
	"Ovary",
	"Esophagus",
	"Endometrium",
	"CNS"
	)
	
legend_actual_tissues = c(
	"BONE",
	"BREAST",
	"LUNG",
	"HEADNECK",
	"PANCREAS",
	"CERVICAL",
	"OVARY",
	"OESOPHAGUS",
	"ENDOMETRIUM",
	"CENTRAL_NERVOUS_SYSTEM"
	)

legend_col=c(
	"yellow",
	"deeppink",
	"darkgrey",
	"firebrick4",
	"purple",
	"blue",
	"cadetblue",
	"green",
	"orange",
	"darkgoldenrod4"
	)
names(legend_col) <- legend_actual_tissues


# This is used at various points where we want to
# know which values in a vector are ≤ -2...
get_ltneg2_per_target <- function(x){
	length(which(x <= -2))
}


# Function to read in data, find the intersecting
# rownames and return a list of dataframes with 
# the common cell lines.

read_rnai_mutations <- function(
	rnai_file,
	func_muts_file,
	all_muts_file,
	mut_classes_file,
	tissues_file
	){
	
	rnai <- read.table(
		file=rnai_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	rnai_qn <- t(normalize.quantiles(t(rnai)))
	rownames(rnai_qn) <- rownames(rnai)
	colnames(rnai_qn) <- colnames(rnai)
	
	func_muts <- read.table(
		file=func_muts_file,
		sep="\t",
		header=TRUE,
		row.names=1
		)

	all_muts <- read.table(
		file=all_muts_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)

	mut_classes <- read.table(
		file=mut_classes_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	tissues <- read.table(
		file=tissues_file,
		header=TRUE,
		sep="\t",
		row.names=1
		)
	
	common_celllines <- intersect(
		rownames(rnai),
		rownames(func_muts)
		)

	common_celllines <- intersect(
		common_celllines,
		rownames(tissues)
		)

	i <- NULL
	row.index <- NULL
	rnai_muts_cmn <- NULL
	rnai_qn_muts_cmn <- NULL
	func_muts_rnai_cmn <- NULL
	all_muts_rnai_cmn <- NULL
	mut_classes_rnai_cmn <- NULL
	tissues_rnai_cmn <- NULL
	for(i in seq(1:length(common_celllines))){
		# rnai subset
		row.index <- NULL
		row.index <- which(rownames(rnai) == common_celllines[i])
		rnai_muts_cmn <- rbind(
			rnai_muts_cmn,
			rnai[row.index,]
		)
		# rnai_qn subset
		row.index <- NULL
		row.index <- which(rownames(rnai_qn) == common_celllines[i])
		rnai_qn_muts_cmn <- rbind(
			rnai_qn_muts_cmn,
			rnai_qn[row.index,]
		)
		# func_muts subset
		row.index <- NULL
		row.index <- which(rownames(func_muts) == common_celllines[i])
		func_muts_rnai_cmn <- rbind(
			func_muts_rnai_cmn,
			func_muts[row.index,]
		)
		# all_muts subset
		row.index <- NULL
		row.index <- which(rownames(all_muts) == common_celllines[i])
		all_muts_rnai_cmn <- rbind(
			all_muts_rnai_cmn,
			all_muts[row.index,]
		)
		# mut_classes subset
		row.index <- NULL
		row.index <- which(rownames(mut_classes) == common_celllines[i])
		mut_classes_rnai_cmn <- rbind(
			mut_classes_rnai_cmn,
			mut_classes[row.index,]
		)
		# tissues_rnai subset
		row.index <- NULL
		row.index <- which(rownames(tissues) == common_celllines[i])
		tissues_rnai_cmn <- rbind(
			tissues_rnai_cmn,
			tissues[row.index,]
		)
	}
	rownames(rnai_muts_cmn) <- common_celllines
	rownames(func_muts_rnai_cmn) <- common_celllines
	rownames(all_muts_rnai_cmn) <- common_celllines
	rownames(mut_classes_rnai_cmn) <- common_celllines
	rownames(tissues_rnai_cmn) <- common_celllines
	
	# calculate inter-quartile range (iqr) lower fence values
	# as a threshold for sensitivity. Do for the rnai and rnai_qn
	# data sets
	
	rnai_iqr_thresholds <- NULL
	i <- NULL
	for(i in 1:ncol(rnai)){
		rnai_iqr_stats <- quantile(rnai[,i], na.rm=TRUE)
		rnai_iqr_thresholds[i] <- rnai_iqr_stats[2] - 	((rnai_iqr_stats[4] - rnai_iqr_stats[2]) * 1.5)
	}
	names(rnai_iqr_thresholds) <- colnames(rnai)
	
	rnai_qn_iqr_thresholds <- NULL
	i <- NULL
	for(i in 1:ncol(rnai_qn)){
		rnai_qn_iqr_stats <- quantile(rnai_qn[,i], na.rm=TRUE)
		rnai_qn_iqr_thresholds[i] <- rnai_qn_iqr_stats[2] - 	((rnai_qn_iqr_stats[4] - rnai_qn_iqr_stats[2]) * 1.5)
	}
	names(rnai_qn_iqr_thresholds) <- colnames(rnai_qn)

	return(
		list(
			rnai=rnai_muts_cmn,
			rnai_qn=rnai_qn_muts_cmn,
			func_muts=func_muts_rnai_cmn,
			all_muts=all_muts_rnai_cmn,
			mut_classes=mut_classes_rnai_cmn,
			rnai_qn_iqr_thresholds=rnai_qn_iqr_thresholds,
			rnai_iqr_thresholds=rnai_iqr_thresholds,
			tissues=tissues_rnai_cmn
			)
		)
}


#
# Run a set of statistical tests and
# collect descriptive stats for kinase
# dependency associations
#

run_univariate_tests <- function(
	zscores,
	mutations,
	all_variants,
	sensitivity_thresholds=NULL,
	nperms=1000000,
	alt="less"
	){
	
	zscores <- as.matrix(zscores)
	mutations <- as.matrix(mutations)
	all_variants <- as.matrix(all_variants)
	
	results <- NULL
	i <- NULL
	for(i in seq(1:length(colnames(mutations)))){
		
		# grpA is the mutants/altered group
		grpA <- which(mutations[,i] > 0)
		
		gene <- colnames(mutations)[i]
				
		# grpB includes cell lines with no 
		# reported mutations at all in a gene
		grpB <- which(all_variants[,gene] == 0)
		
		# skip if nA < 3
		if(length(grpA) < 3 | length(grpB) < 3){
			next
		}
		
		# this is used for spearman correlation...
		mut.status <- rep(NA,nrow(zscores))
		mut.status[which(mutations[,i] == 1)] <- 1
		mut.status[which(all_variants [,i] == 0)] <- 0

		j <- NULL
		for(j in seq(1:length(colnames(zscores)))){
				
			# skip if we have no viability measurements
			# in one or other group
			if(length(na.omit(zscores[grpA,j])) < 3){
				next
			}
			if(length(na.omit(zscores[grpB,j])) < 3){
				next
			}
			
			# calc the median permutation test p-value
			real.med.diff <- median(na.omit(zscores[grpA,j])) - median(na.omit(zscores[grpB,j]))
			
			# permute the groups nperms times, sampling with group sizes equal to grpA
			sample.size <- length(grpA)

			mpt.pval <- "NA"

			# run permutation tests. If the p-value
			# is small after a reasonable number of
			# permutations, increase the number of
			# permutations and repeat
			k <- NULL
			permuted.med.diffs <- NULL
			grpAB <- c(grpA,grpB) # join up the actual cell lines we used so we can sample them.

			for(k in 1:1000){
				index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
				permuted.grpA <- grpAB[index]
				permuted.grpB <- grpAB[-index]
				permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
			}
			mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 1000		
			if(mpt.pval < 0.010){
				for(k in 1001:10000){
					index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
					permuted.grpA <- grpAB[index]
					permuted.grpB <- grpAB[-index]
					permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
				}
				mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 10000
			}

			if(mpt.pval < 0.0010){
				for(k in 10001:100000){
					index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
					permuted.grpA <- grpAB[index]
					permuted.grpB <- grpAB[-index]
					permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
				}
				mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 100000
			}

			if(mpt.pval < 0.00010){
				for(k in 100001:1000000){
					index <- sample(1:length(grpAB), size=sample.size, replace=FALSE)
					permuted.grpA <- grpAB[index]
					permuted.grpB <- grpAB[-index]
					permuted.med.diffs[k] <- median(na.omit(zscores[permuted.grpA,j]))-median(na.omit(zscores[permuted.grpB,j]))
				}
				mpt.pval <- length(which(permuted.med.diffs <= real.med.diff)) / 1000000
			}

			wilcox.p <- NA
			try(
				test <- wilcox.test(
					zscores[grpA,j],
					zscores[grpB,j],
					alternative=alt
				)
			)
			wilcox.p <- test$p.value
			
			# get the Spearman r value as an effect size estimate
			spearman <- NULL
			try(
				spearman <- cor.test(zscores[,j], mut.status, method="spearman", use="complete.obs", alternative=alt)
			)
			
			this_threshold <- 0
			if(!is.null(sensitivity_thresholds)){
				this_threshold <- sensitivity_thresholds[colnames(zscores)[j]]
			}
			
			med.grpA <- median(zscores[grpA,j], na.rm=TRUE)
			med.grpB <- median(zscores[grpB,j], na.rm=TRUE)
			mad.grpA <- mad(zscores[grpA,j], na.rm=TRUE)
			mad.grpB <- mad(zscores[grpB,j], na.rm=TRUE)
			countA.sens <- length(which(zscores[grpA,j] <= this_threshold))
			countB.sens <- length(which(zscores[grpB,j] <= this_threshold))
			med.diff <- med.grpA - med.grpB
			marker <- colnames(mutations)[i]
			target <- colnames(zscores)[j]
			nA <- length(grpA)
			nB <- length(grpB)
			nMin <- min(nA, nB)
			pcnt.grpA.sens <- 100 * countA.sens / nA
			pcnt.grpB.sens <- 100 * countB.sens / nB
			min.grpA <- min(zscores[grpA,j], na.rm=TRUE)
			min.grpB <- min(zscores[grpB,j], na.rm=TRUE)
			spearman.r <- spearman$estimate
			spearman.p <- spearman$p.value
			
			# Output the result if min sample size is 2 or more
			# The condition is redundant as we skip tests with
			# fewer than 3 samples anyway
			if(nMin > 2){
				results <- rbind(
					results,
					c(
						marker,
						target,
						nA,
						nB,
						this_threshold,
						med.grpA,
						med.grpB,
						med.diff,
						countA.sens,
						countB.sens,
						pcnt.grpA.sens,
						pcnt.grpB.sens,
						min.grpA,
						min.grpB,
						spearman.r,
						spearman.p,
						wilcox.p,
						mpt.pval
					)
				)
			}
		}
	}
	
	if(is.null(nrow(results))){
		return(NULL)
	}
	
	colnames(results) <- c(
		"marker",
		"target",
		"nA",
		"nB",
		"sensitive.threshold",
		"med.grpA",
		"med.grpB",
		"med.grpA-med.grpB",
		"count.grpA.sens",
		"count.grpB.sens",
		"percent.grpA.sens",
		"percent.grpB.sens",
		"min.grpA",
		"min.grpB",
		"spearman.r",
		"spearman.p",
		"wilcox.p",
		"mptest.p"
	)
	return(results)
} # end run_univariate_tests


#
# Wrap run_univariate_tests() in
# a loop to test specific histotypes
#

run_univariate_test_bytissue <- function(x){
	tissue_types <- colnames(x$tissues)
	uv_results_bytissue <- NULL
	tissue <- NULL
	for(tissue in tissue_types){	
		cellline_count <- sum(
			x$tissues[,tissue]
			)
		if(cellline_count < 5){
			next
		}
		tissue_rows <- which(
			x$tissues[,tissue] == 1
			)
		temp_results <- NULL
		temp_results <- run_univariate_tests(
			zscores=x$rnai[tissue_rows,],
			mutations=x$func_muts[tissue_rows,],
			all_variants=x$all_muts[tissue_rows,],
			sensitivity_thresholds=x$rnai_iqr_thresholds
			)
		
		if(is.null(nrow(temp_results))){
			print(paste("Skipping ", tissue, " - no results", sep=""))
			next
		}
		
		temp_results <- cbind(
			temp_results,
			rep(tissue, times=nrow(temp_results))
			)
		
		uv_results_bytissue <- rbind(
			uv_results_bytissue,
			temp_results
			)
	}
	colnames(
		uv_results_bytissue
		)[ncol(
			uv_results_bytissue
			)
			] <- "tissue"
	return(uv_results_bytissue)
}

#
# make a colour scale bar for
# the cheese and pineapple plots
#

make_rho_legend <- function(
	min=-1,
	max=0,
	steps=10,
	output_file="rho_legend.pdf",
	alpha_multiplier=1
	){
	alphas <- seq(from=min, to=max, length=steps)
	cols <- NULL
	for(alpha in alphas){
		cols <- c(
			cols,
			rgb(1,0,0,-(alpha*alpha_multiplier))
			)
	}
	pdf(output_file, width=2)
	par(mar=c(1,6,5,1))
	image(
		t(matrix(data=c(1:10))),
		col=cols,
		xaxt="n",
		yaxt="n",
		main="rho",
		cex.main=2
		)
	axis(
	2,
	at=seq(0,1,by=0.2),
	labels=seq(min,max,by=abs(min-max)/5),
	las=2,
	cex.axis=2
	)
	dev.off()
}


#
# Make cheese and pineapple plots
# from a results table
#

make_cheese_pinapple_plot <- function(
	data,
	output_file="test.pdf",
	cex_adjust=1,
	centerx=12,
	centery=12,
	pval_colname="mptest.p"
	){
	
	if(nrow(data) == 0){
		break
	}
	
	pval_col <- which(colnames(data) == pval_colname)
	
	# see if we are in pancan mode or tissue specific
	this_tissue <- "across histoypes"
	if(!is.null(data$tissue[1])){
		this_tissue <- data$tissue[1]
	}
	if(this_tissue == "BONE"){
		this_tissue <- "OSTEOSARCOMA"
	}
	
	# get the marker we are dealing with, check we have
	# only one and clean up the name for printing
	marker <- levels(as.factor(data$marker))
	if(length(marker) != 1){
		stop("please supply only one marker in the results")
	}
	marker_name <- strsplit(marker, "_", fixed=TRUE)[[1]][1]
	
	# order the rows by target permutation test p-value
	data <- data[order(data[,pval_col], decreasing=FALSE),]
	neglogps <- -log10(data[,pval_col])
	data <- cbind(data, neglogps)
	
	# clean up the target names
	targets <- data$target
	i <- NULL
	target_names <- NULL
	for(i in 1:length(targets)){
		target_names[i] <- strsplit(
			targets[i], "_", fixed=TRUE
			)[[1]][1]
	}
	
	# find where to place the points (as a circle)
	# we calculate one more point than there are targets and then throw
	# away the extra position to leave space for the scale markers
	target.count <- nrow(data)
	center.x <- centerx
	center.y <- centery
	t <- seq(0,2*pi,length=target.count+2)
	t <- t[-1] # we skip the first rad value to leave space for the scale bar
	
	number_of_labels <- length(t)
	half_labels <- number_of_labels / 2
	if((number_of_labels %% 2) == 1){
		half_labels <- (number_of_labels - 1) / 2
	}
	label_extensions <- abs(t[1:half_labels] - 1.5) #
	
	label_extensions_for_plot <- c(
		label_extensions,
		1.55,
		label_extensions
		)
	if((number_of_labels %% 2) == 0){
		label_extensions_for_plot <- c(
			label_extensions,
			label_extensions
			)
	}
	
	coords <- matrix(NA, nrow=length(neglogps), ncol=2)
	coords_labs <- matrix(NA, nrow=length(neglogps), ncol=2)
	i <- NULL
	for(i in 1:target.count){ 
		coords[i,1] <- center.x + (sin(t[i])*(neglogps[i]))
		coords[i,2] <- center.y + (cos(t[i])*(neglogps[i]))
		coords_labs[i,1] <- center.x + (sin(t[i])* (6.5 + label_extensions_for_plot[i]))
		coords_labs[i,2] <- center.y + (cos(t[i])* (6.5 + label_extensions_for_plot[i]))
	}

	# open a PDF for output.
	# If the radius multiplier is != 1 then adjust the size
	pdf(file=output_file,7,7)
	par(oma=c(0,0,0,0), mar=c(0,0,0,0))
	plot(
		NULL,
		NULL,
		xlim=c(0,24),
		ylim=c(0,24),
		xlab="",
		ylab="",
		xaxt="n",
		yaxt="n",
		bty="n"
		)
	
	text(
		center.x,
		center.y+9.3,
		tolower(this_tissue),
		cex=1.5
		)
	
	# draw circles at positions 1:6
	symbols(
		x=rep(center.x, times=6),
		y=rep(center.y, times=6),
		circles=seq(1,6,by=1),
		lty=2,
		add=TRUE,
		inches=FALSE
		)
	
	# draw values at the 12 o'clock position
	rect(
		center.x-1,
		center.y,
		center.x+1,
		center.y+(6),
		col="white",
		border="white"
		)
	
	# draw lines from target points to lables
	# and write label
	i <- NULL

	for(i in 1:(nrow(coords))){
		text_pos=2
		text_x_shift <- 0.2
		if(i <= nrow(coords)/2){
			text_pos=4
			text_x_shift <- -0.2
		}
		lines(
			c(coords[i,1],coords_labs[i,1]),
			c(coords[i,2],coords_labs[i,2]),
			lwd=2,
			col="grey"
			)
		text(
			(coords_labs[i,1] + text_x_shift),
			coords_labs[i,2],
			target_names[i],
			cex=1*cex_adjust,
			pos=text_pos
			)
	}
	
	# blank out the positions where we will later put the
	# target red points. This lets the scale bar show through
	i <- NULL
	for(i in (nrow(coords):1)){
		
		points(
			coords[i,1],
			coords[i,2],
			pch=19,
			col="white",
			cex=2*cex_adjust
			)		
	}
	
	# mark the scale on the circles
	i <- NULL
	for(i in 2:6){
		text(
			center.x,
			center.y+(i),
			paste(
				"1e-",
				i,
				sep=""
				)
			)
	}
	
	# draw a ball the the center
	points(
		centerx,
		centery,
		pch=19,
		cex=12,
		col="white"
		)
	points(
		centerx,
		centery,
		pch=19,
		cex=12,
		col=rgb(0,0,1,0.5)
		)

	# draw points for the targets
	i <- NULL
	for(i in (nrow(coords):1)){
		
		points(
			coords[i,1],
			coords[i,2],
			pch=19,
			col=rgb(1,0,0, (-data$spearman.r[i])),
			cex=2*cex_adjust
			)
		points(
			coords[i,1],
			coords[i,2],
			col=rgb(1,0,0),# put a red border around the spoke
			cex=2*cex_adjust
			)
	}

	# add the marker name last so not obscured by targets
	text(
		centerx,
		centery,
		marker_name,
		cex=1.4
		)

	dev.off()

} # end make_cheese_pinapple_plot


#
# Make boxplots with jitted points over
# the top. Colour points by histotypes
#

make_box_dot_plots_col_by_tissue <- function(
	results,
	zscores,
	mutation.classes,
	mutations,
	exclusions,
	tissues,
	filename,
	tissue_pretty_names,
	tissue_actual_names,
	tissue_cols,
	response_type="Z-score"
	){

	pdf(file=filename, width=2, height=3)

	par(bty="n", tcl=-0.2, mai=c(0.75, 0.7, 0.1, 0.1)) # turn off boxes for plots
	
	i <- NULL
	for(i in 1:nrow(results)){
		if(is.na(results$med.grpA.med.grpB[i])){
			next
		}

		if(results$nA[i] > 2){
			marker_gene <- strsplit(results$marker[i], "_")[[1]][1]
			target_gene <- strsplit(results$target[i], "_")[[1]][1]
						
			# start by setting all cell lines to wt
			wt_mut_grps_strings <- rep(
				"wt",
				times=length(mutations[,results$marker[i]])
				)
			
			# set the non-recurrent mutations
			wt_mut_grps_strings[which(exclusions[,results$marker[i]] == 1)] <- "non-rec. mut."
			
			# set the recurrent/functional mutations
			wt_mut_grps_strings[which(mutations[,results$marker[i]] == 1)] <- "rec. mut."
			
			wt_grp_rows <- which(wt_mut_grps_strings == "wt")
			nonfunc_mut_grp_rows <- which(wt_mut_grps_strings == "non-rec. mut.")
			func_mut_grp_rows <- which(wt_mut_grps_strings == "rec. mut.")
									
		# boxplot based on all data (wt and mut groups)
			boxplot(
				zscores[wt_grp_rows,results$target[i]],
				zscores[func_mut_grp_rows,results$target[i]],
				pch="",
				names=c("wt", "mutant"),
				border="white"
				)
			
			marker_gene_xlab <- marker_gene
			if(nchar(marker_gene_xlab) > 25){
				lastbit <- sub(".{25}", "", marker_gene_xlab)
				firstbit <- sub(last_bit, "", marker_gene_xlab)
				marker_gene_xlab <- paste(
					firstbit,
					"\n",
					lastbit,
					sep="\t"
					)
				mtext(paste(marker_gene_xlab, "status"), 1, line=3)
			}else{
				mtext(paste(marker_gene_xlab, "status"), 1, line=2)
			}
			mtext(paste(target_gene, response_type), 2, line=2)

			# points for each tissue type
			j <- NULL
			tissue_colnames <- colnames(tissues)
			for(j in 1:length(tissue_colnames)){
				tissue <- tissue_colnames[j]				
				wt_rows_by_tissue <- which(
					wt_mut_grps_strings == "wt" &
					tissues[,tissue] == 1
					)
				mutant_rows_by_tissue <- which(
					wt_mut_grps_strings == "rec. mut." &
					tissues[,tissue] == 1
					)
							
				if(length(wt_rows_by_tissue) > 0){
					# plot at 1
					points(
						jitter(rep(1,times=length(wt_rows_by_tissue)), amount=0.33),
						zscores[wt_rows_by_tissue,results$target[i]],
						col=tissue_cols[tissue],
						pch=19
						)
				}
				if(length(mutant_rows_by_tissue) > 0){
					# plot at 2
					points(
						jitter(rep(2,times=length(mutant_rows_by_tissue)), amount=0.33),
						zscores[mutant_rows_by_tissue,results$target[i]],
						col=tissue_cols[tissue],
						pch=19
						)
				}
			}
				# boxplot based on all data (wt and mut groups)
			boxplot(
				zscores[wt_grp_rows,results$target[i]],
				zscores[func_mut_grp_rows,results$target[i]],
				pch="",
				names=c("wt", "mutant"),
				border=rgb(0,0,0,0.8),
				add=TRUE
				)
			abline(
				-2,0,col="red",lty=2
				)
		}
	}
	
	dev.off()

}


# function from http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
# to draw a scale bar
color.bar <- function(
	lut,
	min,
	max=-min,
	nticks=11,
	ticks=round(seq(min, max, len=nticks), 2),
	title=''){
		scale = (length(lut)-1)/(max-min)
		plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title, cex.main=2)
		axis(2, ticks, las=1, cex=2, cex.axis=2)
		for (i in 1:(length(lut)-1)) {
			y = (i-1)/scale + min
			rect(0,y,10,y+1/scale, col=lut[i], border=NA)
		}
}


plot_z_and_mutations3 <- function(
	kinome,
	mutations,
	threshold,
	target,
	tree_genes, # list of genes used in tree
	file="z_and_mutations_plot.pdf",
	bottom_margin=20,
	right_margin=6
	){
	
	# missing values (in z-scores) are a massive pain
	# for this reason, I'm joining the rnai and mutation
	# data and na.omit()ing the lot from the start...
	
	target_mutations <- cbind(
		kinome[,target],
		mutations[,tree_genes]
		)
	target_mutations <- na.omit(target_mutations)
	colnames(target_mutations)[1] <- target
	
	pathway_status <- apply(target_mutations[,-1],1,max)
	target_mutations <- cbind(pathway_status, target_mutations)
	colnames(target_mutations)[1] <- "pathway_status"
	
	# get the sort order for cell lines by z-score
	zscore_sort_order <- sort(
		target_mutations[,target],
		index.return=TRUE
		)
	
	# any cell lines with missing values for z-scores
	# above will have been removed from zscore_sort_order
	# need to be able to exclude these cell lines from the 
	# mutation data too.
	
	celllines_ordered <- rownames(target_mutations)[zscore_sort_order$ix]
	celllines_ordered_axis_labels <- sub("_.+", "", celllines_ordered, perl=TRUE)
	
	# z-score heatmap colour palette
	breaks=seq(-4, 4, by=0.2) 
	breaks=append(breaks, 100)
	breaks=append(breaks, -100, 0)
	mycol <- colorpanel(n=length(breaks)-1,low="cyan",mid="black",high="yellow")
	
	pdf(
		file=file,
		width=14,
		height=5
		)
	
	layout(matrix(c(1,2,3,3,3,3,3,3,3,3),10,1,byrow=TRUE), respect=FALSE)
	par(oma=c(0,0,0.1,0))
	
	# add z-score colour bar
	par(mar=c(0.5,9,0.5,right_margin))
	
	image(
		as.matrix(target_mutations[celllines_ordered,target]),
		xaxt="n",
		yaxt="n",
		col=mycol,
		breaks=breaks
		)
	axis(
		side=2,
		at=0.5,
		labels=paste("si", strsplit(target,"_")[[1]][1], "\nz-score", sep=""),
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)

	# add pathway status
	par(mar=c(0.5,9,0.5,right_margin))
	
	image(
		as.matrix(target_mutations[celllines_ordered,"pathway_status"]),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#303030")
		)
	axis(
		side=2,
		at=0.5,
		labels="Pathway\nstatus",
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	
	# decide bottom margin size based on number of marker genes
	#bottom_margin_size <- round(30 - (length(tree_genes_to_print) * 1.429), digits=0)
	
	# add mutation status by gene
	par(mar=c(bottom_margin,9,0.5,right_margin))
	image(
		as.matrix(
			target_mutations[
				celllines_ordered,
				tree_genes
				]
			),
		xaxt="n",
		yaxt="n",
		col=c("#FFFFFF","#A7A7A7")
		)
	# axis for cell line names
	axis(
		side=1,
		at=seq(from=0, to=1, by=1/(nrow(target_mutations)-1)),
		labels=celllines_ordered_axis_labels,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
		
	# clean up tree_gene names for printing
	tree_genes_to_print <- NULL
	i <- NULL
	for(i in 1:length(tree_genes)){
		tree_genes_to_print[i] <- strsplit(tree_genes[i],"_")[[1]][1]
	}
	
	#axis for gene names
	axis(
		side=2,
		at=seq(from=0, to=1, by=1/(length(tree_genes)-1)),
		labels=tree_genes_to_print,
		tick=FALSE,
		lwd=0,
		las=2,
		cex.axis=1.5
		)
	dev.off()
}





