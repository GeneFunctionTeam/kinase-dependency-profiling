# ================================== #
# Analysis code as described in
# "Large scale profiling of kinase
# dependencies in cancer cell lines"
# jamesc@icr.ac.uk, 3rd March 2016
# ================================== #

# To run these scripts you will need to
# need to have preprocessCore, gplots 
# and mixtools installed. It is also 
# assumed that you are starting the 
# script from the folder called 'scripts'
# in the github repository at:
# https://github.com/GeneFunctionTeam/kinase-dependency-profiling
# If not, you may need to adjust the setwd()
# argument below to point to a different path

setwd("../analyses/")
source("../scripts/Intercell_analysis_functions.R")


# ---------------------- #
# Define the input files
# ---------------------- #
 

kinome_file <- "../data_sets/siRNA_Zscores/Intercell_v18_rc4_kinome_zp0_for_publication.txt"

tissues_file <- "../data_sets/siRNA_Zscores/Intercell_v18_rc4_tissues_zp0_for_publication.txt"

pathways_func_file <- "../data_sets/func_mut_calls/pathways_combined_exome_cnv_func_muts_150203.txt"

pathways_all_file <- "../data_sets/func_mut_calls/pathways_combined_exome_cnv_all_muts_150203.txt"

combmuts_func_file <- "../data_sets/func_mut_calls/combined_exome_cnv_func_muts_150225.txt"

combmuts_all_file <- "../data_sets/func_mut_calls/combined_exome_cnv_all_muts_150225.txt"

combmuts_classes_file <- "../data_sets/func_mut_calls/combined_exome_cnv_mut_classes_150225.txt"


# ------------------------ #
# define the results files
# ------------------------ #

uv_results_kinome_tissue_file="uv_results_v26_kinome_tissue_as_predictor_150202.txt"

uv_results_kinome_combmuts_file <- "univariate_results_v26_pancan_kinome_combmuts_150202.txt"

uv_results_kinome_pathways_file <- "univariate_results_v26_pancan_kinome_pathways_150202.txt"  

uv_results_kinome_combmuts_bytissue_file <- "univariate_results_v26_bytissue_kinome_combmuts_150202.txt"

uv_results_kinome_pathways_bytissue_file <- "univariate_results_v25_bytissue_kinome_pathways_150202.txt"



# ---------------- #
# Read in the data
# ---------------- #


kinome_combmuts <- read_rnai_mutations(
	rnai_file=kinome_file,
	func_muts_file=combmuts_func_file,
	all_muts_file=combmuts_all_file,
	mut_classes_file=combmuts_classes_file,
	tissues_file=tissues_file
	)

kinome_pathways <- read_rnai_mutations(
	rnai_file=kinome_file,
	func_muts_file=pathways_func_file,
	all_muts_file=pathways_all_file,
	mut_classes_file=pathways_func_file,
	tissues_file=tissues_file
	)



# ----------------------------------- #
# pie of cell line histotypes
# included in the analysis
# (Figure 1B in Campbell et al. 2016)
# ----------------------------------- #

tissues_classes <- read.table(
	file=tissues_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

tissue_counts <- apply(tissues_classes,2,sum)
tissue_counts_with_names <- tissue_counts
i <- NULL
for(i in 1:length(tissue_counts_with_names)){
	names(tissue_counts_with_names)[i] <- paste(
		legend_pretty_tissues[which(legend_actual_tissues == names(tissue_counts)[i])],
		" (",
		tissue_counts[i],
		")",
		sep=""
		)
}
tissue_counts_with_names_sort_order <- sort(tissue_counts_with_names, index.return=TRUE)

pdf(file="piechart_of_tissue_types_150701.pdf", width=11, height=8)
pie(
	sort(tissue_counts_with_names),
	cex=2.5,
	col=legend_col[names(tissue_counts)[tissue_counts_with_names_sort_order$ix]],
	init.angle=180
	)

dev.off()


# ----------------------------------- #
# how many cell lines have z ≤ -2
# in the combmuts kinome set? Plot
# ordered number of dependent cell
# lines per kinase
# (Figure 1C in Campbell et al. 2016)
# ----------------------------------- #

kinome <- read.table(
	kinome_file,
	sep="\t",
	header=TRUE,
	row.names=1,
	stringsAsFactors=FALSE
	)

kinome_combmuts_senstarg_counts <- apply(
	kinome,
	2,
	get_ltneg2_per_target # defined in the library
	)

# plot counts of dependent cell lines per kinase
pdf("number_dependent_cell_lines_per_kinase_150408.pdf", width=4, height=3)
par(mar=c(4,6,1,1))
plot(
	sort(kinome_combmuts_senstarg_counts, decreasing=TRUE),
	pch=19,
	cex=0.5,
	xlab="kinases",
	ylab="number of dependent\ncell lines",
	las=2
	)
abline(
	(0.1*nrow(kinome)),
	0,
	lty=2,
	)
abline(
	(0.5*nrow(kinome)),
	0,
	lty=2,
	)
abline(
	(0.7*nrow(kinome)),
	0,
	lty=2,
	)
abline(
	(0.9*nrow(kinome)),
	0,
	lty=2,
	)
text(
	680,
	(5+0.1*nrow(kinome)),
	"10%",
	)
text(
	680,
	(5+0.5*nrow(kinome)),
	"50%",
	)
text(
	680,
	(5+0.7*nrow(kinome)),
	"70%",
	)
text(
	680,
	(5+0.9*nrow(kinome)),
	"90%",
	)
dev.off()


# -------------------------------- #
# Test for kinase dependencies
# associated with histotype.
# Table S1D in Campbell et al. 2016
# -------------------------------- #

# joing the kinome data with
# the histotype data
kinome_tissue <- read_rnai_mutations(
	rnai_file=kinome_file,
	func_muts_file=tissues_file,
	all_muts_file=tissues_file,
	mut_classes_file=tissues_file,
	tissues_file=tissues_file
	)

# dependencies associated with histotype
uv_results_kinome_tissue <- run_univariate_tests(
	zscores=kinome_tissue$rnai,
	mutations=kinome_tissue$func_muts,
	all_variants=kinome_tissue$all_muts,
	sensitivity_thresholds=kinome_tissue$rnai_iqr_thresholds
	)

# Table S1D in Campbell et al. 2016
write.table(
	uv_results_kinome_tissue,
	file=uv_results_kinome_tissue_file,
	sep="\t",
	quote=FALSE,
	row.names=FALSE
	)

# re-load the results
uv_results_kinome_tissue <- read.table(
	file=uv_results_kinome_tissue_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)

#
# TO DO: Add the histotype results with
# mptest.p values to the ./analysis dir
#

# FDR correction of p-values
uv_results_kinome_tissue_fdr <- data.frame(
	uv_results_kinome_tissue,
	fdr=p.adjust(uv_results_kinome_tissue$wilcox.p)
	)

#
# Read in a filtered version of
# uv_results_kinome_tissue with
# only FDR ≤ 0.1 included
#

uv_results_kinome_tissue_filtered <- read.table(
	file="histotype_as_markers_FDR10percent_for_pineapples_150630.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)


histos <- levels(as.factor(uv_results_kinome_tissue_filtered$marker))

for(histo in histos){
	marker_results <- uv_results_kinome_tissue_filtered[which(
		uv_results_kinome_tissue_filtered$marker == histo
		),]

	make_cheese_pinapple_plot(
		marker_results,
		output_file=paste(
			histo,
			"_targets_cheese_pineapple_150630.pdf",
			sep=""
			),
		pval_colname="PermutationP"
		)
}

# make a legend to show
# rho values for cheese
# and pineapples
make_box_dot_plots_rho_legend()


# ------------------------------- #
# Test and visualise associations
# between osteosarcoma hitostype
# and FGFRi responses
# ------------------------------- #


intercell_drug_tissues_file <- "../data_sets/drug_responses/Intercell_FGFRi_drugscreens_tissues.txt"

intercell_drug_auc_file <- "../data_sets/drug_responses/Intercell_FGFRi_drugscreens_auc.txt"

intercell_drug_fgfr_file <- "../data_sets/drug_responses/Intercell_FGFRi_drugscreens_FGFR_status.txt"



# read in tissues for all cell lines screened
intercell_drug_tissues <- read.table(
	file=intercell_drug_tissues_file,
	header=TRUE,
	sep="\t",
	row.names=1
	)

# read in AUC for all cell lines screened
intercell_drug_auc <- read.table(
	file=intercell_drug_auc_file,
	header=TRUE,
	sep="\t",
	row.names=1
	)

# read in FGFR status for all cell lines screened
intercell_drug_fgfr <- read.table(
	file=intercell_drug_fgfr_file,
	header=TRUE,
	sep="\t",
	row.names=1
	)
	
combined_fgfr_status <- apply(
	intercell_drug_fgfr,
	1,
	max
	)


# Find the rows corresponding to bone and those not
# corresponding to bone - avoiding '_rep' rows

intercell_drug_auc_bone_rows <- which(
	intercell_drug_tissues$BONE == 1
	)
intercell_drug_auc_nonbone_rows <- which(
	intercell_drug_tissues$BREAST == 1 |
	intercell_drug_tissues$CENTRAL_NERVOUS_SYSTEM == 1 |
	intercell_drug_tissues$CERVICAL == 1 |
	intercell_drug_tissues$HEADNECK == 1 |
	intercell_drug_tissues$LUNG == 1	)

# subset the data to remove FGFR1/2 amp cell lines

intercell_drug_auc_bone_rows_no_FGFRamp <- which(
	intercell_drug_tissues$BONE == 1 &
	combined_fgfr_status == 0
	)

intercell_drug_auc_nonbone_rows_no_FGFRamp <- which(
	intercell_drug_tissues$BONE == 0 &
	combined_fgfr_status == 0
	)


# test AZ4547 association with BONE
wilcox.test(
	intercell_drug_auc$AZ4547[intercell_drug_auc_bone_rows],
	intercell_drug_auc$AZ4547[intercell_drug_auc_nonbone_rows],
	alternative="less"
	)

# test AZ4547 association with BONE with no FGFR amp (bone or non-bone)
wilcox.test(
	intercell_drug_auc$AZ4547[intercell_drug_auc_bone_rows_no_FGFRamp],
	intercell_drug_auc$AZ4547[intercell_drug_auc_nonbone_rows_no_FGFRamp],
	alternative="less"
	)


# visualise AZ4547 association with BONE
pdf(file="Intercell_drug_FGFRi_AZ4547_AUC_dep_on_BONE_151007.pdf", width=2.5, height=3)
par(bty="n", tcl=-0.2, mai=c(0.75, 0.8, 0.1, 0.1))

# The plots below were updated 26th June to re-scale
# the AUC values to 1/1000 so as to match the GDSC

# boxplot based on all data (wt and mut groups)
boxplot(
	intercell_drug_auc$AZ4547[intercell_drug_auc_nonbone_rows],
	intercell_drug_auc$AZ4547[intercell_drug_auc_bone_rows],
	pch="",
	names=c("other", "osteo"),
	ylim=c(0.5,1),
	las=1
	)
mtext("AZD4547 AUC", 2, line=3)

# plot bone points at 2
points(
	jitter(rep(2,times=length(intercell_drug_auc_bone_rows)), amount=0.33),
	intercell_drug_auc$AZ4547[intercell_drug_auc_bone_rows],
	col=legend_col["BONE"],
	pch=19
	)

# points for each tissue type
for(tissue in c("BREAST","CENTRAL_NERVOUS_SYSTEM","CERVICAL","HEADNECK","LUNG")){
	non_bone_rows_by_tissue <- which(
		intercell_drug_tissues[,tissue] == 1
		)
	if(length(non_bone_rows_by_tissue) > 0){
		# plot at 1
		points(
			jitter(rep(1,times=length(non_bone_rows_by_tissue)), amount=0.33),
			intercell_drug_auc$AZ4547[non_bone_rows_by_tissue],
			col=legend_col[tissue],
			pch=19
			)
	}
}		
dev.off()



# test PD173074 association with BONE
wilcox.test(
	intercell_drug_auc$PD173074[intercell_drug_auc_bone_rows_no_FGFRamp],
	intercell_drug_auc$PD173074[intercell_drug_auc_nonbone_rows_no_FGFRamp],
	alternative="less"
	)

# test PD173074 association with BONE
wilcox.test(
	intercell_drug_auc$PD173074[intercell_drug_auc_bone_rows],
	intercell_drug_auc$PD173074[intercell_drug_auc_nonbone_rows],
	alternative="less"
	)

# visualise PD173074 association with BONE
pdf(file="Intercell_drug_FGFRi_PD173074_AUC_dep_on_BONE_151007.pdf", width=2.5, height=3)
par(bty="n", tcl=-0.2, mai=c(0.75, 0.8, 0.1, 0.1))

# boxplot based on all data (wt and mut groups)
boxplot(
	intercell_drug_auc$PD173074[intercell_drug_auc_nonbone_rows],
	intercell_drug_auc$PD173074[intercell_drug_auc_bone_rows],
	pch="",
	names=c("other", "osteo"),
	ylim=c(0.5,1),
	las=1
	)
mtext("PD173074 AUC", 2, line=3)

# plot bone points at 2
points(
	jitter(rep(2,times=length(intercell_drug_auc_bone_rows)), amount=0.33),
	intercell_drug_auc$PD173074[intercell_drug_auc_bone_rows],
	col=legend_col["BONE"],
	pch=19
	)

# points for each tissue type
for(tissue in c("BREAST","CENTRAL_NERVOUS_SYSTEM","CERVICAL","HEADNECK","LUNG")){
	non_bone_rows_by_tissue <- which(
		intercell_drug_tissues[,tissue] == 1
		)
	if(length(non_bone_rows_by_tissue) > 0){
		# plot at 1
		points(
			jitter(rep(1,times=length(non_bone_rows_by_tissue)), amount=0.33),
			intercell_drug_auc$PD173074[non_bone_rows_by_tissue],
			col=legend_col[tissue],
			pch=19
			)
	}
}		
dev.off()


# plot the sanger osteo FGFRi data 
sanger_osteo_fgfr_data <- read.table(
	file="../data_sets/drug_responses/GDSC_FGFRi_osteo_data.txt",
	header=TRUE,
	sep="\t"
	)

pdf(file="Sanger_drug_FGFRi_IC50_dep_on_BONE_150626.pdf", width=2.5, height=3)
par(bty="n", tcl=-0.2, mai=c(0.75, 0.8, 0.1, 0.1))
boxplot(
	sanger_osteo_fgfr_data[,c(2,1)],
	pch="",
	names=c("other",'osteo'),
#	main="PD-173074 Sensitivity",
	ylab="PD-173074 log10 IC50 (µM)"
	)
stripchart(
	sanger_osteo_fgfr_data[,c(2,1)],
	vertical=TRUE,
	method="jitter",
	jitter=0.3,
	add=TRUE,
	pch=c(19,19),
	col=c(rgb(0,0,0,0.25), "yellow")
	)
dev.off()


# -------------------------------- #
# Test for kinase dependencies
# associated with driver mutations
# in the combined histotypes set.
# Table S1G in Campbell et al. 2016
# -------------------------------- #

uv_results_kinome_combmuts <- run_univariate_tests(
	zscores=kinome_combmuts$rnai,
	mutations=kinome_combmuts$func_muts,
	all_variants=kinome_combmuts$all_muts,
	sensitivity_thresholds=kinome_combmuts$rnai_iqr_thresholds
	)
write.table(
	uv_results_kinome_combmuts,
	file=uv_results_kinome_combmuts_file,
	sep="\t",
	quote=FALSE,
	row.names=FALSE
	)


# Read in a filtered version of
# uv_results_kinome_combmuts with
# only FDR ≤ 0.5 included

# TO DO: Get the file listed below and add to ./analyses

uv_results_kinome_combmuts_fdr <- read.table(
	file="../analyses/univariate_results_v26_pancan_kinome_combmuts_150202_permuted_filtered_vogelstein.txt",
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)


# Define the set of 21 genes with
# good represention (≥ 7 mutants).
# This list can be used to filter
# the complete set of tests
cgc_vogel_genes_with_n7 <- c(
	"CCND1_595_ENSG00000110092",
	"CDKN2A_1029_ENSG00000147889",
	"EGFR_1956_ENSG00000146648",
	"ERBB2_2064_ENSG00000141736",
	"GNAS_2778_ENSG00000087460",
	"KRAS_3845_ENSG00000133703",
	"SMAD4_4089_ENSG00000141646",
	"MDM2_4193_ENSG00000135679",
	"MYC_4609_ENSG00000136997",
	"NF1_4763_ENSG00000196712",
	"NOTCH2_4853_ENSG00000134250",
	"NRAS_4893_ENSG00000213281",
	"PIK3CA_5290_ENSG00000121879",
	"PTEN_5728_ENSG00000171862",
	"RB1_5925_ENSG00000139687",
	"MAP2K4_6416_ENSG00000065559",
	"SMARCA4_6597_ENSG00000127616",
	"STK11_6794_ENSG00000118046",
	"TP53_7157_ENSG00000141510",
	"ARID1A_8289_ENSG00000117713",
	"FBXW7_55294_ENSG00000109670"
	)


# --------------------------------------- #
# make cheese and pineapple plots for
# kinase dependencies associcated with
# driver gene mutations in the combined
# histotypes data set. As used for Figure
# 3B in Campbell et al. 2016
# --------------------------------------- #

markers <- levels(
	as.factor(
		uv_results_kinome_combmuts_fdr$marker
		)
	)
for(marker in markers){
	marker_results <- uv_results_kinome_combmuts_fdr[which(
		uv_results_kinome_combmuts_fdr$marker == marker
		),]
	make_cheese_pinapple_plot(
		marker_results,
		output_file=paste(
			marker,
			"_targets_cheese_pineapple_150630.pdf",
			sep=""
			)
		)
}


# --------------------------------------- #
# make boxplots for each kinase depdency
# associated with a driver gene mutation
# in the combined set of histotypes as 
# used in Figure 3 of Campbell et al 2016
# --------------------------------------- #

make_box_dot_plots_col_by_tissue(
	results=as.data.frame(
		uv_results_kinome_combmuts_fdr[which(
			uv_results_kinome_combmuts_fdr[,"FDR"] <= 0.5
			),]
		),
	zscores=kinome_combmuts$rnai,
	mutation.classes=kinome_combmuts$mut_classes,
	mutations=kinome_combmuts$func_muts,
	exclusions=kinome_combmuts$all_muts,
	tissues=kinome_combmuts$tissues,
	filename="boxplot_kinome_combmuts_allhists_FDR50_coloured_by_tissue_150701.pdf",
	tissue_pretty_names=legend_pretty_tissues,
	tissue_actual_names=legend_actual_tissues,
	tissue_cols=legend_col
	)



# ----------------------------------- #
# Test for kinase dependencies
# associated with driver mutations
# in the separate histotypes.
# Table S1J/K in Campbell et al. 2016
# ----------------------------------- #


# Driver gene mutations in separate histotypes
uv_results_kinome_combmuts_bytissue <- run_univariate_test_bytissue(
	kinome_combmuts
	)
write.table(
	uv_results_kinome_combmuts_bytissue,
	file=uv_results_kinome_combmuts_bytissue_file,
	sep="\t",
	quote=FALSE,
	row.names=FALSE
)

#
# TO DO: Add the FDR filetered data to the 
# ./analyses folder
#


# ------------------------------------- #
# Make boxplot of kinase dependencies
# associated with driver gene mutations
# within specific histotypes. As used
# in Figure 3 for Campbell et al. 2016
# ------------------------------------- #

fdr_tissues <- levels(as.factor(uv_results_kinome_combmuts_bytissue_fdr$tissue))
for(this_tissue in fdr_tissues){
	rows_to_plot <- which(
		kinome_combmuts$tissues[,this_tissue] == 1
		)
	make_box_dot_plots_col_by_tissue(
		results=as.data.frame(
			uv_results_kinome_combmuts_bytissue_fdr[which(
				uv_results_kinome_combmuts_bytissue_fdr[,"PermutationP"] <= 0.5 &
				uv_results_kinome_combmuts_bytissue_fdr[,"tissue"] == this_tissue
				),]
			),
		zscores=kinome_combmuts$rnai[rows_to_plot,],
		mutation.classes=kinome_combmuts$mut_classes[rows_to_plot,],
		mutations=kinome_combmuts$func_muts[rows_to_plot,],
		exclusions=kinome_combmuts$all_muts[rows_to_plot,],
		tissues=kinome_combmuts$tissues[rows_to_plot,],
		filename=paste(
			this_tissue, 			"_boxplots_kinome_combmuts_FDR50_coloured_by_tissue_150701.pdf", sep=""),
		tissue_pretty_names=legend_pretty_tissues,
		tissue_actual_names=legend_actual_tissues,
		tissue_cols=legend_col
		)
}



#
# TO DO: Add osteo/DYRK1A validation plots
#


# --------------------------------- #
# Test for kinase dependencies
# associated with driver pathways
# in the combined histotypes.
# Table S1N in Campbell et al. 2016
# --------------------------------- #

# pathway analysis in combined histotypes
uv_results_kinome_pathways <- run_univariate_tests(
	zscores=kinome_pathways$rnai,
	mutations=kinome_pathways$func_muts,
	all_variants=kinome_pathways$all_muts,
	sensitivity_thresholds=kinome_pathways$rnai_iqr_thresholds
	)
write.table(
	uv_results_kinome_pathways,
	file=uv_results_kinome_pathways_file,
	sep="\t",
	quote=FALSE,
	row.names=FALSE
)

uv_results_kinome_pathways <- read.table(
	file=uv_results_kinome_pathways_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)



# --------------------------------- #
# Test for kinase dependencies
# associated with driver pathways
# in the separate histotypes.
# Table S1O in Campbell et al. 2016
# --------------------------------- #

# pathway analysis in separate histotypes
uv_results_kinome_pathways_bytissue <- run_univariate_test_bytissue(kinome_pathways)
write.table(
	uv_results_kinome_pathways_bytissue,
	file=uv_results_kinome_pathways_bytissue_file,
	sep="\t",
	quote=FALSE,
	row.names=FALSE
)

uv_results_kinome_pathways_bytissue <- read.table(
	file=uv_results_kinome_pathways_bytissue_file,
	header=TRUE,
	sep="\t",
	stringsAsFactors=FALSE
	)


# ------------------------------ #
# Create multipanel heatmaps
# representing siRNA Z-score
# pathway status and specific
# gene mutation status. Some
# are for the combined histotype
# and others are for specific
# histotypes. As used in Figure
# 5 of Campbell et al. 2016
# ------------------------------ #


# Make a colour bar for the Z-scores
breaks=seq(-4, 4, by=0.2) 
breaks=append(breaks, 10)
breaks=append(breaks, -10, 0)
mycol <- colorpanel(n=length(breaks)-1,low="cyan",mid="black",high="yellow")
pdf("zscore_colour_bar_cyan_black_yellow_minus4toPlus4.pdf", width=2, height=6)
	color.bar(mycol, -4, title="z-score", nticks=9)
	color.bar(mycol, -4, title="z-score", nticks=9)
dev.off()

# define tissue rows for later
lung_rows <- which(
	kinome_combmuts$tissues[,"LUNG"] == 1
	)
bone_rows <- which(
	kinome_combmuts$tissues[,"BONE"] == 1
	)
ovary_rows <- which(
	kinome_combmuts$tissues[,"OVARY"] == 1
	)
oesophagus_rows <- which(
	kinome_combmuts$tissues[,"OESOPHAGUS"] == 1
	)


# TWF2 ~ SWI/SNF full
plot_z_and_mutations3(
	kinome=kinome_combmuts$rnai,
	mutations=kinome_combmuts$func_muts,
	threshold=-2,
	target="TWF2_ENSG00000247596",
	tree_genes=c(
		"PBRM1_55193_ENSG00000163939",
		"ARID2_196528_ENSG00000189079",
		"ARID1B_57492_ENSG00000049618",
		"SMARCA4_6597_ENSG00000127616",
		"ARID1A_8289_ENSG00000117713"
		),
	file="TWF2_z_and_Chromatin_SWI-SNF_mutations_plot_150622_for_Fig7.pdf",
	bottom_margin=23
	)

# CDK6 ~ MAPK.signaling..colm
plot_z_and_mutations3(
	kinome=kinome_combmuts$rnai,
	mutations=kinome_combmuts$func_muts,
	threshold=-2,
	target="CDK6_ENSG00000105810",
	tree_genes=c(
		"HRAS_3265_ENSG00000174775",
		"BRAF_673_ENSG00000157764",
		"NRAS_4893_ENSG00000213281",
		"KRAS_3845_ENSG00000133703"
		),
	file="CDK6_z_and_MAPK.signaling..colm_mutations_plot_150622_for_Fig7.pdf",
	bottom_margin=24
	)


# Chromatin SWI/SNF complex - UCK2 - OVARY
plot_z_and_mutations3(
	kinome=kinome_combmuts$rnai[ovary_rows,],
	mutations=kinome_combmuts$func_muts[ovary_rows,],
	threshold=-2,
	target="UCK2_ENSG00000143179",
	tree_genes=c(
#	"PBRM1_55193_ENSG00000163939", # not mutated in Ov lines
#	"ARID2_196528_ENSG00000189079",# not mutated in Ov lines
#	"SMARCA1_6594_ENSG00000102038",# not mutated in Ov lines
	"SMARCA4_6597_ENSG00000127616",
	"ARID1B_57492_ENSG00000049618",
	"ARID1A_8289_ENSG00000117713"
	),
	file="pathway_muts_z_plots/UCK2_z_and_SwiSnfCmplx_muts_in_Ovary_plot_150701.pdf",
	bottom_margin=25,
	right_margin=85
	)


# Chromatin SWI/SNF complex - DAPK1 - OESOPHAGUS
plot_z_and_mutations3(
	kinome=kinome_combmuts$rnai[oesophagus_rows,],
	mutations=kinome_combmuts$func_muts[oesophagus_rows,],
	threshold=-2,
	target="DAPK1_ENSG00000196730",
	tree_genes=c(
#	"PBRM1_55193_ENSG00000163939", # not mutated in eso lines
	"ARID2_196528_ENSG00000189079",
#	"SMARCA1_6594_ENSG00000102038",# not mutated in eso lines
	"SMARCA4_6597_ENSG00000127616",
#	"ARID1B_57492_ENSG00000049618",# not mutated in eso lines
	"ARID1A_8289_ENSG00000117713"
	),
	file="pathway_muts_z_plots/DAKP1_z_and_SwiSnfCmplx_muts_in_Oesophagus_plot_150701.pdf",
	bottom_margin=25,
	right_margin=85
	)




