##########################################################################################
# Script for Fiona
# load 'Weber_AML_sim' dataset, run diffcyt pipeline, create hierarchical clustering tree
# 
# Lukas Weber, Aug 2019
##########################################################################################


# there are 3 simulations: spike-in thresholds 5%, 1%, and 0.1%
# each simulation contains 3 conditions: healthy, CN, and CBF

# here we use data from one simulation (1% by default) and 2 conditions (CN and healthy); 
# code can be adjusted below to select the other simulations or conditions

# datasets are available in HDCytoData package from Bioconductor (version 3.10 onwards)


library(diffcyt)
library(HDCytoData)
library(SummarizedExperiment)
library(dplyr)
library(ape)


dir_output <- "../../outputs"
if (!dir.exists(dir_output)) dir.create(dir_output)



# ---------
# load data
# ---------

# note: adjust code here if you want to use one of the other simulations

# loading data from HDCytoData package

d_SE <- Weber_AML_sim_main_1pc_SE()



# -------------------
# show data structure
# -------------------

d_SE

rowData(d_SE)
colData(d_SE)

names(assays(d_SE))
dim(assay(d_SE))
head(assay(d_SE), 2)
summary(assay(d_SE))

names(metadata(d_SE))
metadata(d_SE)$experiment_info
str(metadata(d_SE)$experiment_info)



# -------------------------------------------------------------
# select data from one diseased condition (e.g. CN) and healthy
# -------------------------------------------------------------

# note: adjust code here if you want to select the other diseased condition (CBF)

# subset data
d_SE <- d_SE[rowData(d_SE)$group_id %in% c("healthy", "CN"), ]

dim(d_SE)
table(rowData(d_SE)$group_id)

# remove empty levels
rowData(d_SE) <- droplevels(rowData(d_SE))

# also manually update 'experiment_info' data frame to remove rows corresponding to removed samples
experiment_info <- metadata(d_SE)$experiment_info[1:10, ]
experiment_info <- droplevels(experiment_info)
experiment_info



# -------------------
# prepare for diffcyt
# -------------------

# split expression table into one table per sample

# (note: will update diffcyt 'prepareData' function so this can be done automatically)

# keep original data
d_SE_original <- d_SE

# check levels of 'sample_id' column in row data are in correct order
rowData(d_SE)
levels(rowData(d_SE)$sample_id)

# split data (note: using 'split.data.frame' to split into list of matrices; see '?split')
d_input <- split.data.frame(assay(d_SE), rowData(d_SE)$sample_id)

length(d_input)
names(d_input)
class(d_input)
sapply(d_input, class)
sapply(d_input, dim)
head(d_input[[1]], 2)



# -----------------------------------
# run diffcyt pipeline: initial steps
# -----------------------------------

# see also original script:
# https://github.com/lmweber/diffcyt-evaluations/blob/master/AML_sim/2_run_methods/main/run_AML_sim_diffcyt_DA_edgeR_main.R


# prepare data

# get experiment data and marker info
experiment_info  # from above

marker_info <- as.data.frame(colData(d_SE))
rownames(marker_info) <- NULL  # (note: will update diffcyt 'prepareData' so this line is not required)
marker_info

d_se <- prepareData(d_input, experiment_info, marker_info)

# check marker classes
colnames(d_se)[colData(d_se)$marker_class == "type"]
colnames(d_se)[colData(d_se)$marker_class == "state"]


# transform data (using arcsinh transform with cofactor = 5)
d_se <- transformData(d_se, cofactor = 5)


# run clustering (using FlowSOM with 20x20 grid, i.e. 400 clusters)
seed <- 1234
d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed_clustering = seed)

# check clustering
length(table(rowData(d_se)$cluster_id))  # number of clusters
nrow(rowData(d_se))                      # number of cells
sum(table(rowData(d_se)$cluster_id))
min(table(rowData(d_se)$cluster_id))     # size of smallest cluster
max(table(rowData(d_se)$cluster_id))     # size of largest cluster


# calculate cluster cell counts
d_counts <- calcCounts(d_se)

dim(d_counts)
rowData(d_counts)
length(assays(d_counts))

# calculate cluster medians
d_medians <- calcMedians(d_se)

dim(d_medians)
rowData(d_medians)
length(assays(d_medians))
names(assays(d_medians))

# calculate medians by cluster and marker
d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)

dim(d_medians_by_cluster_marker)
length(assays(d_medians_by_cluster_marker))

# calculate medians by sample and marker
d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)

dim(d_medians_by_sample_marker)
length(assays(d_medians_by_sample_marker))



# --------------------------------------------------------------------
# run diffcyt pipeline: test for differentially abundant (DA) clusters
# --------------------------------------------------------------------

# see also original script:
# https://github.com/lmweber/diffcyt-evaluations/blob/master/AML_sim/2_run_methods/main/run_AML_sim_diffcyt_DA_edgeR_main.R


# note: testing condition 'CN' vs. healthy
# to instead test 'CBF' vs. healthy, change the data subsetting above, and check contrast is still correct


# set up design matrix
# note: include fixed effects for 'patient_id'
design <- createDesignMatrix(experiment_info, cols_design = 1:2)
design

# set up contrast (entries corresponding to columns of design matrix)
contrast <- createContrast(c(0, 1, 0, 0, 0, 0))
contrast


# run tests
# note: using default filtering, which is appropriate for a 2-group comparison (since we
# have already subsetted the data)
out_DA <- testDA_edgeR(d_counts, design, contrast)


# show results
rowData(out_DA)

# top DA clusters
topTable(out_DA, format_vals = TRUE)

# number of significant DA clusters at 10% FDR
print(table(rowData(out_DA)$p_adj < 0.1))



# --------------------
# true spike-in status
# --------------------

# calculate number and proportion of true spike-in cells per cluster


# number of cells per sample (including spike-in cells)
n_cells <- table(rowData(d_SE_original)$sample_id)

# add spike-in status to rowData of 'd_se' object
stopifnot(nrow(d_se) == nrow(d_SE_original))
stopifnot(all(rowData(d_se)$sample_id == rowData(d_SE_original)$sample_id))

rowData(d_se)$spikein <- rowData(d_SE_original)$spikein
rowData(d_se)

# calculate number of spike-in cells per cluster
df_truth <- 
  rowData(d_se) %>% 
  as.data.frame %>% 
  group_by(cluster_id) %>% 
  summarize(n_cells = n(), n_spikein = sum(spikein), prop_spikein = mean(spikein)) %>% 
  as.data.frame

# define true spike-in clusters as containing >50% true spike-in cells
df_truth$truth <- df_truth$prop_spikein > 0.5

head(df_truth)
table(df_truth$truth)



# ---------------------
# combined output table
# ---------------------

# create combined output table of results

stopifnot(all(df_truth$cluster_id == rowData(out_DA)$cluster_id))
stopifnot(all(df_truth$cluster_id == rowData(d_counts)$cluster_id))
stopifnot(nrow(df_truth) == nrow(assay(d_counts)))

res_table <- data.frame(
  cluster_id = df_truth$cluster_id, 
  p_val = rowData(out_DA)$p_val, 
  p_adj = rowData(out_DA)$p_adj, 
  n_cells = df_truth$n_cells, 
  n_spikein = df_truth$n_spikein, 
  prop_spikein = df_truth$prop_spikein, 
  truth = df_truth$truth
)

res_table <- cbind(res_table, assay(d_counts))

# save output table
save(res_table, file = file.path(dir_output, "res_table.RData"))



# ----------------------
# create clustering tree
# ----------------------

# generate hierarchical clustering tree using cluster-level median expression values


# keep columns from 'cell type' markers only
meds_clustering <- assay(d_medians_by_cluster_marker)[, metadata(d_medians_by_cluster_marker)$id_type_markers]
colnames(meds_clustering)


# hierarchical clustering using 'hclust'
# note: not using 'members' argument for number of cells, since we want rare cell
# populations to be treated equally
# note: using Euclidean distance measure (i.e. distance between 16-dimensional cluster
# medians), and average agglomeration method
hclust_out <- hclust(dist(meds_clustering, method = "euclidean"), method = "average")

plot(hclust_out)


# can also display as fan plot using 'ape' package
library(ape)
plot(as.phylo(hclust_out), type = "fan")


# 'hclust' output object: 'merge' contains the merging scheme, which can be manipulated
# together with the 'd_counts' object (from above) to get per-sample cell counts for the
# internal nodes
str(hclust_out)

# 'as.phylo' converts 'hclust' output object to 'phylo' class
str(as.phylo(hclust_out))


# save clustering tree object

save(hclust_out, file = file.path(dir_output, "hclust_out.RData"))


