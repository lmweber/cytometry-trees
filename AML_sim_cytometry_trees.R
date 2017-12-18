##########################################################################################
# Script to generate tree of similarities and count table for AML-sim data set
#
# Lukas Weber, December 2017
##########################################################################################

# Note: The AML-sim data set contains several thresholds (percentage of spiked in AML
# blast cells) and two conditions (cytogenetically normal, CN; core binding factor
# translocation, CBF); each condition is compared with healthy (H). Here we extract the
# data for one threshold (5%) and one condition (CN) vs. healthy.


library(diffcyt)
library(flowCore)
library(SummarizedExperiment)
library(dplyr)
library(magrittr)


DIR_BENCHMARK <- "../../benchmark_data/AML_sim/data/main"
DIR_OUT <- "../outputs"


# spike-in thresholds
thresholds <- c("5pc", "1pc", "0.1pc", "0.01pc")
th <- 1

# condition names
cond_names <- c("CN", "CBF")
j <- 1

# contrasts (to compare each of 'CN' and 'CBF' vs. 'healthy')
# note: include zeros for patient_IDs fixed effects
#contrasts_list <- list(CN = c(0, 1, 0, 0, 0, 0, 0), CBF = c(0, 0, 1, 0, 0, 0, 0))
contrasts_list <- list(CN = c(0, 1, 0, 0, 0, 0))




###########
# Load data
###########

# filenames
files_healthy <- list.files(file.path(DIR_BENCHMARK, "healthy"), 
                            pattern = "\\.fcs$", full.names = TRUE)
files_CN <- list.files(file.path(DIR_BENCHMARK, "CN"), 
                       pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)
files_CBF <- list.files(file.path(DIR_BENCHMARK, "CBF"), 
                        pattern = paste0("_", thresholds[th], "\\.fcs$"), full.names = TRUE)

# load data
#files_load <- c(files_healthy, files_CN, files_CBF)
files_load <- c(files_healthy, files_CN)
files_load

d_input <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs, group IDs, patient IDs
sample_IDs <- gsub("(_[0-9]+pc$)|(_0\\.[0-9]+pc$)", "", 
                   gsub("^AML_sim_", "", 
                        gsub("\\.fcs$", "", basename(files_load))))
sample_IDs

#group_IDs <- factor(gsub("_.*$", "", sample_IDs), levels = c("healthy", "CN", "CBF"))
group_IDs <- factor(gsub("_.*$", "", sample_IDs), levels = c("healthy", "CN"))
group_IDs

patient_IDs <- factor(gsub("^.*_", "", sample_IDs))
patient_IDs

# check
data.frame(sample_IDs, group_IDs, patient_IDs)

# indices of all marker columns, lineage markers, and functional markers
# (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
# Information, p. 4)
cols_markers <- 11:41
cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
cols_func <- setdiff(cols_markers, cols_lineage)


# ------------------------------------
# choose markers to use for clustering
# ------------------------------------

cols_clustering <- cols_lineage




##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

d_se <- prepareData(d_input, sample_IDs, group_IDs, 
                    cols_markers, cols_clustering, cols_func)

colnames(d_se)[cols_clustering]
colnames(d_se)[cols_func]

# transform data
d_se <- transformData(d_se, cofactor = 5)

# clustering
# (runtime: ~60 sec for 30x30 clusters)
# (note: clustering all samples together)
seed <- 1234
d_se <- generateClusters(d_se, xdim = 30, ydim = 30, seed = seed)

length(table(rowData(d_se)$cluster))  # number of clusters
nrow(rowData(d_se))                   # number of cells
sum(table(rowData(d_se)$cluster))
min(table(rowData(d_se)$cluster))     # size of smallest cluster
max(table(rowData(d_se)$cluster))     # size of largest cluster

# calculate cluster cell counts
d_counts <- calcCounts(d_se)

dim(d_counts)
rowData(d_counts)
length(assays(d_counts))

# calculate cluster medians by sample
d_medians <- calcMedians(d_se)

dim(d_medians)
rowData(d_medians)
length(assays(d_medians))
names(assays(d_medians))

# calculate cluster medians across all samples
d_medians_all <- calcMediansAll(d_se)

dim(d_medians_all)


# ----------------------------------------------
# test for differentially abundant (DA) clusters
# ----------------------------------------------

# set up design matrix
# - note: include 'patient_IDs' as fixed effects in design matrix
design <- createDesignMatrix(group_IDs, block_IDs = patient_IDs)
design

# set up contrast matrix
contrast <- createContrast(group_IDs, contrast = contrasts_list[[j]])
contrast

# run tests
res <- testDA_edgeR(d_counts, design, contrast)

# show results
rowData(res)

# sort to show top (most highly significant) clusters first
res_sorted <- rowData(res)[order(rowData(res)$FDR), ]
print(head(res_sorted, 10))
#View(as.data.frame(res_sorted))

# number of significant DA clusters
print(table(res_sorted$FDR <= 0.05))


# ---------------------------------------------
# store results at cluster level (for plotting)
# ---------------------------------------------

res_clusters <- as.data.frame(rowData(res))




######################
# True spike-in status
######################

# calculate number and proportion of true spike-in cells per cluster


# number of cells per sample (including spike-in cells)
n_cells <- sapply(d_input, nrow)

# spike-in status for each cell
is_spikein <- unlist(sapply(d_input, function(d) exprs(d)[, "spikein"]))

stopifnot(length(is_spikein) == sum(n_cells), 
          length(is_spikein) == nrow(d_se))

# add to rowData of 'd_se' object
rowData(d_se)$spikein <- is_spikein

# number of spike-in cells per cluster
df_truth <- 
  as.data.frame(rowData(d_se)) %>% 
  group_by(cluster) %>% 
  summarize(n_cells = n(), n_spikein = sum(spikein), prop_spikein = mean(spikein)) %>% 
  as.data.frame

# define 'true' clusters as containing >50% true spike-in cells
df_truth$true <- df_truth$prop_spikein > 0.5




###################
# Save output table
###################

# create combined output table of results

stopifnot(all(df_truth$cluster == res_clusters$cluster), 
          all(df_truth$cluster == rowData(d_counts)$cluster), 
          nrow(df_truth) == nrow(assay(d_counts)))

res_table <- data.frame(cluster = df_truth$cluster, 
                        p_vals = res_clusters$PValue, 
                        p_adj = res_clusters$FDR, 
                        n_cells = df_truth$n_cells, 
                        n_spikein = df_truth$n_spikein, 
                        prop_spikein = df_truth$prop_spikein, 
                        true = df_truth$true)

res_table <- cbind(res_table, assay(d_counts))

save(res_table, file = file.path(DIR_OUT, "res_table.RData"))




########################
# Create clustering tree
########################

# generate hierarchical clustering tree using cluster-level median expression values


# use clustering columns only
meds_clustering <- assay(d_medians_all)[, cols_clustering]

str(meds_clustering)
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



