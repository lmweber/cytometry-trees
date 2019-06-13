# cytometry-trees


## Contents

This repository contains a script to perform differential abundance testing and generate a hierarchical clustering tree on our AML-sim dataset.


## AML-sim dataset

The AML-sim dataset is a semi-simulated dataset generated for our publication introducing the `diffcyt` framework:

- Weber L. M. et al. (2019), *diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering*, Communications Biology, 2:183. [Available here.](https://www.ncbi.nlm.nih.gov/pubmed/31098416)

The AML-sim dataset contains 3 simulations (5%, 1%, and 0.1% spike-in thresholds for blast cells), and 2 conditions (CN and CBF) vs. healthy. For more details, see Supplementary Note 1 in our publication above.


## Script

The main script `AML_sim_clustering_trees.R` contains code to:
- load the AML-sim dataset
- run diffcyt pipeline (clustering and testing for differential abundance)
- generate a hierarchical clustering tree

The following are some possible options that can be adjusted within the script, e.g., to change the 'difficulty' of the dataset:
- select a different simulation: i.e., 5%, 1%, or 0.1% spike-in threshold for the rare population of blast cells
- select a different condition: i.e., CBF or CN (either of these can then be compared against healthy)
- use a different number of clusters
- use a different threshold for the proportion of spike-in cells required to call a cluster a 'true' spike-in cluster (default = 50%)
- use different parameter settings for generating the hierarchical tree

