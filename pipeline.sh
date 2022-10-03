#!/bin/bash 

# analysis is meant to be run in the unzipped repository
R_version="/opt/R/4.0.3/bin/Rscript"
script_dir="./R"


# setup R virtual enviornmment for reproducibility
${R_version} "${script_dir}/setup_renv.R"

# Normalized gene expression already doone by Dr. Ricardo Zamel who has provided annotation files. Lets do the differential gene expression
${R_version} "${script_dir}/limma_moderated_t_test.R"

# generate heatmaps 
${R_version} "${script_dir}/generate_heatmaps.R"

#generate Venn Diagrams
${R_version} "${script_dir}/venn_de_genes.R"

#generate PCA plots
${R_version} "${script_dir}/make_pca_plots.R"


