####################################################################################################
# Imports
####################################################################################################

# Python libs
import pandas as pd
import scanpy as sc
import dspin.plot as dsp

# Local functions
from eda.inspect_anndata import inspect_annotated_data
from dspin_tools.networks import build_program_network


####################################################################################################
# Setup
####################################################################################################

file_path = "../data/norman_2019.h5ad"
save_path = "../results/models"


####################################################################################################
# Read and inspect data
####################################################################################################

adata = sc.read_h5ad(file_path)
inspect_annotated_data(adata)


####################################################################################################
# Build network
####################################################################################################

# Adapt column names to match `dpsin` package expectations
adata.obs["sample_id"] = (adata.obs["guide_merged"].astype(str))
adata.obs["if_control"] = adata.obs["guide_merged"].eq("ctrl")
adata.obs["batch"] = adata.obs["gemgroup"].astype(str)

# Build network
model = build_program_network(adata=adata, save_path=save_path)


####################################################################################################
# Analysis steps
####################################################################################################


############################################################
# 1. Validate model object
from program_analysis.model_validation import validate_model_object
validate_model_object(model)


############################################################
# 2. Validate oNMF
from program_analysis.onmf_validation import validate_onmf
validate_onmf(model)


############################################################
# 3. Top genes for programs
from program_analysis.get_top_genes import get_top_genes, get_gene_labels
get_top_genes(adata, model, n=10)
get_gene_labels(adata, model, n=10)


############################################################
# 4. Correlation between programs
from program_analysis.get_correlation import (
    get_correlation_matrix,
    plot_correlation_matrix,
    get_highest_correlation_duplets,
)
cm_programs = get_correlation_matrix(model)
plot_correlation_matrix(model)
get_highest_correlation_duplets(model)
get_highest_correlation_duplets(model, direction="up")


############################################################
# 5. Full gene list
from program_analysis.get_full_gene_list import get_full_gene_list
prg = get_full_gene_list(adata, model) # Program-ranked-genes


############################################################
# 6. Enrichment analysis
from program_analysis.enrich import ora_enrich_programs, summarize_ora_results
ora = ora_enrich_programs(adata, prg)
ora_summary = summarize_ora_results(ora)


############################################################
# 7. Relative responses
from program_analysis.relative_responses import (
    get_ranked_responses, 
    check_reciprocal_responses,
    get_program_activation_summary,
    chech_program_activation_selectivity,
    plot_relative_responses,
    plot_program_dendogram
)
relative_responses = pd.DataFrame(
    model.relative_responses,
    index=[f"Program_{i}" for i in range(1, model.relative_responses.shape[0] + 1)],
    columns=model.sample_list
)
top_responses = get_ranked_responses(relative_responses, ascending=False, n=5)
bottom_responses = get_ranked_responses(relative_responses, ascending=True, n=5)
reciprocal_responses = check_reciprocal_responses(relative_responses)
program_activation_summary = get_program_activation_summary(relative_responses)
program_activation_selectivity = chech_program_activation_selectivity(relative_responses)
plot_relative_responses(relative_responses)
program_correlation = relative_responses.T.corr(method="pearson")
print(program_correlation.round(2).to_string())
plot_program_dendogram(relative_responses)


############################################################
# 8. D-SPIN results

node_names = [str(i) for i in range(model.network.shape[0])]
network_graph, network_matrix = dsp.create_directed_network(
    model.network,
    node_names=node_names,
    thres_strength=0.05,
    thres_direction=0.05
)
print(type(network_graph))
print(type(network_matrix))

modules = dsp.compute_modules(
    network_graph,
    resolution=1.0,
    seed=1
)
print(modules)

spin_name_list_short = [f"P{i}" for i in range(1, 21)]
dsp.plot_network_diagram(
    network_matrix,
    modules,
    directed=True,
    weight_thres=0.25,
    spin_name_list_short=spin_name_list_short    
)