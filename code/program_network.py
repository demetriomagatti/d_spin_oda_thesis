####################################################################################################
# Imports
####################################################################################################

# Python libs
import scanpy as sc

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

# 1. Validate model object
from program_analysis.model_validation import validate_model_object
validate_model_object(model)

# 2. Validate oNMF
from program_analysis.onmf_validation import validate_onmf
validate_onmf(model)

# 3. Top genes for programs
from program_analysis.get_top_genes import get_top_genes, get_gene_labels
get_top_genes(adata, model, n=10)
get_gene_labels(adata, model, n=10)
