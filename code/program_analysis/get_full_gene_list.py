import numpy as np
import pandas as pd


def get_full_gene_list(adata, model):
    gene_symbols = adata.var["gene_name"].to_numpy()
    components = model.onmf_decomposition.components_

    program_ranked_genes = {}

    for program_idx in range(components.shape[0]):
        order = np.argsort(components[program_idx])[::-1]

        program_ranked_genes[f"Program_{program_idx + 1}"] = pd.DataFrame({
            "gene": gene_symbols[order],
            "weight": components[program_idx, order]
        })

    for program_name, genes in program_ranked_genes.items():
        print(
            program_name,
            "n_genes =", len(genes),
            "top =", genes.iloc[0]["gene"],
            "top_weight =", genes.iloc[0]["weight"]
        )
    
    return program_ranked_genes