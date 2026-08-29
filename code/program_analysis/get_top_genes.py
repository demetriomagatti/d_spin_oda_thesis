import numpy as np


def get_top_genes(adata, model, n=10):

    gene_names = np.asarray(adata.var_names)
    components = model.onmf_decomposition.components_
    top_n = n

    for program_idx in range(components.shape[0]):
        top_gene_indices = np.argsort(components[program_idx])[::-1][:top_n]
        top_genes = gene_names[top_gene_indices]
        top_weights = components[program_idx, top_gene_indices]

        print(f"\nProgram {program_idx + 1}")
        for gene, weight in zip(top_genes, top_weights):
            print(f"  {gene}: {weight:.4f}")

    return 


def get_gene_labels(adata, model, n=10):
    gene_names = adata.var["gene_name"].to_numpy()
    components = model.onmf_decomposition.components_

    top_n = n

    program_genes = {}

    for program_idx in range(components.shape[0]):
        top_gene_indices = np.argsort(components[program_idx])[::-1][:top_n]

        program_genes[program_idx + 1] = [
            (gene_names[i], components[program_idx, i])
            for i in top_gene_indices
        ]

        print(f"\nProgram {program_idx + 1}")
        print(", ".join(gene for gene, _ in program_genes[program_idx + 1]))    
    return 