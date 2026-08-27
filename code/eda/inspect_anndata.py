import numpy as np
# import pandas as pd
import scanpy as sc

def inspect_annotated_data(file_path):
    """
    Inspect the structure, expression matrix, metadata, and basic quality
    metrics of an AnnData file.

    Parameters
    ----------
    file_path : str
        Path to the .h5ad file to inspect.

    Notes
    -----
    This function performs inspection only. It does not modify the AnnData
    object or apply any filtering, normalization, or other preprocessing.
    """
        
    adata = sc.read_h5ad(file_path)

    print("=== DATASET ===")
    print(f"shape: {adata.shape}")
    print(f"X dtype: {adata.X.dtype}")

    print("\n=== X SUMMARY ===")
    x = adata.X

    if hasattr(x, "toarray"):
        x_sample = x[:1000].toarray()
    else:
        x_sample = x[:1000]

    print(f"sample shape: {x_sample.shape}")
    print(f"min: {x_sample.min():.4f}")
    print(f"max: {x_sample.max():.4f}")
    print(f"mean: {x_sample.mean():.4f}")
    print(f"median: {np.median(x_sample):.4f}")
    print(f"fraction exactly zero: {(x_sample == 0).mean():.4f}")

    print("\n=== OBS COLUMNS ===")
    for column in adata.obs.columns:
        print(f"{column}: {adata.obs[column].dtype}")

    print("\n=== GUIDE IDENTITY ===")
    print(f"unique guide identities: {adata.obs['guide_identity'].nunique()}")
    print(adata.obs["guide_identity"].value_counts().head(30))

    print("\n=== GEM GROUP ===")
    print(adata.obs["gemgroup"].value_counts().sort_index())

    print("\n=== GOOD COVERAGE ===")
    print(adata.obs["good_coverage"].value_counts(dropna=False))

    print("\n=== QC SUMMARY ===")
    qc_columns = [
        "UMI_count",
        "read_count",
        "coverage",
        "number_of_cells"
    ]

    print(adata.obs[qc_columns].describe().T)

    print("\n=== GENES ===")
    print(f"number of genes: {adata.n_vars}")
    print(
        adata.var[
            ["gene_id", "gene_name", "highly_variable", "highly_variable_rank"]
        ].head(20)
    )

    print("\n=== GENE PROGRAM ===")
    print(adata.obs["gene_program"].value_counts(dropna=False).head(30))
    
    print("\n=== PERTURBATION STRUCTURE ===")

    guide_series = adata.obs["guide_identity"].astype(str)

    print("Examples:")
    print(guide_series.drop_duplicates().head(50).to_string(index=False))

    print("\nNumber of cells per guide:")
    print(guide_series.value_counts().describe())

    print("\nSmallest populations:")
    print(guide_series.value_counts().sort_values().head(30))

    print("\nLargest populations:")
    print(guide_series.value_counts().sort_values(ascending=False).head(30))    
    
    return 

