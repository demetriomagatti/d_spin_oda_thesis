# Repository structure

## Code

The code sub-repository is constructed with a main-file-does-everything logic, with subfolders organized in package-like style and main file importing from subfiles.

# LaTeX

LaTeX code and a compiled PDF file for a brief thesis-like document are available. 

The first part of the document gives theoretical introduction to algorithms used for data analysis. It covers basic concepts behind Orthogonal Nonnegative Matrix Factorization (oNMF) and Dimension-scalable Single-cell Perturbation Integration Network (D-SPIN).

The second part covers data analysis and results. 


# Requirements

Source data: `https://figshare.com/articles/dataset/Norman_et_al_2019_Perturb-seq/27766323?file=50533869`. <br>
Data is meant to be downloaded from the source and placed inside the gitignored `data` directory, which is placed at the hing same level of `code` at the top level of the repository.

The complete list of installed Python packages is listed inside file `code/requirements.txt` generated via `pip freeze`. 

## Hardcoded fixes
Package `dspin` ensures compatibility with Python version 3.9.18. Here I am using the more recent version 3.11.14. In order to fix some dtype mismatches, I introduced a couple of hardocoded modification in the package source files

___

<br>
<b>dspin.py</b>
<br><br>

Code block 

```Python
if issparse(gene_matrix):
    gene_matrix = np.asarray(gene_matrix.toarray()).astype(np.float64)
# Transform the original matrix by the oNMF summary components and normalize by standard deviation
self._onmf_rep_ori = onmf_summary.transform(
    gene_matrix / self.matrix_std)
```

replaced by

```Python
if issparse(gene_matrix):
    gene_matrix = np.asarray(
        gene_matrix.toarray(),
        dtype=onmf_summary.components_.dtype
    )            
gene_matrix_normalized = (
    gene_matrix / self.matrix_std
).astype(onmf_summary.components_.dtype)

self._onmf_rep_ori = onmf_summary.transform(
    gene_matrix_normalized
) 
```

___

Code line 

```Python
program_names = [f'P{ii}-' + ','.join(top_gene_list_filtered[ii]) for ii in range(num_spin)]
```

replaced by

```Python
program_names = [
    f'P{ii}-' + ','.join(
        str(gene).removesuffix('.0')
        for gene in top_gene_list_filtered[ii]
    )
    for ii in range(num_spin)
]     
```

___