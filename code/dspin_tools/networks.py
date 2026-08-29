# import scanpy as sc
from dspin.dspin import DSPIN


def build_gene_network(
    adata,
    save_path,
    sample_id_key="sample_id",
    if_control_key="if_control",
    batch_key="batch",
    num_spin='auto',
    method="pseudo_likelihood",
    directed=True,
    run_with_matlab=False,
    params=None,
    discretize_params=None
):
    """
    Load a gene-expression dataset and perform D-SPIN network inference
    followed by calculation of responses relative to control samples.

    Parameters
    ----------
    adata : anndata.AnnData
    Annotated data matrix containing the single-cell expression data,
    typically loaded from an `.h5ad` file using `scanpy.read_h5ad()`.

    save_path : str
        Directory in which D-SPIN stores intermediate and final results
        generated during network inference.

    sample_id_key : str, default="sample_id"
        Name of the column in ``adata.obs`` containing the identifier of
        each sample or perturbation. D-SPIN uses this column to determine
        which cells belong to the same sample when constructing the
        perturbation-level representation.

    if_control_key : str, default="if_control"
        Name of the column in ``adata.obs`` indicating whether a cell or
        sample belongs to a control condition. This information is used by
        ``response_relative_to_control`` to calculate responses relative to
        the appropriate control population.

    batch_key : str, default="batch"
        Name of the column in ``adata.obs`` identifying experimental batches.
        D-SPIN uses batch information when estimating responses relative to
        controls, allowing differences between batches to be accounted for.
        
    num_spin : int or str, default="auto"
        Number of D-SPIN spins to perform during network inference. Each spin
        corresponds to an independent network-inference run used to assess the
        stability of the inferred network.

        If set to an integer, that exact number of spins is performed. Larger
        values generally provide a more stable estimate of the inferred network
        but increase computational cost.

        If set to ``"auto"``, the number of spins is automatically set to the
        number of genes in the input expression matrix, i.e. ``adata.shape[1]``.
        For example, an expression matrix containing 2,000 genes results in
        ``num_spin=2000``.

    method : str, default="pseudo_likelihood"
        Statistical method used for network inference. The
        ``"pseudo_likelihood"`` method estimates the relationships between
        genes using a pseudo-likelihood-based approach.

    directed : bool, default=True
        Whether the inferred network should be directed. If ``True``, an
        edge from gene A to gene B can be distinguished from an edge from
        gene B to gene A. If ``False``, the inferred relationships are
        treated as undirected.

    run_with_matlab : bool, default=False
        Whether to run the network-inference procedure using the MATLAB
        implementation instead of the Python implementation. When ``False``,
        the Python implementation is used.

    params : dict, default={"stepsz": 0.05, "lam_l1_j": 0.01}
        Parameters controlling the network-inference optimization.

        ``stepsz`` : float
            Optimization step size used during parameter estimation.

        ``lam_l1_j`` : float
            L1 regularization strength applied to the inferred network
            parameters. Increasing this value generally promotes a sparser
            network by penalizing small edge weights.

    Returns
    -------
    None
        The function currently does not return the D-SPIN model. Results are
        written to ``save_path`` by D-SPIN.

    Examples
    --------
    >>> build_gene_network(
    ...     file_path="expression_data.h5ad",
    ...     save_path="dspin_results"
    ... )    
    """
    
    if params is None:
        params = {
            "stepsz": 0.05,
            "lam_l1_j": 0.01
        }  
        
    if discretize_params is None:
        discretize_params = {
            "clip_percentile": 100
        }
        
    if num_spin == 'auto':
        num_spin = adata.shape[1]
    
    model = DSPIN(
        adata, 
        save_path, 
        num_spin=num_spin,
        discretize_params=discretize_params
    )  
    
    model.network_inference(
        sample_id_key=sample_id_key,
        method=method,
        directed=directed,
        run_with_matlab=run_with_matlab,
        params=params
    )
    
    model.response_relative_to_control(
        sample_id_key=sample_id_key,
        if_control_key=if_control_key,
        batch_key=batch_key
    )    

    return 


def build_program_network(
    adata,
    save_path,
    num_programs=20,
    num_repeat=10,
    seed=0,
    cluster_key=None,
    sample_id_key="sample_id",
    if_control_key="if_control",
    batch_key="batch",
    # num_spin='auto',
    method="pseudo_likelihood",
    directed=True,
    run_with_matlab=False,
    params_discovery=None,
    params_inference=None
):
    """
    Load a gene-expression dataset, discover gene-expression programs using
    repeated oNMF, perform D-SPIN network inference on the programs, and
    calculate responses relative to control samples.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix containing the single-cell expression data,
        typically loaded from an `.h5ad` file using `scanpy.read_h5ad()`.

    save_path : str
        Directory in which D-SPIN stores intermediate and final results
        generated during gene-program discovery and network inference.

    num_programs : int, default=20
        Number of gene-expression programs to identify during the oNMF
        gene-program discovery step.

    num_repeat : int, default=10
        Number of repeated oNMF runs used to identify stable
        gene-expression programs. Because oNMF is sensitive to
        initialization, repeated runs help identify reproducible programs.

        Increasing this value generally improves the robustness of the
        program discovery step at the cost of additional computation.

    seed : int, default=0
        Random seed used for the stochastic gene-program discovery procedure.
        Setting a fixed value makes the oNMF analysis reproducible, subject
        to the numerical environment and implementation.

    cluster_key : str or None, default=None
        Name of a column in ``adata.obs`` used to define cell clusters for
        balancing during gene-program discovery.

        If ``None``, no cluster-based balancing is requested. If specified,
        D-SPIN can use the corresponding cell-level annotation to avoid
        over-representing highly abundant cell populations during oNMF.

    sample_id_key : str, default="sample_id"
        Name of the column in ``adata.obs`` containing the identifier of
        each sample or perturbation. D-SPIN uses this column to determine
        which cells belong to the same perturbation when constructing the
        perturbation-level representation.

    if_control_key : str, default="if_control"
        Name of the column in ``adata.obs`` indicating whether a cell or
        sample belongs to a control condition. This information is used by
        ``response_relative_to_control`` to calculate responses relative to
        control populations.

    batch_key : str, default="batch"
        Name of the column in ``adata.obs`` identifying experimental batches.
        D-SPIN uses batch information when estimating responses relative to
        controls, allowing differences between batches to be accounted for.

    num_spin : int, default=10
        Number of D-SPIN spins to perform during network inference. Each
        spin corresponds to an independent network-inference run used to
        assess the stability of the inferred network.

    method : str, default="pseudo_likelihood"
        Statistical method used for network inference.

        ``"pseudo_likelihood"`` performs network inference using a
        pseudo-likelihood-based approach and is suitable for relatively
        large networks.

    directed : bool, default=True
        Whether the inferred network should be directed. If ``True``, an
        edge from feature A to feature B can be distinguished from an edge
        from feature B to feature A. If ``False``, the inferred relationships
        are treated as undirected.

    run_with_matlab : bool, default=False
        Whether to run the network-inference procedure using the MATLAB
        implementation instead of the Python implementation. When ``False``,
        the Python implementation is used.

    params : dict or None, default=None
        Parameters controlling the network-inference optimization.

        If ``None``, the following defaults are used:

        ``stepsz`` : float
            Optimization step size. Default is ``0.05``.

        ``lam_l1_j`` : float
            L1 regularization strength applied to the inferred network
            parameters. Increasing this value generally promotes a sparser
            network. Default is ``0.01``.

    Returns
    -------
    None
        The function currently does not return the D-SPIN model. Results are
        written to ``save_path`` by D-SPIN.

    Examples
    --------
    >>> build_program_network(
    ...     file_path="expression_data.h5ad",
    ...     save_path="dspin_results",
    ...     num_programs=20,
    ...     num_repeat=10,
    ...     seed=0
    ... )
    """

    if params_inference is None:
        params_inference = {
            "stepsz": 0.05,
            "lam_l1_j": 0.01
        }
        
    if params_discovery is None:
        params_discovery = {
            "balance_method": None
        }    
        cluster_key = "dspin_dummy_cluster"
        adata.obs["dspin_dummy_cluster"] = 0        

    model = DSPIN(
        adata,
        save_path,
        num_spin=num_programs
    )

    model.gene_program_discovery(
        num_repeat=num_repeat,
        seed=seed,
        cluster_key=cluster_key,
        params=params_discovery
    )

    model.network_inference(
        sample_id_key=sample_id_key,
        method=method,
        directed=directed,
        run_with_matlab=run_with_matlab,
        params=params_inference
    )

    model.response_relative_to_control(
        sample_id_key=sample_id_key,
        if_control_key=if_control_key,
        batch_key=batch_key
    )
    
    print("Program network built successfully.")

    return model