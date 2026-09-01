import pandas as pd 
import gseapy as gp


def ora_enrich_programs(adata, prg, gene_sets="GO_Biological_Process_2023"):
    gene_symbols = adata.var["gene_name"].to_numpy()
    ora = {}
    for key in prg.keys():
        program_genes = prg[key]["gene"].tolist()        
        ora_key = gp.enrich(
            gene_list=program_genes[:100],
            gene_sets=gene_sets,
            background=gene_symbols.tolist(),
            cutoff=0.05,
            no_plot=True,
            verbose=False
        )
        ora[key] = ora_key
    
    return ora


def summarize_ora_results(ora, n_top=5):
    summary_rows = []
    for program_name, result in ora.items():
        results = result.results.copy()

        significant = results[
            results["Adjusted P-value"] < 0.05
        ].sort_values("Adjusted P-value")

        top_terms = significant.head(n_top)
        for _, row in top_terms.iterrows():
            summary_rows.append({
                "program": program_name,
                "term": row["Term"],
                "adjusted_p_value": row["Adjusted P-value"],
                "odds_ratio": row["Odds Ratio"],
                "combined_score": row["Combined Score"]
            })

    ora_summary = pd.DataFrame(summary_rows)    
    return ora_summary