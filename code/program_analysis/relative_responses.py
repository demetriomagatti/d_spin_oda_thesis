import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform


def get_ranked_responses(relative_responses, ascending=True, n=5):
    responses = {}

    for sample in relative_responses.columns:
        ranked = relative_responses[sample].sort_values(ascending=ascending)

        responses[sample] = ranked.head(n)

    for sample, values in responses.items():
        print(f"\n{sample}")
        print(values.to_string())    
    return responses


def check_reciprocal_responses(relative_responses):
    reciprocal_pairs = []

    for sample in relative_responses.columns:
        parts = sample.split("+")
        if len(parts) != 2:
            continue
        gene_1, gene_2 = parts
        if gene_2 == "ctrl":
            reciprocal = f"ctrl+{gene_1}"

            if reciprocal in relative_responses.columns:
                x = relative_responses[sample]
                y = relative_responses[reciprocal]

                correlation = x.corr(y)
                rmse = np.sqrt(np.mean((x - y) ** 2))
                max_abs_difference = np.max(np.abs(x - y))

                reciprocal_pairs.append({
                    "pair_1": sample,
                    "pair_2": reciprocal,
                    "correlation": correlation,
                    "rmse": rmse,
                    "max_abs_difference": max_abs_difference
                })
    reciprocal_reproducibility = pd.DataFrame(reciprocal_pairs)
    print(
        reciprocal_reproducibility
        .sort_values("correlation")
        .to_string(index=False)
    )    
    return reciprocal_reproducibility


def get_program_activation_summary(relative_responses, n=10):
    
    program_summary = []
    for program in relative_responses.index:
        values = relative_responses.loc[program]

        top_samples = values.nlargest(n)

        program_summary.append({
            "program": program,
            "max_response": values.max(),
            "mean_response": values.mean(),
            "sd_response": values.std(),
            "n_positive": (values > 0).sum(),
            "top_1": top_samples.index[0],
            "top_1_score": top_samples.iloc[0],
            "top_2": top_samples.index[1],
            "top_2_score": top_samples.iloc[1],
            "top_3": top_samples.index[2],
            "top_3_score": top_samples.iloc[2]
        })

    program_summary = pd.DataFrame(program_summary)

    print(program_summary.to_string(index=False))
    return program_summary


def chech_program_activation_selectivity(relative_responses):
    program_selectivity = []

    for perturbation in relative_responses.columns:
        values = relative_responses[perturbation].sort_values(ascending=False)

        top_score = values.iloc[0]
        second_score = values.iloc[1]

        program_selectivity.append({
            "perturbation": perturbation,
            "top_program": values.index[0],
            "top_score": top_score,
            "second_program": values.index[1],
            "second_score": second_score,
            "top_minus_second": top_score - second_score
        })

    program_selectivity = pd.DataFrame(program_selectivity)

    program_selectivity = program_selectivity.sort_values(
        "top_minus_second",
        ascending=False
    )

    print(program_selectivity.head(10).to_string(index=False))
    return program_selectivity


def plot_relative_responses(relative_responses, figsize=(16, 10)):
    plt.figure(figsize=figsize)

    plt.imshow(
        relative_responses.values,
        aspect="auto",
        interpolation="none"
    )

    plt.colorbar(label="Relative response")
    plt.xlabel("Perturbation")
    plt.ylabel("Program")
    plt.yticks(
        range(len(relative_responses.index)),
        relative_responses.index
    )

    plt.xticks(
        range(len(relative_responses.columns)),
        relative_responses.columns,
        rotation=90,
        fontsize=7
    )

    plt.tight_layout()
    plt.show()
    return


def plot_program_dendogram(relative_responses):
    program_correlation = relative_responses.T.corr(method="pearson")
    program_distance = 1 - program_correlation

    program_linkage = linkage(
        squareform(program_distance.values, checks=False),
        method="average"
    )

    plt.figure(figsize=(12, 6))

    dendrogram(
        program_linkage,
        labels=program_correlation.index,
        leaf_rotation=90
    )

    plt.ylabel("Distance")
    plt.xlabel("Program")
    plt.tight_layout()
    plt.show()    
    return 