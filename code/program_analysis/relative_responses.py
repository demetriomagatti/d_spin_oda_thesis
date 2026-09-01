import pandas as pd
import numpy as np


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