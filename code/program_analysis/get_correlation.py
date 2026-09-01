import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def get_correlation_matrix(model):

    program_activity = pd.DataFrame(
        model.program_representation,
        columns=[f"Program_{i}" for i in range(1, 21)]
    )
    program_correlation = program_activity.corr()
    print(program_correlation.round(2))
    return program_correlation


def plot_correlation_matrix(model):
    program_correlation = get_correlation_matrix(model)
    plt.figure(figsize=(10, 8))
    plt.imshow(program_correlation, vmin=-1, vmax=1, cmap="coolwarm")
    plt.colorbar(label="Pearson correlation")
    plt.xticks(range(20), program_correlation.columns, rotation=90)
    plt.yticks(range(20), program_correlation.index)
    plt.title("Correlation between gene-program activities")
    plt.tight_layout()
    plt.show()
    return 


def get_highest_correlation_duplets(model, direction = "down"):
    program_correlation = get_correlation_matrix(model)
    correlation_pairs = (
        program_correlation
        .where(np.triu(np.ones(program_correlation.shape), k=1).astype(bool))
        .stack()
        .sort_values(ascending=False)
    )
    if direction == "down":
        duplets = correlation_pairs.head(program_correlation.shape[0])
    elif direction == "up":
        duplets = correlation_pairs.sort_values().head(program_correlation.shape[0])
    return duplets