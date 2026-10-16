from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


def plot_molecular_weight_distribution(
    df: pd.DataFrame,
    mw_column: str = "MW",
    bins: int = 50,
    plot_range: tuple = (0, 2000),
    save_path: Optional[str] = None,
):
    """Plot molecular weight distribution and optionally save the figure."""
    if mw_column not in df.columns:
        raise KeyError(f"Column '{mw_column}' not found in dataframe.")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(df[mw_column], bins=bins, range=plot_range)
    ax.set_xlabel("Molecular Weight (Da)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Molecular Weights (<= 2000 Da shown)")

    if save_path:
        output_path = Path(save_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved figure to: {output_path}")

    return fig, ax
