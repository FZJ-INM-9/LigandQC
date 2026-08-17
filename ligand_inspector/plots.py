import pandas as pd


class ParetoPlot:
    def __init__(self, data: pd.DataFrame):
        """
        Initializes the ParetoPlot class with data.

        Args:
            data (pd.DataFrame): DataFrame containing the data to plot. It is assumed that the DataFrame represents an instance table with a single column of interest, where each value will be treated as categorical.
        """
        self.data = data

    def plot(
        self,
        ax,
        column: str,
        truncate: int = None,
        sort_by_category: bool = False,
        ascending: bool = True,
    ):
        """
        Generates and saves a Pareto chart for the specified column.

        Args:
            ax (matplotlib.pyplot.axis): Axis of plot.
            column (str): Column name to create the Pareto chart for.
            truncate (int, optional): Number of top entries to show. Remaining entries
                                    will be summarized as "Others". If None, all entries are shown.
            sort_by_category (bool, optional): sorts the values of the x-axis by category name/value.
            ascending (bool, optional): if sort_by_category is True, this boolean is applied to the sorting.

        Returns:
            None
        """
        if column not in self.data.columns:
            raise ValueError(f"Column '{column}' not found in data")

        # Remove NaN values and calculate counts
        series_clean = self.data[column].dropna()
        value_counts = series_clean.value_counts()

        # Handle truncation
        if (truncate is not None) and (len(value_counts) > truncate):
            # Keep top (truncate-1) entries and sum the rest as "Others"
            top_counts = value_counts.iloc[: truncate - 1]
            others_count = value_counts.iloc[truncate - 1 :].sum()
            n_others = len(value_counts) - (truncate - 1)
            others_label = f"Others ({n_others})"

            # Create new series with "Others" category
            truncated_counts = pd.concat(
                [top_counts, pd.Series([others_count], index=[others_label])]
            )
            value_counts = truncated_counts

        # Sort the counts if required
        if sort_by_category:
            if "others_label" in locals():
                others = value_counts.pop(others_label)
            else:
                others = None
            value_counts = value_counts.sort_index(ascending=ascending)
            if others is not None:
                value_counts[others_label] = others

        # Calculate cumulative percentages
        cumulative_pct = (value_counts.cumsum() / value_counts.sum()) * 100

        # Plot the bars (counts) on primary axis
        bars = ax.bar(
            range(len(value_counts)),
            value_counts.values,
            color="steelblue",
            alpha=0.7,
            label="Count",
        )

        # Customize primary axis
        ax.set_xlabel(column, fontsize=12)
        ax.set_ylabel("Count", fontsize=12, color="steelblue")
        ax.tick_params(axis="y", labelcolor="steelblue")
        ax.set_xticks(range(len(value_counts)))
        ax.set_xticklabels(
            value_counts.index,
            rotation=25,
            ha="right",
            rotation_mode="anchor",
            va="top",
        )
        ax.set_xlim(-0.5, len(value_counts) - 0.5)

        # Create secondary axis for cumulative percentage
        ax2 = ax.twinx()

        # Plot the cumulative percentage line
        line = ax2.plot(
            range(len(cumulative_pct)),
            cumulative_pct.values,
            color="red",
            marker="o",
            linewidth=2,
            markersize=4,
            label="Cumulative %",
        )

        # Customize secondary axis
        ax2.set_ylabel("Cumulative Percentage (%)", fontsize=12, color="red")
        ax2.tick_params(axis="y", labelcolor="red")
        ax2.set_ylim(0, 100)

        # Add horizontal line at 80% for Pareto principle
        ax2.axhline(y=80, color="red", linestyle="--", alpha=0.7, linewidth=1)
        ax2.text(len(value_counts) * 0.7, 82, "80% Line", color="red", fontsize=10)

        # Add grid for better readability
        ax.grid(True, alpha=0.3)

        # Create combined legend
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(
            lines1 + lines2,
            labels1 + labels2,
            ncols=2,
            bbox_to_anchor=(0.5, 1.0),
            loc="lower center",
        )
