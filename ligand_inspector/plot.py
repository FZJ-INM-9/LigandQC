import pandas as pd
import matplotlib.colors as mcolors


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
        add_axis_indicators: bool = False,
        rotate_xticks: bool = True,
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
            add_axis_indicators (bool, optional): if True, adds arrows to the bar and line plots to indicate to which axis they belong to.

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
            label="Number of " + column,
        )

        # Customize primary axis
        ax.set_xlabel(column, fontsize=12)
        ax.set_ylabel(f"Number of observations", fontsize=12, color="steelblue")
        ax.tick_params(axis="y", labelcolor="steelblue")
        ax.set_xticks(range(len(value_counts)))
        ax.set_xticklabels(
            value_counts.index,
            rotation=25 if rotate_xticks else 0,
            ha="right" if rotate_xticks else "center",
            rotation_mode="anchor" if rotate_xticks else None,
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
        # ax2.text(len(value_counts) * 0.7, 82, "80% Line", color="red", fontsize=10)

        # Add grid for better readability
        ax.grid(True, alpha=0.3)

        if add_axis_indicators:
            # colorful arrows (single arrow per axis)

            # --- parameters ---
            x_span = len(value_counts) - 1
            arrow_len = x_span * 0.15

            # darker colors
            bar_color = mcolors.to_rgb("steelblue")
            bar_color_dark = tuple(max(0, c - 0.15) for c in bar_color)

            line_color = mcolors.to_rgb("red")
            line_color_dark = tuple(max(0, c - 0.15) for c in line_color)

            # --- bar arrow: 1/3 position, pointing right ---
            bar_idx = len(value_counts) // 3
            bar_val = value_counts.values[bar_idx]

            ax.annotate(
                "",
                xytext=(bar_idx, 0.5 * bar_val),
                xy=(bar_idx - arrow_len, 0.5 * bar_val),
                arrowprops=dict(arrowstyle="-|>", color=bar_color_dark, lw=2),
                zorder=10,
            )

            # --- line arrow: smallest index in upper third, pointing right toward y-axis ---
            upper_third_start = len(cumulative_pct) * 2 // 3
            line_idx = upper_third_start
            line_val = cumulative_pct.values[line_idx]

            ax2.annotate(
                "",
                xytext=(line_idx, line_val),
                xy=(line_idx + arrow_len, line_val),
                arrowprops=dict(arrowstyle="-|>", color=line_color_dark, lw=2),
                zorder=10,
            )
