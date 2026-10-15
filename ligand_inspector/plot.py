import pandas as pd
import matplotlib.colors as mcolors
from matplotlib_venn import venn2, venn3
import seaborn as sns
from ligand_inspector.const import NUM_METRICS

class VennPlot:
    def __init__(self, datasets, labels):
        """
        Initializes the VennPlot class with datasets and labels.

        Args:
            datasets (list): List of datasets to compare.
            labels (list): List of labels corresponding to the datasets.
        """
        self.datasets = datasets
        self.labels = labels

    def plot(self, ax=None):
        """
        Plots a Venn diagram for the provided datasets.

        Args:
            ax (matplotlib.axes.Axes, optional): Axis to plot on. If None, a new figure and axis are created.

        Returns:
            None
        """
        if len(self.datasets) == 2:
            venn2([set(dataset) for dataset in self.datasets], set_labels=self.labels, ax=ax)
        elif len(self.datasets) == 3:
            venn3([set(dataset) for dataset in self.datasets], set_labels=self.labels, ax=ax)

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


# class ComparativeParetoPlot:
#     def __init__(self, df: pd.DataFrame):
#         """
#         Initializes the ComparativeParetoPlot class with data.

#         Args:
#             df (pd.DataFrame): DataFrame containing the data to plot.
#         """
#         self.df = df

#     def plot(
#         self,
#         ax,
#         column: str,
#         truncate: int = None,
#         sort_by_category: bool = False,
#         ascending: bool = True,
#         rotate_xticks: bool = True,
#     ):
#         """
#         plots a comparative Pareto chart for the specified column across multiple datasets. It shows the distribution of values in the specified column for each dataset in a grouped bar format, along with cumulative percentage lines for each dataset.
#         """

#         vc = []
#         cpct = []

#         dataset_labels = sorted(list(self.df['Dataset'].unique()))
#         for i, label in enumerate(dataset_labels):
#             dataset = self.df[self.df['Dataset'] == label]
#             if column not in dataset.columns:
#                 raise ValueError(f"Column '{column}' not found in dataset {i+1}")

#             # Remove NaN values and calculate counts
#             series_clean = dataset[column].dropna()
#             value_counts = series_clean.value_counts()

#             # Handle truncation
#             if (truncate is not None) and (len(value_counts) > truncate):
#                 # Keep top (truncate-1) entries and sum the rest as "Others"
#                 top_counts = value_counts.iloc[: truncate - 1]
#                 others_count = value_counts.iloc[truncate - 1 :].sum()
#                 n_others = len(value_counts) - (truncate - 1)
#                 others_label = f"Others ({n_others})"

#                 # Create new series with "Others" category
#                 truncated_counts = pd.concat(
#                     [top_counts, pd.Series([others_count], index=[others_label])]
#                 )
#                 value_counts = truncated_counts

#             # Sort the counts if required
#             if sort_by_category:
#                 if "others_label" in locals():
#                     others = value_counts.pop(others_label)
#                 else:
#                     others = None
#                 value_counts = value_counts.sort_index(ascending=ascending)
#                 if others is not None:
#                     value_counts[others_label] = others

#             vc.append(value_counts)

#             # Calculate cumulative percentages
#             cumulative_pct = (value_counts.sort_index(ascending=False).cumsum() / value_counts.sum()) * 100
#             cpct.append(cumulative_pct)

#         dfs = [vc[i].to_frame(name=dataset_labels[i]) for i in range(len(vc))]
#         df = pd.concat(dfs, axis=1).reset_index()
#         # df = df.melt(id_vars=dataset_labels, var_name="dataset label")
#         df = pd.melt(
#             df,
#             id_vars=[df.columns[0]],  # Keep this column as is
#             value_vars=dataset_labels,  # Columns to melt into rows
#             var_name="dataset label",  # Name for the new column holding column names
#             value_name="instance count",  # Name for the new column holding values
#         )
#         df.fillna(0, inplace=True)  # Fill NaN values with 0 for plotting
#         df.sort_values(by=[df.columns[0], "dataset label"], ascending=[False, True], inplace=True)  # Sort by the first column and then by dataset label

#         # Create grouped bar plot (counts) on primary axis
#         sns.barplot(df,
#             color="steelblue",
#             alpha=0.7,
#             hue="dataset label",
#             x=df.columns[0],
#             y="instance count",
#             ax=ax,
#             order=sorted(df[df.columns[0]].unique(), reverse=True)
#         )


#         # Customize primary axis
#         ax.set_xlabel(column, fontsize=12)
#         ax.set_ylabel(f"Number of observations", fontsize=12, color="steelblue")
#         ax.tick_params(axis="y", labelcolor="steelblue")
#         ax.set_xticks(range(len(value_counts)))
#         ax.set_xticklabels(
#             ax.get_xticklabels(),
#             rotation=25 if rotate_xticks else 0,
#             ha="right" if rotate_xticks else "center",
#             rotation_mode="anchor" if rotate_xticks else None,
#             va="top",
#         )
#         ax.set_xlim(-0.5, len(value_counts) - 0.5)

#         # Create secondary axis for cumulative percentage
#         ax2 = ax.twinx()
        
#         # create cumulative percentage DataFrame for seaborn lineplot
#         filled_dfs = []
#         for i, label in enumerate(dataset_labels):
#             g = df.groupby("dataset label").get_group(label)
#             g = g.copy()
#             # fill missing values for 'instance count' with 0
#             g = g.set_index(df.columns[0]).reindex(range(0, NUM_METRICS + 1), fill_value=0).reset_index()
#             g["dataset label"] = label
#             g.sort_values(by=df.columns[0], ascending=False, inplace=True)
#             g["cumulative percentage"] = (g["instance count"].cumsum() / g["instance count"].sum()) * 100
#             filled_dfs.append(g)
#         df = pd.concat(filled_dfs, ignore_index=True)

#         sns.lineplot(
#             data=df,
#             x=df.columns[0],
#             y="cumulative percentage",
#             hue="dataset label",
#             ax=ax2,
#             marker="o",
#         )

#         # for i, label in enumerate(dataset_labels):
#         #     print(label)
#         #     print(cpct[i])
#         #     # Plot the cumulative percentage line for each dataset
#         #     ax2.plot(
#         #         range(len(cpct[i])),
#         #         cpct[i].values,
#         #         marker="o",
#         #         linewidth=2,
#         #         markersize=4,
#         #         label=f"Cumulative % ({label})",
#         #     )

#         # # Plot the cumulative percentage line
#         # line = ax2.plot(
#         #     range(len(cumulative_pct)),
#         #     cumulative_pct.values,
#         #     color="red",
#         #     marker="o",
#         #     linewidth=2,
#         #     markersize=4,
#         #     label="Cumulative %",
#         # )

#         # # Customize secondary axis
#         # ax2.set_ylabel("Cumulative Percentage (%)", fontsize=12, color="red")
#         # ax2.tick_params(axis="y", labelcolor="red")
#         # ax2.set_ylim(0, 100)

#         # # Add horizontal line at 80% for Pareto principle
#         # ax2.axhline(y=80, color="red", linestyle="--", alpha=0.7, linewidth=1)
#         # # ax2.text(len(value_counts) * 0.7, 82, "80% Line", color="red", fontsize=10)

#         # # Add grid for better readability
#         # ax.grid(True, alpha=0.3)


from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


class ComparativeParetoPlot:
    def __init__(self, data: pd.DataFrame, dataset_column: str = "dataset"):
        """
        Initializes the GroupedParetoPlot class with data.

        Args:
            data (pd.DataFrame): Instance table containing the column of interest
                and a column identifying which dataset each row belongs to.
            dataset_column (str): Name of the column that distinguishes datasets.
        """
        if dataset_column not in data.columns:
            raise ValueError(f"Column '{dataset_column}' not found in data")
        self.data = data
        self.dataset_column = dataset_column

    def plot(
        self,
        ax,
        column: str,
        truncate: Optional[int] = None,
        sort_by_category: bool = False,
        ascending: bool = True,
        order_by: Optional[str] = None,
        normalize: bool = False,
        dataset_order: Optional[Sequence[str]] = None,
        colors: Optional[Sequence] = None,
        rotate_xticks: bool = True,
    ):
        """
        Draws a grouped Pareto chart: one group of bars per category (one bar per
        dataset) plus one cumulative-percentage line per dataset, left to right.

        Args:
            ax (matplotlib.axes.Axes): Axis to draw on.
            column (str): Column containing the categories.
            truncate (int, optional): Number of entries to show. The remaining
                categories are summed into "Others (n)". If None, all are shown.
            sort_by_category (bool): Sort x-axis by category name instead of count.
            ascending (bool): Sort direction if sort_by_category is True.
            order_by (str, optional): Dataset name whose counts define the
                category order. If None, the total across all datasets is used.
            normalize (bool): If True, bars show the percentage of each dataset
                (useful when datasets differ in size) instead of raw counts.
            dataset_order (sequence, optional): Order of the datasets within each group.
            colors (sequence, optional): One color per dataset.
            rotate_xticks (bool): Rotate x tick labels.

        Returns:
            None
        """
        if column not in self.data.columns:
            raise ValueError(f"Column '{column}' not found in data")

        ds_col = self.dataset_column
        df = self.data[[ds_col, column]].dropna()

        # Count table: rows = categories, columns = datasets
        counts = df.groupby([column, ds_col]).size().unstack(ds_col, fill_value=0)
        if dataset_order is not None:
            counts = counts.reindex(columns=list(dataset_order), fill_value=0)
        if order_by is not None and order_by not in counts.columns:
            raise ValueError(f"order_by '{order_by}' is not a dataset")

        # Order categories by size (descending), shared across all datasets
        sort_key = counts[order_by] if order_by is not None else counts.sum(axis=1)
        counts = counts.loc[sort_key.sort_values(ascending=False, kind="stable").index]

        # Truncate: keep top (truncate-1) and merge the rest into "Others"
        others_label = None
        if (truncate is not None) and (len(counts) > truncate):
            n_others = len(counts) - (truncate - 1)
            others_label = f"Others ({n_others})"
            others_row = counts.iloc[truncate - 1 :].sum().to_frame(name=others_label).T
            counts = pd.concat([counts.iloc[: truncate - 1], others_row])

        # Optionally sort by category name, keeping "Others" last
        if sort_by_category:
            if others_label is not None:
                others_row = counts.loc[[others_label]]
                counts = counts.drop(index=others_label)
            counts = counts.sort_index(ascending=ascending)
            if others_label is not None:
                counts = pd.concat([counts, others_row])

        # Per-dataset cumulative percentage along the shared category order
        totals = counts.sum()
        cumulative_pct = counts.cumsum() / totals * 100
        bar_values = counts / totals * 100 if normalize else counts

        n_cat, n_ds = counts.shape
        x = np.arange(n_cat)
        width = 0.8 / n_ds
        offsets = (np.arange(n_ds) - (n_ds - 1) / 2) * width
        if colors is None:
            cmap = plt.get_cmap("tab10")
            colors = [cmap(i % 10) for i in range(n_ds)]

        # Grouped bars on primary axis
        for i, ds in enumerate(counts.columns):
            ax.bar(
                x + offsets[i],
                bar_values[ds].values,
                width=width,
                color=colors[i],
                alpha=0.7,
                label=str(ds),
            )

        ax.set_xlabel(column, fontsize=12)
        ax.set_ylabel(
            "Percentage of observations (%)" if normalize else "Number of observations",
            fontsize=12,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(
            counts.index,
            rotation=25 if rotate_xticks else 0,
            ha="right" if rotate_xticks else "center",
            rotation_mode="anchor" if rotate_xticks else None,
            va="top",
        )
        ax.set_xlim(-0.5, n_cat - 0.5)

        # Cumulative lines on secondary axis, one per dataset, at group centers
        ax2 = ax.twinx()
        for i, ds in enumerate(counts.columns):
            ax2.plot(
                x,
                cumulative_pct[ds].values,
                color=colors[i],
                marker="o",
                linewidth=2,
                markersize=4,
            )
        ax2.set_ylabel("Cumulative Percentage (%)", fontsize=12)
        ax2.set_ylim(0, 100)
        ax2.axhline(y=80, color="gray", linestyle="--", alpha=0.7, linewidth=1)

        # Legend: one entry per dataset + one generic entry for the lines
        handles, labels = ax.get_legend_handles_labels()
        handles.append(Line2D([0], [0], color="black", marker="o", markersize=4, linewidth=2))
        labels.append("Cumulative %")
        ax.legend(handles, labels)

        ax.grid(True, alpha=0.3)