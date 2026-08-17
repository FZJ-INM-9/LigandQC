"""
RMSZ
====

RMSZ is a Python package to explore a rationale for categorizing the RMSZ value of bond lengths and angles in a ligand into the categories "good", "ok" and "bad".
In this script, we will explore the RMSZ of a standard deviation by
1) generating a dataset that adheres to a normal distribution with a given mean and standard deviation called "reference"
2) calculating the mean and standard deviation of the dataset "reference"
3) calculating the Z-score of each value in the dataset "reference"
4) calculating the RMSZ of the dataset "reference"
5) repeat steps 1-4 for 3 different datasets with the same mean but standard deviations of 0.5, 2.0 and 3.0, calculate the RMSZ with respect to the reference dataset
6) plot the results

"""

from scipy.stats import f
import numpy as np
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import seaborn as sns
from .LigandDownloader import LigandDownloader
from scipy.stats import gaussian_kde


def gaussian(x, mean, std_dev):
    """Calculate the Gaussian function."""
    return (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(
        -0.5 * ((x - mean) / std_dev) ** 2
    )


def generate_dataset(mean, std_dev, size=1000):
    """Generate a dataset with a normal distribution."""
    return np.random.normal(mean, std_dev, size)


def calculate_zscore(dataset, reference_mean=0, reference_std_dev=1):
    """Calculate the Z-score of a dataset."""
    z_scores = (dataset - reference_mean) / reference_std_dev
    return z_scores


def calculate_rmsz(dataset, reference_std_dev=1.0, reference_mean=0):
    """Calculate the RMSZ of a dataset."""
    z_scores = calculate_zscore(dataset, reference_mean, reference_std_dev)
    rmsz = np.sqrt(np.mean(z_scores**2))
    return rmsz


def calculate_percentage(dataset, reference_std_dev=1.0, reference_mean=0):
    """Calculate the percentage of values within 2 reference standard deviations of the reference mean."""
    lower_bound = reference_mean - 2 * reference_std_dev
    upper_bound = reference_mean + 2 * reference_std_dev
    percentage_within_bounds = (
        np.sum((dataset >= lower_bound) & (dataset <= upper_bound)) / len(dataset) * 100
    )
    return percentage_within_bounds


class RMSZInterpreter:
    def __init__(self):
        self.reference_mean = 0
        self.reference_std_dev = 1.0
        self.reference_dataset_size = 1e7
        self.reference_dataset = np.random.normal(
            self.reference_mean,
            self.reference_std_dev,
            int(self.reference_dataset_size),
        )
        self.ok_bad_threshold = 2.0
        self.precomputed_rmsz_thresholds = {None: 1.0}

    def interpret_rmsz(self, n, rmsz):
        """
        Interpret the RMSZ value of a distribution.

        Parameters:
            n (int): sample size
            rmsz (float): The RMSZ value of the distribution.
            mode (str): The mode of interpretation ("bond" or "angle").

        Returns:
            str: Interpretation of the RMSZ value.
        """
        rmsz_threshold = self.find_rmsz_threshold(n)

        if rmsz <= rmsz_threshold:
            return "good"
        elif rmsz <= self.ok_bad_threshold:
            return "ok"
        else:
            return "bad"

    def perform_f_test(self, rmsz, n, alpha=0.05):
        """
        Perform an F-test to compare variances of two samples.

        Parameters:
            rmsz (float): RMSZ value of the dataset.
            n (int): Sample size of the dataset.
            alpha (float): Significance level (default is 0.05).

        Returns:
            dict: F-statistic, p-value, and test result.
        """
        # Calculate the F-statistic
        var1 = rmsz**2  # Variance of the first sample (RMSZ)
        var2 = (
            self.reference_std_dev**2
        )  # Variance of the second sample (reference std dev)
        F_statistic = var1 / var2 if var1 > var2 else var2 / var1

        # Degrees of freedom for each sample
        df1 = n - 1
        df2 = self.reference_dataset_size - 1

        # Calculate the p-value
        if var1 > var2:
            p_value = f.sf(F_statistic, df1, df2)  # Survival function (1 - CDF)
        else:
            p_value = f.sf(F_statistic, df2, df1)

        # Determine if the result is statistically significant
        significant = p_value < alpha

        return {
            "F-statistic": F_statistic,
            "p-value": p_value,
            "significant": significant,
        }

    def find_rmsz_threshold(self, n, alpha=0.05):
        """
        Find the rmsz threshold for a given sample size n for the F-test to be significant using the bisection method.

        Args:
            n (int): Sample size.
            alpha (float): Significance level (default is 0.05).

        Returns:
            float: The rmsz threshold.
        """
        rmsz_threshold = self.precomputed_rmsz_thresholds.get(n, None)
        if rmsz_threshold is not None:
            return rmsz_threshold

        # Initialize bisection method
        rmsz_start = 3.0
        rmsz_end = 1.01
        if n >= 13526:
            raise ValueError(
                "Sample size n is too large. Please use a smaller sample size."
            )

        # Perform bisection method to find the rmsz threshold
        rmsz_values = []
        while rmsz_start - rmsz_end > 1e-6:
            mid_rmsz = (rmsz_start + rmsz_end) / 2
            result = self.perform_f_test(mid_rmsz, n, alpha)

            if not result["significant"]:
                rmsz_end = mid_rmsz
            else:
                rmsz_start = mid_rmsz
            rmsz_values.append(mid_rmsz)

        # Save the result and return
        rounded_value = round(mid_rmsz, 5)
        self.precomputed_rmsz_thresholds[n] = rounded_value
        return rounded_value

    def plot_rationale_categorization(self, output_path=None):
        """
        Plot the rationale for flagging ligands as "good" based on the topology and RMSZ value.

        Args:
            output_path (str): Path to save the plot. If None, the plot will be shown.

        Returns:
            None
        """
        n_values = range(2, 201)

        rmsz_values = [self.find_rmsz_threshold(n) for n in n_values]

        plt.figure(figsize=(7, 6))
        plt.plot(n_values, rmsz_values, label="RMSZ Threshold", color="blue")

        opacity = 0.3
        plt.axhline(
            y=self.ok_bad_threshold,
            color="red",
            linestyle="--",
            label=f"RMSZ = {self.ok_bad_threshold}",
        )
        # color between the rmsz values and 3.0
        plt.fill_between(
            n_values,
            y1=4.0,
            y2=self.ok_bad_threshold,
            color="red",
            alpha=opacity,
            label="bad",
        )
        plt.fill_between(
            n_values,
            rmsz_values,
            self.ok_bad_threshold,
            where=(np.array(rmsz_values) < self.ok_bad_threshold),
            color="yellow",
            alpha=opacity,
            label="ok",
        )
        # color between the rmsz values and 0
        plt.fill_between(
            n_values,
            rmsz_values,
            0,
            where=(np.array(rmsz_values) > 0),
            color="green",
            alpha=opacity,
            label="good",
        )

        # fill above

        plt.ylim(0, self.ok_bad_threshold + 0.99)
        plt.xlim(0, 200)
        plt.title("Rationale for flagging Ligands as 'Good', 'Ok', or 'Bad'")
        plt.xlabel("Sample Size (n)")
        plt.ylabel("RMSZ Threshold")
        plt.legend()
        plt.grid()
        plt.tight_layout()
        if output_path:
            plt.savefig(output_path)
        else:
            plt.show()

    def plot_rationale_ok_threshold(self, output_path=None):
        """
        Plot the rationale for flagging ligands as "ok" based on the topology and RMSZ value.

        Args:
            output_path (str): Path to save the plot. If None, the plot will be shown.

        Returns:
            None
        """

        fig, axs = plt.subplots(2, 1, figsize=(7, 6))
        fig.suptitle(
            f"Comparison of two Normal Distributions with RMSZ 1.0 and {self.ok_bad_threshold}"
        )
        mean = 0
        reference_std_dev = 1.0
        dataset_std_dev = self.ok_bad_threshold
        x = np.linspace(-7, 7, 1_000_000)
        gaussian1 = gaussian(x, mean, reference_std_dev)
        gaussian2 = gaussian(x, mean, dataset_std_dev)

        inters_x = np.minimum(gaussian1, gaussian2)
        opacity = 0.3

        axs[0].plot(x, gaussian1, color="b", label="RMSZ 1.0")
        axs[0].fill_between(x, gaussian1, 0, color="b", alpha=opacity)
        axs[0].plot(x, gaussian2, color="orange", label=f"RMSZ {self.ok_bad_threshold}")
        axs[0].fill_between(x, gaussian2, 0, color="orange", alpha=opacity)
        axs[0].plot(x, inters_x, color="r")
        axs[0].fill_between(
            x,
            inters_x,
            0,
            facecolor="none",
            edgecolor="r",
            hatch="xx",
            label="intersection",
        )
        axs[0].set_xlim(-7, 7)
        axs[0].set_ylim([0, 0.5])

        area_inters_x = np.trapezoid(inters_x, x)

        handles, labels = axs[0].get_legend_handles_labels()
        labels[2] += f": {area_inters_x * 100:.1f} %"
        axs[0].legend(handles, labels)

        x2 = np.array([xi for xi in x if abs(xi) <= self.ok_bad_threshold])
        gaussian1_x2 = gaussian(x2, mean, reference_std_dev)
        gaussian2_x2 = gaussian(x2, mean, dataset_std_dev)

        axs[1].plot(x, gaussian1, color="b", label="RMSZ 1.0")
        axs[1].fill_between(
            x2, gaussian1_x2, 0, color="b", alpha=opacity, label="RMSZ 1.0"
        )
        axs[1].plot(x, gaussian2, color="orange", label=f"RMSZ {self.ok_bad_threshold}")
        axs[1].fill_between(
            x2,
            gaussian2_x2,
            0,
            color="orange",
            alpha=opacity,
            label=f"RMSZ {self.ok_bad_threshold}",
        )
        axs[1].vlines(
            [self.ok_bad_threshold, -self.ok_bad_threshold],
            0,
            0.5,
            linestyles="dashed",
            colors="black",
        )
        axs[1].set_xlim(-7, 7)
        axs[1].set_ylim([0, 0.5])

        area_inters_1 = np.trapezoid(gaussian1_x2, x2)
        area_inters_2 = np.trapezoid(gaussian2_x2, x2)

        handles, labels = axs[1].get_legend_handles_labels()
        labels[1] += f": {area_inters_1 * 100:.1f} %"
        labels[3] += f": {area_inters_2 * 100:.1f} %"
        axs[1].legend(handles, labels)
        # axs[1].title('RMSZ Distribution Comparison')
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path)
        else:
            plt.show()
