"""
interpret
=========

The `interpret` module is a Python module to explore a rationale for categorizing the RMSZ value of bond lengths and angles in a ligand into the categories 'good', 'ok' and 'bad'.
The rationale is based on the F-test, which compares the variance of the ligand's bond lengths or angles to a reference distribution with a known variance. The exact rationale is described in the documentation.

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import f

from ligandqc.const import *
from ligandqc.log import logger


def gaussian(x, mean, std_dev):
    """Calculate the Gaussian function."""
    return (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std_dev) ** 2)


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
    percentage_within_bounds = np.sum((dataset >= lower_bound) & (dataset <= upper_bound)) / len(dataset) * 100
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
        self.ok_bad_threshold = THRESH_RMSZ_BAD_UPPER
        self.precomputed_rmsz_thresholds = {None: 1.0}
        self.mode = "single-sided"

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
        if rmsz_threshold is None:
            return CATEGORY_UNKNOWN

        if rmsz <= rmsz_threshold:
            return CATEGORY_GOOD
        elif rmsz <= self.ok_bad_threshold:
            return CATEGORY_OK
        else:
            return CATEGORY_BAD

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
        factor = (n / (n - 1)) if n > 1 else 2
        var1 = factor * rmsz**2  # Variance of the first sample (RMSZ)
        var2 = self.reference_std_dev**2  # Variance of the second sample (reference std dev)
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

    def find_rmsz_threshold(self, n, alpha=0.05, tol=1e-4):

        rmsz_threshold = self.precomputed_rmsz_thresholds.get(n, None)
        if rmsz_threshold is not None:
            return rmsz_threshold

        try:
            calculated_threshold = self._calculate_rmsz_threshold(n, alpha, tol, lb_start=1, ub_start=20)
        except Exception as e:
            logger.warning(f"Could not calculate RMSZ threshold for sample size {n}: {e}")
            # log the exception text
            logger.info(f"Exception text: {e}")
            calculated_threshold = None
        self.precomputed_rmsz_thresholds[n] = calculated_threshold
        return calculated_threshold

    def _calculate_rmsz_threshold(self, n, alpha=0.05, tol=1e-4, lb_start=1, ub_start=20):
        """
        Find the rmsz threshold for a given sample size n for the F-test to be significant using the bisection method.

        Args:
            n (int): Sample size.
            alpha (float): Significance level (default is 0.05).
            tol (float): Tolerance for the bisection method (default is 1e-4).

        Returns:
            float: The rmsz threshold.
        """
        assert n, f"Sample size n must be provided. Got n={n}"
        assert n > 1, f"Sample size n must be greater than 1. Got n={n}"

        # Initialize bisection method
        if n >= 13526:
            raise ValueError("Sample size n is too large. Please use a smaller sample size.")

        # Perform bisection method to find the rmsz threshold
        lb = lb_start
        ub = ub_start
        lb_significant = self.perform_f_test(lb, n, alpha)["significant"]
        ub_significant = self.perform_f_test(ub, n, alpha)["significant"]

        assert (
            lb_significant != ub_significant
        ), f"Bisection method requires one bound to be significant and the other not. Please check the input parameters. lb: {lb}, ub: {ub}, lb_significant: {lb_significant}, ub_significant: {ub_significant}, n: {n}, alpha: {alpha}"

        # accuracy determined by tol
        while abs(ub - lb) > tol:
            mid = (lb + ub) / 2
            result = self.perform_f_test(mid, n, alpha)

            # assign the new bound such that always lb_significant != ub_significant.
            if (result["significant"] and lb_significant and not ub_significant) or (
                not result["significant"] and not lb_significant and ub_significant
            ):
                lb = mid
                lb_significant = result["significant"]
                continue
            elif (result["significant"] and not lb_significant and ub_significant) or (
                not result["significant"] and lb_significant and not ub_significant
            ):
                ub = mid
                ub_significant = result["significant"]
                continue

            assert False, "Bisection method failed to converge. Please check the input parameters."

        # Save the result and return
        rounded_value = round(mid, int(np.log10(1 / tol)) - 1)
        return rounded_value

    def plot_rationale_categorization(self, output_path=None, alpha=0.05):
        """
        Plot the rationale for flagging ligands as 'good', 'ok', or 'bad' based on the topology and RMSZ value.

        Args:
            output_path (str): Path to save the plot. If None, the plot will be shown.
            alpha (float): Significance level for the F-test (default is 0.05).

        Returns:
            None
        """
        n_values = range(2, 201)

        rmsz_values = [self.find_rmsz_threshold(n, alpha=alpha) for n in n_values]

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
            color=COLORS_MATPLOTLIB[CATEGORY_BAD],
            alpha=opacity,
            label=CATEGORY_BAD,
        )
        plt.fill_between(
            n_values,
            rmsz_values,
            self.ok_bad_threshold,
            where=(np.array(rmsz_values) < self.ok_bad_threshold),
            color=COLORS_MATPLOTLIB[CATEGORY_OK],
            alpha=opacity,
            label=CATEGORY_OK,
        )
        # color between the rmsz values and 0
        plt.fill_between(
            n_values,
            rmsz_values,
            0,
            where=(np.array(rmsz_values) > 0),
            color=COLORS_MATPLOTLIB[CATEGORY_GOOD],
            alpha=opacity,
            label=CATEGORY_GOOD,
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
        fig.suptitle(f"Comparison of two Normal Distributions with RMSZ 1.0 and {self.ok_bad_threshold}")
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
        axs[1].fill_between(x2, gaussian1_x2, 0, color="b", alpha=opacity, label="RMSZ 1.0")
        axs[1].plot(x, gaussian2, color="orange", label=f"RMSZ {self.ok_bad_threshold}")
        axs[1].fill_between(
            x2,
            gaussian2_x2,
            0,
            color="orange",
            alpha=opacity,
            label=f"RMSZ {self.ok_bad_threshold}",
        )
        for i in [-self.ok_bad_threshold, self.ok_bad_threshold]:
            axs[1].axvline(
                x=i,
                linestyle="dashed",
                color="black",
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


class RMSZInterpreterDoubleSided(RMSZInterpreter):
    def __init__(self):
        self.reference_mean = 0
        self.reference_std_dev = 1.0
        self.reference_dataset_size = 1e7
        self.reference_dataset = np.random.normal(
            self.reference_mean,
            self.reference_std_dev,
            int(self.reference_dataset_size),
        )
        self.ok_bad_upper_threshold = THRESH_RMSZ_BAD_UPPER
        self.ok_bad_lower_threshold = THRESH_RMSZ_BAD_LOWER
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
        rmsz_lower_threshold, rmsz_upper_threshold = self.find_rmsz_threshold(n)

        if rmsz <= rmsz_upper_threshold and rmsz >= rmsz_lower_threshold:
            return CATEGORY_GOOD
        elif (rmsz < self.ok_bad_upper_threshold and rmsz > rmsz_upper_threshold) or (
            rmsz > self.ok_bad_lower_threshold and rmsz < rmsz_lower_threshold
        ):
            return CATEGORY_OK
        else:
            return CATEGORY_BAD

    def find_rmsz_threshold(self, n, alpha=0.05, tol=1e-4):
        rmsz_threshold = self.precomputed_rmsz_thresholds.get(n, None)
        if rmsz_threshold is not None:
            return rmsz_threshold

        # If not precomputed, calculate thresholds
        lower_threshold = self._calculate_rmsz_threshold(n, alpha, tol, lb_start=0.01, ub_start=1)
        upper_threshold = self._calculate_rmsz_threshold(n, alpha, tol, lb_start=1, ub_start=20)

        self.precomputed_rmsz_thresholds[n] = (lower_threshold, upper_threshold)
        return (lower_threshold, upper_threshold)

    def plot_rationale_categorization(self, output_path=None, alpha=0.05):
        """
        Plot the rationale for flagging ligands as 'good', 'ok', or 'bad' based on the topology and RMSZ value.

        Args:
            output_path (str): Path to save the plot. If None, the plot will be shown.
            alpha (float): Significance level for the F-test (default is 0.05).

        Returns:
            None
        """
        n_values = range(2, 201)

        rmsz_values = [self.find_rmsz_threshold(n, alpha=alpha) for n in n_values]
        rmsz_values_upper = [rmsz[1] for rmsz in rmsz_values]
        rmsz_values_lower = [rmsz[0] for rmsz in rmsz_values]

        plt.figure(figsize=(7, 6))
        plt.plot(n_values, rmsz_values_upper, label="RMSZ Upper Threshold", color="blue")
        plt.plot(n_values, rmsz_values_lower, label="RMSZ Lower Threshold", color="green")

        opacity = 0.3
        plt.axhline(
            y=self.ok_bad_upper_threshold,
            color="red",
            linestyle="--",
            label=f"RMSZ = {self.ok_bad_upper_threshold}",
        )
        # color between the rmsz values and 4.0
        plt.fill_between(
            n_values,
            y1=4.0,
            y2=self.ok_bad_upper_threshold,
            color=COLORS_MATPLOTLIB[CATEGORY_BAD],
            alpha=opacity,
            label=CATEGORY_BAD,
        )
        plt.fill_between(
            n_values,
            rmsz_values_upper,
            self.ok_bad_upper_threshold,
            where=(np.array(rmsz_values_upper) < self.ok_bad_upper_threshold),
            color=COLORS_MATPLOTLIB[CATEGORY_OK],
            alpha=opacity,
            label=CATEGORY_OK,
        )
        # color between upper and lower thresholds
        plt.fill_between(
            n_values,
            rmsz_values_upper,
            rmsz_values_lower,
            where=(np.array(rmsz_values_upper) > 0),
            color=COLORS_MATPLOTLIB[CATEGORY_GOOD],
            alpha=opacity,
            label=CATEGORY_GOOD,
        )

        # color between the lower rmsz threshold and self.ok_bad_lower_threshold
        plt.fill_between(
            n_values,
            rmsz_values_lower,
            self.ok_bad_lower_threshold,
            where=(np.array(rmsz_values_lower) > self.ok_bad_lower_threshold),
            color=COLORS_MATPLOTLIB[CATEGORY_OK],
            alpha=opacity,
        )

        # color below the lower rmsz threshold
        lower_line = [min(self.ok_bad_lower_threshold, rmsz) for rmsz in rmsz_values_lower]

        # fill below the lower line
        plt.fill_between(
            n_values,
            0,
            lower_line,
            color=COLORS_MATPLOTLIB[CATEGORY_BAD],
            alpha=opacity,
        )

        plt.ylim(0, self.ok_bad_upper_threshold + 0.99)
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

        fig, axs = plt.subplots(2, 2, figsize=(9, 9))
        fig.suptitle(
            f"Comparison of two Normal Distributions with RMSZ 1.0 and {self.ok_bad_upper_threshold} and {self.ok_bad_lower_threshold} respectively"
        )

        # common plot parameters
        xlim = (-7, 7)
        ylim = (0, 0.9)
        opacity = 0.3

        x = np.linspace(-7, 7, 1_000_000)
        gaussian_ref = gaussian(x, self.reference_mean, self.reference_std_dev)

        threshold_specs = [
            (0, self.ok_bad_upper_threshold),
            (1, self.ok_bad_lower_threshold),
        ]

        for col, threshold in threshold_specs:
            gaussian_threshold = gaussian(x, self.reference_mean, threshold)
            inters_x = np.minimum(gaussian_ref, gaussian_threshold)

            # Full-distribution overlap
            axs[0, col].plot(x, gaussian_ref, color="b", label="RMSZ 1.0")
            axs[0, col].fill_between(x, gaussian_ref, 0, color="b", alpha=opacity)
            axs[0, col].plot(x, gaussian_threshold, color="orange", label=f"RMSZ {threshold}")
            axs[0, col].fill_between(x, gaussian_threshold, 0, color="orange", alpha=opacity)
            axs[0, col].plot(x, inters_x, color="r")
            axs[0, col].fill_between(
                x,
                inters_x,
                0,
                facecolor="none",
                edgecolor="r",
                hatch="xx",
                label="intersection",
            )
            axs[0, col].set_xlim(xlim)
            axs[0, col].set_ylim(ylim)

            area_inters_x = np.trapezoid(inters_x, x)
            handles, labels = axs[0, col].get_legend_handles_labels()
            labels[2] += f": {area_inters_x * 100:.1f} %"
            axs[0, col].legend(handles, labels)

            # In-range overlap within +/- threshold
            r = 2
            x2 = np.array([xi for xi in x if abs(xi) <= r])
            gaussian_ref_x2 = gaussian(x2, self.reference_mean, self.reference_std_dev)
            gaussian_threshold_x2 = gaussian(x2, self.reference_mean, threshold)

            axs[1, col].plot(x, gaussian_ref, color="b", label="RMSZ 1.0")
            axs[1, col].fill_between(x2, gaussian_ref_x2, 0, color="b", alpha=opacity, label="RMSZ 1.0")
            axs[1, col].plot(x, gaussian_threshold, color="orange", label=f"RMSZ {threshold}")
            axs[1, col].fill_between(
                x2,
                gaussian_threshold_x2,
                0,
                color="orange",
                alpha=opacity,
                label=f"RMSZ {threshold}",
            )
            for i in [-r, r]:
                axs[1, col].axvline(
                    x=i,
                    linestyle="dashed",
                    color="black",
                )
            axs[1, col].set_xlim(xlim)
            axs[1, col].set_ylim(ylim)

            area_inters_1 = np.trapezoid(gaussian_ref_x2, x2)
            area_inters_2 = np.trapezoid(gaussian_threshold_x2, x2)

            handles, labels = axs[1, col].get_legend_handles_labels()
            labels[1] += f": {area_inters_1 * 100:.1f} %"
            labels[3] += f": {area_inters_2 * 100:.1f} %"
            axs[1, col].legend(handles, labels)
        # axs[1].title('RMSZ Distribution Comparison')
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path)
        else:
            plt.show()
