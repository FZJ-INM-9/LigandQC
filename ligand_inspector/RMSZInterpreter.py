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
        self.reference_dataset_size = 1E7
        self.reference_dataset = np.random.normal(self.reference_mean, self.reference_std_dev, int(self.reference_dataset_size))
        self.ok_bad_threshold = 2.0
        self.precomputed_rmsz_thresholds = {}
    
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
            "significant": significant
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
            raise ValueError("Sample size n is too large. Please use a smaller sample size.")
        
        # Perform bisection method to find the rmsz threshold
        rmsz_values = []
        while rmsz_start - rmsz_end > 1E-6:
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

    def plot_rationale(self):
        """
        Plot the rationale for flagging ligands as "good" based on the topology and RMSZ value.
        """
        n_values = range(2, 201)
        
        rmsz_values = [self.find_rmsz_threshold(n) for n in n_values]

        plt.figure(figsize=(12, 6))
        plt.plot(n_values, rmsz_values, label="RMSZ Threshold", color='blue')

        plt.axhline(y=2.0, color='red', linestyle='--', label="RMSZ = 2.0")
        # color between the rmsz values and 3.0
        plt.fill_between(n_values, y1=4.0, y2=2.0, color='red', alpha=0.2, label="bad")
        plt.fill_between(n_values, rmsz_values, 2.0, where=(np.array(rmsz_values) < 2.0), color='yellow', alpha=0.2, label="ok")
        # color between the rmsz values and 0
        plt.fill_between(n_values, rmsz_values, 0, where=(np.array(rmsz_values) > 0), color='green', alpha=0.2, label="good")

        # fill above

        plt.ylim(0, 4)
        plt.xlim(0, 200)
        plt.title("Rationale for Flagging Ligands as 'Good'")
        plt.xlabel("Sample Size (n)")
        plt.ylabel("RMSZ Threshold")
        plt.legend()
        plt.grid()
        plt.show()
    
    def plot_explanation(self):
        mean = 0
        std_devs = [0.5, 2.0, 3.0]
        size = 1_000_000
        reference_std_dev = 1.0
        reference_dataset = generate_dataset(mean, reference_std_dev, size)
        reference_rmsz = calculate_rmsz(reference_dataset, reference_std_dev, mean)
        print(f"Reference RMSZ: {reference_rmsz:.4f}")

        rmsz_values = []
        datasets = []
        for std_dev in std_devs:
            dataset = generate_dataset(mean, std_dev, size)
            rmsz = calculate_rmsz(dataset, reference_std_dev, mean)
            rmsz_values.append(rmsz)
            datasets.append(dataset)
            percentage = calculate_percentage(dataset, reference_std_dev, mean)
            print(f"Dataset with std_dev {std_dev} RMSZ: {rmsz:.4f}")
            print(f"Percentage within 2 std_dev: {percentage:.2f}%")
            # Perform t-test
            t_stat, p_value = ttest_ind(reference_dataset, dataset)

        
        # Plotting
        plt.figure(figsize=(12, 6))
        sns.set(style="whitegrid")
        sns.kdeplot(reference_dataset, label="Reference Dataset", fill=True, alpha=0.5)
        for i, dataset in enumerate(datasets):
            sns.kdeplot(dataset, label=f"mean 0; std {std_devs[i]}", fill=True, alpha=0.5)
        plt.vlines([3.0, -3.0], 0, 1, linestyles="dashed", colors="black")
        plt.xlim([-10, 10])
        plt.ylim([0,0.8])
        plt.title("Density Plot of Datasets")
        plt.xlabel("Value")
        plt.ylabel("Density")
        plt.legend()
        plt.show()
        assert False