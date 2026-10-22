"""
const
=====

This module contains constants used in the ligand inspector package, including categories, metrics, thresholds, and plotting colors. In addition, it provides functions to evaluate the quality of ligands based on these thresholds.


"""

# categories
CATEGORY_GOOD = "good"
CATEGORY_OK = "ok"
CATEGORY_BAD = "bad"
CATEGORY_UNKNOWN = "unknown"
CATEGORY_ORDER = [CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD]

# metrics
METRIC_RSR = "Real Space R-Factor"
METRIC_RSCC = "Real Space Correlation Coefficient"
METRIC_OUTLIERS_BOND_LENGTH = "Outliers of Bond Length"
METRIC_OUTLIERS_BOND_ANGLE = "Outliers of Bond Angle"
METRIC_STRUCTURAL_OUTLIERS = "Structural Outliers"  # deprecated: kept for backwards compatibility with legacy datasets
METRIC_RMSZ_BOND_LENGTH = "RMSZ Bond Length"
METRIC_RMSZ_BOND_ANGLE = "RMSZ Bond Angle"
METRIC_AVG_OCCUPANCY = "Average Occupancy"
METRIC_ATOMIC_CLASHES = "Atomic Clashes"
METRIC_STEREO_CHEMICAL_ERRORS = "Stereochemical Errors"
METRICS = [
    METRIC_RSR,
    METRIC_RSCC,
    METRIC_OUTLIERS_BOND_LENGTH,
    METRIC_OUTLIERS_BOND_ANGLE,
    METRIC_RMSZ_BOND_LENGTH,
    METRIC_RMSZ_BOND_ANGLE,
    METRIC_AVG_OCCUPANCY,
    METRIC_ATOMIC_CLASHES,
    METRIC_STEREO_CHEMICAL_ERRORS,
]

METRIC_RSR_SHORT = "RSR"
METRIC_RSCC_SHORT = "RSCC"
METRIC_OUTLIERS_BOND_LENGTH_SHORT = "Bond Length Outliers"
METRIC_OUTLIERS_BOND_ANGLE_SHORT = "Bond Angle Outliers"
METRIC_STRUCTURAL_OUTLIERS_SHORT = "Outliers"
METRIC_RMSZ_BOND_LENGTH_SHORT = "RMSZ Bond Length"
METRIC_RMSZ_BOND_ANGLE_SHORT = "RMSZ Bond Angle"
METRIC_AVG_OCCUPANCY_SHORT = "Avg Occupancy"
METRIC_ATOMIC_CLASHES_SHORT = "Clashes"
METRIC_STEREO_CHEMICAL_ERRORS_SHORT = "Stereo Errors"
METRICS_SHORT = [
    METRIC_RSR_SHORT,
    METRIC_RSCC_SHORT,
    METRIC_OUTLIERS_BOND_LENGTH_SHORT,
    METRIC_OUTLIERS_BOND_ANGLE_SHORT,
    METRIC_RMSZ_BOND_LENGTH_SHORT,
    METRIC_RMSZ_BOND_ANGLE_SHORT,
    METRIC_AVG_OCCUPANCY_SHORT,
    METRIC_ATOMIC_CLASHES_SHORT,
    METRIC_STEREO_CHEMICAL_ERRORS_SHORT,
]

NUM_METRICS = len(METRICS)

EVAL_TERM = "Evaluation"
EVAL_TERM_PLURAL = "Evaluations"
EVAL_SUFFIX = f" {EVAL_TERM}"
EVAL_COLUMNS = [f"{metric}{EVAL_SUFFIX}" for metric in METRICS]


# thresholds
#
# Real Space R-factor (RSR)
THRESH_RSR_GOOD = 0.4
THRESH_RSR_BAD = 0.4


def is_rsr_good(rsr):
    return rsr <= THRESH_RSR_GOOD


def is_rsr_ok(rsr):
    return False  # There is no "ok" category for RSR


def is_rsr_bad(rsr):
    return rsr > THRESH_RSR_BAD


# Real Space Correlation Coefficient (RSCC)
THRESH_RSCC_GOOD = 0.95
THRESH_RSCC_BAD = 0.8


def is_rscc_good(rscc):
    return rscc >= THRESH_RSCC_GOOD


def is_rscc_ok(rscc):
    return (rscc >= THRESH_RSCC_BAD) & (rscc < THRESH_RSCC_GOOD)


def is_rscc_bad(rscc):
    return rscc < THRESH_RSCC_BAD


# RMSZ-bond-length and RMSZ-bond-angle
# good and ok thresholds are not fix, will be determined in the RMSZInterpreter class.
THRESH_RMSZ_BAD_UPPER = 2.0
THRESH_RMSZ_BAD_LOWER = 0.5

# Outliers of bond length / bond angle as percentage of all bond/angle entries.
THRESH_OUTLIERS_BOND_LENGTH_GOOD = 0.10
THRESH_OUTLIERS_BOND_LENGTH_BAD = 0.25
THRESH_OUTLIERS_BOND_ANGLE_GOOD = 0.10
THRESH_OUTLIERS_BOND_ANGLE_BAD = 0.25


def calculate_outlier_percentage(outliers, total):
    if outliers is None or total in (None, 0):
        return None
    return outliers / total


def is_outliers_bond_length_good(outliers, total):
    proportion = calculate_outlier_percentage(outliers, total)
    return proportion is not None and proportion <= THRESH_OUTLIERS_BOND_LENGTH_GOOD


def is_outliers_bond_length_ok(outliers, total):
    proportion = calculate_outlier_percentage(outliers, total)
    return (
        proportion is not None
        and proportion > THRESH_OUTLIERS_BOND_LENGTH_GOOD
        and proportion <= THRESH_OUTLIERS_BOND_LENGTH_BAD
    )


def is_outliers_bond_length_bad(outliers, total):
    proportion = calculate_outlier_percentage(outliers, total)
    return proportion is not None and proportion > THRESH_OUTLIERS_BOND_LENGTH_BAD


def is_outliers_bond_angle_good(outliers, total):
    proportion = calculate_outlier_percentage(outliers, total)
    return proportion is not None and proportion <= THRESH_OUTLIERS_BOND_ANGLE_GOOD


def is_outliers_bond_angle_ok(outliers, total):
    proportion = calculate_outlier_percentage(outliers, total)
    return (
        proportion is not None
        and proportion > THRESH_OUTLIERS_BOND_ANGLE_GOOD
        and proportion <= THRESH_OUTLIERS_BOND_ANGLE_BAD
    )


def is_outliers_bond_angle_bad(outliers, total):
    proportion = calculate_outlier_percentage(outliers, total)
    return proportion is not None and proportion > THRESH_OUTLIERS_BOND_ANGLE_BAD


# Average occupancy
THRESH_AVG_OCCUPANCY_GOOD = 0.95
THRESH_AVG_OCCUPANCY_BAD = 0.5


def is_avg_occupancy_good(occupancy):
    return occupancy >= THRESH_AVG_OCCUPANCY_GOOD


def is_avg_occupancy_ok(occupancy):
    return (occupancy >= THRESH_AVG_OCCUPANCY_BAD) & (occupancy < THRESH_AVG_OCCUPANCY_GOOD)


def is_avg_occupancy_bad(occupancy):
    return occupancy < THRESH_AVG_OCCUPANCY_BAD


# Atomic clashes
THRESH_ATOMIC_CLASHES_GOOD = 0
THRESH_ATOMIC_CLASHES_BAD = 2


def is_atomic_clashes_good(clashes):
    return clashes <= THRESH_ATOMIC_CLASHES_GOOD


def is_atomic_clashes_ok(clashes):
    return (clashes > THRESH_ATOMIC_CLASHES_GOOD) & (clashes <= THRESH_ATOMIC_CLASHES_BAD)


def is_atomic_clashes_bad(clashes):
    return clashes > THRESH_ATOMIC_CLASHES_BAD


# Stereo-chemical errors
THRESH_STEREO_CHEMICAL_ERRORS_GOOD = 0
THRESH_STEREO_CHEMICAL_ERRORS_BAD = 0


def is_stereo_chemical_errors_good(errors):
    return errors == THRESH_STEREO_CHEMICAL_ERRORS_GOOD


def is_stereo_chemical_errors_ok(errors):
    return False  # There is no "ok" category for stereo-chemical errors


def is_stereo_chemical_errors_bad(errors):
    return errors > THRESH_STEREO_CHEMICAL_ERRORS_BAD


EVAL_FUNS = {
    METRIC_RSR: {
        CATEGORY_GOOD: is_rsr_good,
        CATEGORY_OK: is_rsr_ok,
        CATEGORY_BAD: is_rsr_bad,
    },
    METRIC_RSCC: {
        CATEGORY_GOOD: is_rscc_good,
        CATEGORY_OK: is_rscc_ok,
        CATEGORY_BAD: is_rscc_bad,
    },
    METRIC_AVG_OCCUPANCY: {
        CATEGORY_GOOD: is_avg_occupancy_good,
        CATEGORY_OK: is_avg_occupancy_ok,
        CATEGORY_BAD: is_avg_occupancy_bad,
    },
    METRIC_ATOMIC_CLASHES: {
        CATEGORY_GOOD: is_atomic_clashes_good,
        CATEGORY_OK: is_atomic_clashes_ok,
        CATEGORY_BAD: is_atomic_clashes_bad,
    },
    METRIC_STEREO_CHEMICAL_ERRORS: {
        CATEGORY_GOOD: is_stereo_chemical_errors_good,
        CATEGORY_OK: is_stereo_chemical_errors_ok,
        CATEGORY_BAD: is_stereo_chemical_errors_bad,
    },
}

# plotting
COLORS_MATPLOTLIB = {
    CATEGORY_GOOD: "green",
    CATEGORY_OK: "yellow",
    CATEGORY_BAD: "orangered",
}

# rich
COLORS_RICH = {
    CATEGORY_GOOD: "green",
    CATEGORY_OK: "yellow",
    CATEGORY_BAD: "red",
}
