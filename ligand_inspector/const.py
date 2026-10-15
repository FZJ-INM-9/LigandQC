
# categories
CATEGORY_GOOD = "good"
CATEGORY_OK = "ok"
CATEGORY_BAD = "bad"
CATEGORY_UNKNOWN = "unknown"
CATEGORY_ORDER = [CATEGORY_GOOD, CATEGORY_OK, CATEGORY_BAD]

# metrics
METRIC_RSR = "Real Space R-Factor"
METRIC_RSCC = "Real Space Correlation Coefficient"
METRIC_STRUCTURAL_OUTLIERS = "Structural Outliers"
METRIC_RMSZ_BOND_LENGTH = "RMSZ Bond Length"
METRIC_RMSZ_BOND_ANGLE = "RMSZ Bond Angle"
METRIC_AVG_OCCUPANCY = "Average Occupancy"
METRIC_ATOMIC_CLASHES = "Atomic Clashes"
METRIC_STEREO_CHEMICAL_ERRORS = "Stereochemical Errors"
METRICS = [
    METRIC_RSR,
    METRIC_RSCC,
    METRIC_STRUCTURAL_OUTLIERS,
    METRIC_RMSZ_BOND_LENGTH,
    METRIC_RMSZ_BOND_ANGLE,
    METRIC_AVG_OCCUPANCY,
    METRIC_ATOMIC_CLASHES,
    METRIC_STEREO_CHEMICAL_ERRORS,
]

METRIC_RSR_SHORT = "RSR"
METRIC_RSCC_SHORT = "RSCC"
METRIC_STRUCTURAL_OUTLIERS_SHORT = "Outliers"
METRIC_RMSZ_BOND_LENGTH_SHORT = "RMSZ Bond Length"
METRIC_RMSZ_BOND_ANGLE_SHORT = "RMSZ Bond Angle"
METRIC_AVG_OCCUPANCY_SHORT = "Avg Occupancy"
METRIC_ATOMIC_CLASHES_SHORT = "Clashes"
METRIC_STEREO_CHEMICAL_ERRORS_SHORT = "Stereo Errors"
METRICS_SHORT = [
    METRIC_RSR_SHORT,
    METRIC_RSCC_SHORT,
    METRIC_STRUCTURAL_OUTLIERS_SHORT,
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
is_rsr_good = lambda rsr: rsr <= THRESH_RSR_GOOD
is_rsr_bad = lambda rsr: rsr > THRESH_RSR_BAD

# Real Space Correlation Coefficient (RSCC)
THRESH_RSCC_GOOD = 0.95
THRESH_RSCC_BAD = 0.8
is_rscc_good = lambda rscc: rscc >= THRESH_RSCC_GOOD
is_rscc_ok = lambda rscc: THRESH_RSCC_BAD <= rscc < THRESH_RSCC_GOOD
is_rscc_bad = lambda rscc: rscc < THRESH_RSCC_BAD

# RMSZ-bond-length and RMSZ-bond-angle
# good and ok thresholds are not fix, will be determined in the RMSZInterpreter class.
THRESH_RMSZ_BAD = 2.0
is_rmsz_bad = lambda rmsz: rmsz > THRESH_RMSZ_BAD

# Structural Outliers
THRESH_STRUCTURAL_OUTLIERS_GOOD = 5
THRESH_STRUCTURAL_OUTLIERS_BAD = 15
is_structural_outliers_good = lambda outliers: outliers <= THRESH_STRUCTURAL_OUTLIERS_GOOD
is_structural_outliers_ok = lambda outliers: THRESH_STRUCTURAL_OUTLIERS_GOOD < outliers <= THRESH_STRUCTURAL_OUTLIERS_BAD
is_structural_outliers_bad = lambda outliers: outliers > THRESH_STRUCTURAL_OUTLIERS_BAD

# Average occupancy
THRESH_AVG_OCCUPANCY_GOOD = 0.95
THRESH_AVG_OCCUPANCY_BAD = 0.5
is_avg_occupancy_good = lambda occupancy: occupancy >= THRESH_AVG_OCCUPANCY_GOOD
is_avg_occupancy_ok = lambda occupancy: THRESH_AVG_OCCUPANCY_BAD <= occupancy < THRESH_AVG_OCCUPANCY_GOOD
is_avg_occupancy_bad = lambda occupancy: occupancy < THRESH_AVG_OCCUPANCY_BAD

# Atomic clashes
THRESH_ATOMIC_CLASHES_GOOD = 0
THRESH_ATOMIC_CLASHES_BAD = 2
is_atomic_clashes_good = lambda clashes: clashes <= THRESH_ATOMIC_CLASHES_GOOD
is_atomic_clashes_ok = lambda clashes: THRESH_ATOMIC_CLASHES_GOOD < clashes <= THRESH_ATOMIC_CLASHES_BAD
is_atomic_clashes_bad = lambda clashes: clashes > THRESH_ATOMIC_CLASHES_BAD

# Stereo-chemical errors
THRESH_STEREO_CHEMICAL_ERRORS_GOOD = 0
THRESH_STEREO_CHEMICAL_ERRORS_BAD = 0
is_stereo_chemical_errors_good = lambda errors: errors <= THRESH_STEREO_CHEMICAL_ERRORS_GOOD
is_stereo_chemical_errors_bad = lambda errors: errors > THRESH_STEREO_CHEMICAL_ERRORS_BAD

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