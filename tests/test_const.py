import numpy as np
import pytest
import pandas as pd
from ligandqc.const import *
from random import random


def test_is_rsr_threshold_boundaries():
	assert is_rsr_good(THRESH_RSR_GOOD)
	assert is_rsr_good(THRESH_RSR_GOOD - 1e-12)
	assert not is_rsr_good(THRESH_RSR_GOOD + 1e-12)

	assert not is_rsr_bad(THRESH_RSR_GOOD)
	assert is_rsr_bad(THRESH_RSR_GOOD + 1e-12)


def test_is_rscc_threshold_boundaries():
	assert is_rscc_good(THRESH_RSCC_GOOD)
	assert is_rscc_good(THRESH_RSCC_GOOD + 1e-12)
	assert not is_rscc_good(THRESH_RSCC_GOOD - 1e-12)

	assert is_rscc_ok(THRESH_RSCC_BAD)
	assert is_rscc_ok((THRESH_RSCC_BAD + THRESH_RSCC_GOOD) / 2)
	assert is_rscc_ok(THRESH_RSCC_GOOD - 1e-12)
	assert not is_rscc_ok(THRESH_RSCC_GOOD)
	assert not is_rscc_ok(THRESH_RSCC_BAD - 1e-12)

	assert is_rscc_bad(THRESH_RSCC_BAD - 1e-12)
	assert not is_rscc_bad(THRESH_RSCC_BAD)


def test_is_structural_outliers_threshold_boundaries():
	assert is_structural_outliers_good(THRESH_STRUCTURAL_OUTLIERS_GOOD)
	assert is_structural_outliers_good(THRESH_STRUCTURAL_OUTLIERS_GOOD - 1)
	assert not is_structural_outliers_good(THRESH_STRUCTURAL_OUTLIERS_GOOD + 1)

	assert not is_structural_outliers_ok(THRESH_STRUCTURAL_OUTLIERS_GOOD)
	assert is_structural_outliers_ok(THRESH_STRUCTURAL_OUTLIERS_GOOD + 1)
	assert is_structural_outliers_ok(THRESH_STRUCTURAL_OUTLIERS_BAD)
	assert not is_structural_outliers_ok(THRESH_STRUCTURAL_OUTLIERS_BAD + 1)

	assert is_structural_outliers_bad(THRESH_STRUCTURAL_OUTLIERS_BAD + 1)
	assert not is_structural_outliers_bad(THRESH_STRUCTURAL_OUTLIERS_BAD)


def test_is_avg_occupancy_threshold_boundaries():
	assert is_avg_occupancy_good(THRESH_AVG_OCCUPANCY_GOOD)
	assert is_avg_occupancy_good(THRESH_AVG_OCCUPANCY_GOOD + 1e-12)
	assert not is_avg_occupancy_good(THRESH_AVG_OCCUPANCY_GOOD - 1e-12)

	assert is_avg_occupancy_ok(THRESH_AVG_OCCUPANCY_BAD)
	assert is_avg_occupancy_ok((THRESH_AVG_OCCUPANCY_BAD + THRESH_AVG_OCCUPANCY_GOOD) / 2)
	assert is_avg_occupancy_ok(THRESH_AVG_OCCUPANCY_GOOD - 1e-12)
	assert not is_avg_occupancy_ok(THRESH_AVG_OCCUPANCY_GOOD)
	assert not is_avg_occupancy_ok(THRESH_AVG_OCCUPANCY_BAD - 1e-12)

	assert is_avg_occupancy_bad(THRESH_AVG_OCCUPANCY_BAD - 1e-12)
	assert not is_avg_occupancy_bad(THRESH_AVG_OCCUPANCY_BAD)


def test_is_atomic_clashes_threshold_boundaries():
	assert is_atomic_clashes_good(THRESH_ATOMIC_CLASHES_GOOD)
	assert is_atomic_clashes_good(THRESH_ATOMIC_CLASHES_GOOD - 1)
	assert not is_atomic_clashes_good(THRESH_ATOMIC_CLASHES_GOOD + 1)

	assert not is_atomic_clashes_ok(THRESH_ATOMIC_CLASHES_GOOD)
	assert is_atomic_clashes_ok(THRESH_ATOMIC_CLASHES_GOOD + 1)
	assert is_atomic_clashes_ok(THRESH_ATOMIC_CLASHES_BAD)
	assert not is_atomic_clashes_ok(THRESH_ATOMIC_CLASHES_BAD + 1)

	assert is_atomic_clashes_bad(THRESH_ATOMIC_CLASHES_BAD + 1)
	assert not is_atomic_clashes_bad(THRESH_ATOMIC_CLASHES_BAD)


def test_is_stereochemical_errors_threshold_boundaries():
	assert is_stereo_chemical_errors_good(THRESH_STEREO_CHEMICAL_ERRORS_GOOD)
	assert is_stereo_chemical_errors_good(THRESH_STEREO_CHEMICAL_ERRORS_GOOD - 1)
	assert not is_stereo_chemical_errors_good(THRESH_STEREO_CHEMICAL_ERRORS_GOOD + 1)

	assert is_stereo_chemical_errors_bad(THRESH_STEREO_CHEMICAL_ERRORS_BAD + 1)
	assert not is_stereo_chemical_errors_bad(THRESH_STEREO_CHEMICAL_ERRORS_BAD)


@pytest.mark.parametrize(
	"predicate",
	[
		is_rsr_good,
		is_rsr_bad,
		is_rscc_good,
		is_rscc_ok,
		is_rscc_bad,
		is_structural_outliers_good,
		is_structural_outliers_ok,
		is_structural_outliers_bad,
		is_avg_occupancy_good,
		is_avg_occupancy_ok,
		is_avg_occupancy_bad,
		is_atomic_clashes_good,
		is_atomic_clashes_ok,
		is_atomic_clashes_bad,
		is_stereo_chemical_errors_good,
		is_stereo_chemical_errors_bad,
	],
)
def test_all_predicates_return_false_for_nan(predicate):
	assert not predicate(np.nan)




def test_produce_artificial_data():
	pass