
import numpy as np
import pytest

from ligandqc.const import *


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


def test_is_bond_length_outliers_threshold_boundaries():
	assert is_outliers_bond_length_good(1, 10)
	assert not is_outliers_bond_length_good(2, 10)
	assert not is_outliers_bond_length_good(1, 9)

	assert not is_outliers_bond_length_ok(1, 10)
	assert is_outliers_bond_length_ok(2, 10)
	assert is_outliers_bond_length_ok(2, 8)
	assert not is_outliers_bond_length_ok(3, 8)

	assert is_outliers_bond_length_bad(3, 8)
	assert not is_outliers_bond_length_bad(2, 8)


def test_is_bond_angle_outliers_threshold_boundaries():
	assert is_outliers_bond_angle_good(1, 10)
	assert not is_outliers_bond_angle_good(2, 10)
	assert not is_outliers_bond_angle_good(1, 9)

	assert not is_outliers_bond_angle_ok(1, 10)
	assert is_outliers_bond_angle_ok(2, 10)
	assert is_outliers_bond_angle_ok(2, 8)
	assert not is_outliers_bond_angle_ok(3, 8)

	assert is_outliers_bond_angle_bad(3, 8)
	assert not is_outliers_bond_angle_bad(2, 8)


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
	assert not is_stereo_chemical_errors_good(THRESH_STEREO_CHEMICAL_ERRORS_GOOD + 1)

	assert is_stereo_chemical_errors_bad(THRESH_STEREO_CHEMICAL_ERRORS_BAD + 1)
	assert not is_stereo_chemical_errors_bad(THRESH_STEREO_CHEMICAL_ERRORS_BAD)


@pytest.mark.parametrize(
	"predicate,args",
	[
		(is_rsr_good, (np.nan,)),
		(is_rsr_bad, (np.nan,)),
		(is_rscc_good, (np.nan,)),
		(is_rscc_ok, (np.nan,)),
		(is_rscc_bad, (np.nan,)),
		(is_outliers_bond_length_good, (np.nan, 10)),
		(is_outliers_bond_length_ok, (np.nan, 10)),
		(is_outliers_bond_length_bad, (np.nan, 10)),
		(is_outliers_bond_angle_good, (np.nan, 10)),
		(is_outliers_bond_angle_ok, (np.nan, 10)),
		(is_outliers_bond_angle_bad, (np.nan, 10)),
		(is_avg_occupancy_good, (np.nan,)),
		(is_avg_occupancy_ok, (np.nan,)),
		(is_avg_occupancy_bad, (np.nan,)),
		(is_atomic_clashes_good, (np.nan,)),
		(is_atomic_clashes_ok, (np.nan,)),
		(is_atomic_clashes_bad, (np.nan,)),
		(is_stereo_chemical_errors_good, (np.nan,)),
		(is_stereo_chemical_errors_bad, (np.nan,)),
	],
)
def test_all_predicates_return_false_for_nan(predicate, args):
	assert not predicate(*args)




def test_produce_artificial_data():
	pass