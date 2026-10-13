import pytest

from ligand_inspector.RMSZInterpreter import RMSZInterpreter


def test_rmsz_threshold_finding():
    """
    Test the find_rmsz_threshold method of the RMSZInterpreter class.
    """
    interpreter = RMSZInterpreter()

    # Test with a sample size of 2
    n = 6
    rmsz_threshold = interpreter.find_rmsz_threshold(n)

    # Check if the threshold is within expected bounds
    assert rmsz_threshold == 1.48799


def test_rationale_plots():
    """
    Test the rational_plot method of the RMSZInterpreter class.
    """
    interpreter = RMSZInterpreter()
    interpreter.plot_rationale_categorization(
        "tests/output_files/rmsz_categorization.png", alpha=5e-2
    )
    interpreter.plot_rationale_ok_threshold("tests/output_files/rmsz_ok_threshold.png")


@pytest.mark.parametrize(
    "rmsz_value, n, expected_interpretation",
    [
        (0.2, 13, "good"),
        (0.1, 60, "good"),
        (1.2, 13, "good"),
        (1.1, 60, "good"),
        (1.4, 13, "ok"),
        (1.7, 60, "ok"),
        (2.0, 13, "ok"),
        (2.0, 60, "ok"),
        (2.5, 13, "bad"),
        (2.5, 60, "bad"),
        (3.0, 13, "bad"),
        (3.0, 60, "bad"),
        (4.0, 13, "bad"),
        (4.0, 60, "bad"),
    ],
)
def test_interpretation(rmsz_value, n, expected_interpretation):
    """
    Test the interpretation of a molecule.
    """
    interpreter = RMSZInterpreter()
    interpretation = interpreter.interpret_rmsz(n, rmsz_value)
    rmsz_threshold = interpreter.find_rmsz_threshold(n)
    assert (
        interpretation == expected_interpretation
    ), f"Expected '{expected_interpretation}', but got '{interpretation}' for RMSZ value {rmsz_value}, with threshold {rmsz_threshold}."
