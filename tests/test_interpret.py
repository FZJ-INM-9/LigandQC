import pytest

from ligandqc.interpret import *


class TestRMSZInterpreter:

    def test_rmsz_threshold_finding(self):
        """
        Test the find_rmsz_threshold method of the RMSZInterpreter class.
        """
        interpreter = RMSZInterpreter()

        # Test with a sample size of 2
        n = 6
        rmsz_threshold = interpreter.find_rmsz_threshold(n, tol=1e-6)

        # Check if the threshold is within expected bounds
        assert rmsz_threshold == 1.35834

    def test_f_test(self):
        """
        Test the perform_f_test method of the RMSZInterpreter class.
        """
        interpreter = RMSZInterpreter()

        # Test with a sample size of 2 and a specific RMSZ value
        n = 50
        rmsz_value = 0.9
        alpha = 0.05
        result = interpreter.perform_f_test(rmsz_value, n, alpha)

        # Check if the result is a dictionary with expected keys
        assert isinstance(result, dict)
        assert "F-statistic" in result
        assert "p-value" in result
        assert "significant" in result


    def test_rationale_plots(self):
        """
        Test the plot_rationale_categorization method of the RMSZInterpreter class.
        """
        interpreter = RMSZInterpreter()
        interpreter.plot_rationale_categorization(
            "tests/output_files/rmsz_categorization.png", alpha=0.05
        )
        interpreter.plot_rationale_ok_threshold("tests/output_files/rmsz_thresholds_single_sided.png")


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
    def test_interpretation(self, rmsz_value, n, expected_interpretation):
        """
        Test the interpretation of a molecule.
        """
        interpreter = RMSZInterpreter()
        interpretation = interpreter.interpret_rmsz(n, rmsz_value)
        rmsz_threshold = interpreter.find_rmsz_threshold(n)
        assert (
            interpretation == expected_interpretation
        ), f"Expected '{expected_interpretation}', but got '{interpretation}' for RMSZ value {rmsz_value}, with threshold {rmsz_threshold}."


class TestRMSZInterpreterDoubleSided:

    def test_rmsz_threshold_finding(self):
        """
        Test the find_rmsz_threshold method of the RMSZInterpreter class.
        """
        interpreter = RMSZInterpreterDoubleSided()

        # Test with a sample size of 2
        n = 6
        rmsz_threshold = interpreter.find_rmsz_threshold(n, tol=1e-6)

        # Check if the threshold is within expected bounds
        assert rmsz_threshold == (0.43694, 1.35834)

    def test_f_test(self):
        """
        Test the perform_f_test method of the RMSZInterpreter class.
        """
        interpreter = RMSZInterpreterDoubleSided()

        # Test with a sample size of 2 and a specific RMSZ value
        n = 50
        rmsz_value = 0.9
        alpha = 0.05
        result = interpreter.perform_f_test(rmsz_value, n, alpha)

        # Check if the result is a dictionary with expected keys
        assert isinstance(result, dict)
        assert "F-statistic" in result
        assert "p-value" in result
        assert "significant" in result


    def test_rationale_plots(self):
        """
        Test the plot_rationale_categorization method of the RMSZInterpreter class.
        """
        interpreter = RMSZInterpreterDoubleSided()
        interpreter.plot_rationale_categorization(
            "tests/output_files/rmsz_categorization_double_sided.png", alpha=0.05
        )
        interpreter.plot_rationale_ok_threshold("tests/output_files/rmsz_thresholds_double_sided.png")


    @pytest.mark.parametrize(
        "rmsz_value, n, expected_interpretation",
        [
            (0.25, 3, "good"),
            (0.25, 25, "bad"),
            (0.6, 5, "good"),
            (0.6, 99, "ok"),
            (0.9, 25, "good"),
            (0.9, 99, "good"),
            (0.9, 200, "ok"),
            (1.1, 99, "good"),
            (1.1, 200, "ok"),
            (1.4, 25, "ok"),
            (1.7, 99, "ok"),
            (2.0, 25, "bad"),
            (2.0, 99, "bad"),
            (2.5, 25, "bad"),
            (2.5, 99, "bad"),
            (3.0, 25, "bad"),
            (3.0, 99, "bad"),
            (4.0, 25, "bad"),
            (4.0, 99, "bad"),
        ],
    )
    def test_interpretation(self, rmsz_value, n, expected_interpretation):
        """
        Test the interpretation of a molecule.
        """
        interpreter = RMSZInterpreterDoubleSided()
        interpretation = interpreter.interpret_rmsz(n, rmsz_value)
        rmsz_threshold = interpreter.find_rmsz_threshold(n)
        assert (
            interpretation == expected_interpretation
        ), f"Expected '{expected_interpretation}', but got '{interpretation}' for RMSZ value {rmsz_value}, with threshold {rmsz_threshold}."

