"""Tests of the table preparation without requiring the CAMB shared library."""

import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np
from scipy.interpolate import BPoly, CubicSpline

_spec = importlib.util.spec_from_file_location("_tabulated", Path(__file__).resolve().parents[1] / "_tabulated.py")
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
prepare_table = _module.prepare_table
sample_function = _module.sample_function


class TabulatedInterpolationTest(unittest.TestCase):
    def test_sort_and_layout(self):
        z, y, coefficients = prepare_table([2, 0, 1], [0.2, 1, 0.6])
        np.testing.assert_array_equal(z, [0, 1, 2])
        np.testing.assert_array_equal(y, [1, 0.6, 0.2])
        self.assertEqual(coefficients.shape, (6, 2))
        self.assertEqual(coefficients.dtype, np.dtype("float64"))
        self.assertTrue(coefficients.flags.f_contiguous)

    def test_random_bounds_and_c2(self):
        rng = np.random.default_rng(41)
        for _ in range(30):
            z = np.r_[0, np.cumsum(rng.uniform(0.1, 2, 19))]
            y = rng.uniform(0, 1.2, len(z))
            _, _, c = prepare_table(z, y)
            polynomial = BPoly(c, z)
            probes = z[:-1, None] + np.diff(z)[:, None] * np.linspace(0, 1, 31)
            values = polynomial(probes)
            self.assertTrue(np.all(values >= np.minimum(y[:-1], y[1:])[:, None] - 1e-14))
            self.assertTrue(np.all(values <= np.maximum(y[:-1], y[1:])[:, None] + 1e-14))
            for derivative in (0, 1, 2):
                derived = polynomial.derivative(derivative)
                left = derived.c[-1, :-1]
                right = derived.c[0, 1:]
                np.testing.assert_allclose(left, right, atol=1e-10, rtol=1e-10)
            for derivative in (1, 2):
                np.testing.assert_allclose(polynomial.derivative(derivative)(z[[0, -1]]), 0, atol=1e-12)

    def test_irregular_nonmonotonic_control_hulls_and_node_derivatives(self):
        rng = np.random.default_rng(9271)
        epsilon = np.finfo(np.float64).eps
        for trial in range(80):
            z = np.r_[0, np.cumsum(np.exp(rng.uniform(-3, 1, rng.integers(4, 35))))]
            y = rng.uniform(0, 1, len(z))
            y[:3], y[-1] = [0.1, 0.8, 0.2], 1
            if trial % 10 == 0:
                y[-2] = y[-1]
            # Sorting must not alter the constraints or derivatives.
            order = rng.permutation(len(z))
            sorted_z, sorted_y, c = prepare_table(z[order], y[order])
            np.testing.assert_array_equal(sorted_z, z)
            np.testing.assert_array_equal(sorted_y, y)
            lower, upper = np.minimum(y[:-1], y[1:]), np.maximum(y[:-1], y[1:])
            self.assertTrue(np.all(c >= lower[None, :] - 8 * epsilon))
            self.assertTrue(np.all(c <= upper[None, :] + 8 * epsilon))
            direction = np.sign(np.diff(y))
            self.assertTrue(np.all(direction[None, :] * np.diff(c, axis=0) >= -8 * epsilon))
            midpoint = (y[:-1] + y[1:]) / 2
            self.assertTrue(np.all(direction * (midpoint - c[2]) >= -8 * epsilon))
            self.assertTrue(np.all(direction * (c[3] - midpoint) >= -8 * epsilon))
            polynomial = BPoly(c, z)
            probes = z[:-1, None] + np.diff(z)[:, None] * np.linspace(0, 1, 41)
            sampled = polynomial(probes)
            self.assertTrue(np.all(sampled >= lower[:, None] - 8 * epsilon))
            self.assertTrue(np.all(sampled <= upper[:, None] + 8 * epsilon))
            self.assertTrue(np.all(sampled >= -8 * epsilon))
            self.assertTrue(np.all(sampled <= 1 + 8 * epsilon))
            for derivative in (1, 2):
                derived = polynomial.derivative(derivative)
                left, right = derived.c[-1, :-1], derived.c[0, 1:]
                # Bernstein differences lose O(epsilon / dz**derivative) in
                # floating point; compare against that local arithmetic scale.
                arithmetic_scale = np.maximum(np.diff(z)[:-1] ** -derivative, np.diff(z)[1:] ** -derivative)
                normalized_jump = np.abs(left - right) / arithmetic_scale
                self.assertTrue(np.all(normalized_jump < 256 * epsilon))

    def test_plateau_and_two_nodes(self):
        for z, y in (([0, 1], [0, 1]), ([0, 1, 2, 3, 4], [0, 1, 1, 0, 0])):
            _, _, c = prepare_table(z, y)
            polynomial = BPoly(c, z)
            np.testing.assert_allclose(polynomial(z), y)
        np.testing.assert_allclose(polynomial(np.linspace(1, 2, 25)), 1)
        np.testing.assert_allclose(polynomial(np.linspace(3, 4, 25)), 0)

    def test_nonmonotonic_history_stays_locally_bounded(self):
        z = np.array([0, 1, 2, 3])
        y = np.array([0.1, 0.8, 0.2, 1.0])
        _, _, c = prepare_table(z, y)
        polynomial = BPoly(c, z)
        np.testing.assert_allclose(polynomial(z), y)
        for i in range(len(z) - 1):
            values = polynomial(np.linspace(z[i], z[i + 1], 1001))
            self.assertTrue(np.all(values >= min(y[i], y[i + 1]) - 1e-14))
            self.assertTrue(np.all(values <= max(y[i], y[i + 1]) + 1e-14))
            self.assertTrue(np.all(np.sign(y[i + 1] - y[i]) * np.diff(values) >= -1e-14))
        values = polynomial(np.linspace(0, 3, 3001))
        self.assertTrue(np.all((values >= 0) & (values <= 1)))

    def test_quadratic_second_derivative_converges(self):
        errors = []
        for count in (17, 33, 65):
            z = np.linspace(0, 1, count)
            _, _, c = prepare_table(z, 1 + z**2)
            polynomial = BPoly(c, z)
            errors.append(np.max(np.abs(polynomial.derivative(2)(np.linspace(0.4, 0.6, 31)) - 2)))
        self.assertLess(errors[-1], 1e-5)
        self.assertLess(errors[-1], errors[0] / 100)

    def test_continuity_at_extrema_and_limiter_switches(self):
        # Uniformly scaling both derivatives to zero on a sign conflict used
        # to produce a finite jump at an otherwise infinitesimal perturbation.
        # Include an extremum and an incoming Bernstein-increment sign change
        # in a strictly increasing table (at y[1] = 0.6125).
        cases = (([0, 1, 2], [0, 1, 0], 2), ([0, 1, 2, 3], [0, 0.6125, 0.8, 1], 1))
        for z, y, changed_node in cases:
            probes = np.linspace(z[0], z[-1], 1001)
            _, _, c = prepare_table(z, y)
            baseline = BPoly(c, z)(probes)
            for epsilon in (1e-3, 1e-6, 1e-9, 1e-12):
                for direction in (-1, 1):
                    changed = np.array(y, dtype=float)
                    changed[changed_node] += direction * epsilon
                    if changed[changed_node] < 0:
                        continue
                    _, _, perturbed_c = prepare_table(z, changed)
                    change = np.max(np.abs(BPoly(perturbed_c, z)(probes) - baseline))
                    self.assertLessEqual(change, 4 * epsilon + 1e-14)

    def test_tanh_parameter_difference_converges_on_fixed_grids(self):
        """Check a smooth, limiter-inactive neighbourhood, not global differentiability.

        Limiter switching and an adaptive mesh changing with the parameter are
        outside this test. C2 in redshift alone does not promise differentiability
        with respect to every possible input table or physical parameter.
        """
        center, width = 7.6, 0.8
        probes = np.linspace(center - 1, center + 1, 151)
        exact_derivative = 0.5 / width / np.cosh((probes - center) / width) ** 2

        def parameter_difference(z, step):
            interpolated = []
            interior = (z >= probes[0]) & (z <= probes[-1])
            for shifted_center in (center - step, center + step):
                y = (1 - np.tanh((z - shifted_center) / width)) / 2
                _, _, c = prepare_table(z, y)
                polynomial = BPoly(c, z)
                cubic = CubicSpline(z, y, bc_type="clamped")
                # Confirm the limiter leaves the probed neighbourhood alone.
                for derivative in (1, 2):
                    np.testing.assert_allclose(
                        polynomial(z[interior], nu=derivative),
                        cubic(z[interior], nu=derivative),
                        atol=1e-9,
                        rtol=1e-9,
                    )
                interpolated.append(polynomial(probes))
            return (interpolated[1] - interpolated[0]) / (2 * step)

        fine_grid = np.linspace(0, 16, 321)
        step_errors = [
            np.max(np.abs(parameter_difference(fine_grid, step) - exact_derivative))
            for step in (0.08, 0.04, 0.02, 0.01)
        ]
        # Second-order central differencing before the interpolation error floor.
        self.assertTrue(np.all(np.array(step_errors[:-1]) / step_errors[1:] > 3.5))
        self.assertLess(step_errors[-1], 3.5e-5)
        grid_errors = [
            np.max(np.abs(parameter_difference(np.linspace(0, 16, count), 1e-4) - exact_derivative))
            for count in (81, 161, 321)
        ]
        self.assertTrue(np.all(np.array(grid_errors[:-1]) / grid_errors[1:] > 8))
        self.assertLess(grid_errors[-1], 5e-7)

    def test_bad_tables(self):
        for z, y in (
            ([], []),
            ([0], [0]),
            ([0, 0], [0, 1]),
            ([0, 1], [1]),
            ([-1, 0], [0, 1]),
            ([0, np.nan], [0, 1]),
            ([0, 1], [0, -1]),
            ([0, 1], [0, np.inf]),
            ([[0, 1]], [[0, 1]]),
            (np.array([0j, 1j]), [0, 1]),
        ):
            with self.subTest(z=z, y=y), self.assertRaises(ValueError):
                prepare_table(z, y)

    def test_roundoff_separated_table_nodes_are_rejected_not_removed(self):
        for z in ([0, 1, np.nextafter(1.0, 2.0), 2], [0, 6.1, 6.1000000000000005, 7]):
            with self.subTest(z=z), self.assertRaisesRegex(ValueError, "too close"):
                prepare_table(z, [0, 0.3, 0.7, 1])
        # A genuinely resolved narrow feature is retained, not coalesced.
        z = np.array([0, 1, 1 + 1e-10, 2])
        sorted_z, y, c = prepare_table(z, [0, 0.3, 0.7, 1])
        np.testing.assert_array_equal(sorted_z, z)
        np.testing.assert_allclose(BPoly(c, sorted_z)(z), y)

    def test_sampling_coalesces_only_roundoff_grid_marker_collisions(self):
        marker = 6.1
        z, _, _, _ = sample_function(lambda z: 0.3, 0, 10, max_dz=0.1, breakpoints=[marker])
        self.assertEqual(np.count_nonzero(np.abs(z - marker) < 1e-12), 1)
        self.assertIn(marker, z)
        markers = [6.1, 6.1 + 1e-10]
        z, _, _, _ = sample_function(lambda z: 0.3, 0, 10, max_dz=0.1, breakpoints=markers)
        for resolved_marker in markers:
            self.assertIn(resolved_marker, z)
        z, _, _, _ = sample_function(
            lambda z: 0.3, 0, 10, breakpoints=[np.nextafter(0.0, 1.0), np.nextafter(10.0, 0.0)]
        )
        np.testing.assert_array_equal(z[[0, -1]], [0, 10])

    def test_functions_dense_error(self):
        for function in (lambda z: (1 - np.tanh((z - 5) / 0.4)) / 2, lambda z: 0.001 + np.exp(-(((z - 5) / 0.3) ** 2))):
            z, _, c, diagnostics = sample_function(function, 0, 10, atol=1e-6, rtol=1e-4, max_dz=0.2)
            probes = np.linspace(0, 10, 20001)
            truth = function(probes)
            ratio = np.abs(BPoly(c, z)(probes) - truth) / (1e-6 + 1e-4 * np.abs(truth))
            self.assertLess(np.max(ratio), 1.1)
            self.assertEqual(diagnostics["n_nodes"], len(z))
            self.assertLessEqual(np.max(np.diff(z)), 0.2 + 1e-14)

    def test_constant_and_scalar_callable(self):
        z, y, c, diagnostics = sample_function(lambda z: 0.2, 0, 3)
        np.testing.assert_allclose(y, 0.2)
        np.testing.assert_allclose(c, 0.2)
        self.assertEqual(diagnostics["iterations"], 1)
        z, _, c, _ = sample_function(lambda z: 1 + math.exp(-z), 0, 3)
        probes = np.linspace(0, 3, 1001)
        np.testing.assert_allclose(BPoly(c, z)(probes), 1 + np.exp(-probes), atol=2e-4)

    def test_narrow_feature_requires_resolution(self):
        def function(z):
            return np.exp(-(((z - 0.137) / 0.0001) ** 2))

        z, _, c, _ = sample_function(function, 0, 1, max_dz=1)
        self.assertLess(BPoly(c, z)(0.137), 1e-6)
        z, _, c, _ = sample_function(function, 0, 1, max_dz=1, breakpoints=[0.137])
        self.assertAlmostEqual(BPoly(c, z)(0.137), 1)
        self.assertIn("not a rigorous bound", sample_function.__doc__)

    def test_bad_sampling_requests(self):
        cases = (
            {"z_min": -1},
            {"z_max": 0},
            {"atol": -1},
            {"rtol": np.nan},
            {"max_dz": 0},
            {"atol": 0, "rtol": 0},
            {"max_nodes": 1},
            {"max_nodes": 2.5},
            {"breakpoints": [2]},
            {"breakpoints": [[0.5]]},
            {"max_nodes": 2, "max_dz": 0.1},
        )
        for case in cases:
            arguments = {"z_min": 0, "z_max": 1, **case}
            with self.subTest(case=case), self.assertRaises(ValueError):
                sample_function(lambda z: z, **arguments)
        for function in (None, lambda z: np.nan, lambda z: -1, lambda z: [1, 2, 3], lambda z: z * 1j):
            with self.subTest(function=function), self.assertRaises((TypeError, ValueError)):
                sample_function(function, 0, 1)
        with self.assertRaisesRegex(ValueError, "budget exhausted"):
            sample_function(lambda z: z, 0, 1, atol=1e-12, rtol=0, max_nodes=3)


if __name__ == "__main__":
    unittest.main()
