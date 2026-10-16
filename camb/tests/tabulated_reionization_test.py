"""Real-backend checks of history semantics and independently resolved spectra."""

import json
import unittest

import numpy as np
from scipy.integrate import simpson

import camb
from camb import constants
from camb.reionization import ExpReionization, TabulatedReionization


def parameters():
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.12, mnu=0.06, YHe=0.245)
    p.InitPower.set_params(As=2.1e-9, ns=0.965)
    p.set_for_lmax(100, lens_potential_accuracy=0)
    p.DoLensing = False
    return p


def integrated_tau(results, zmax=80.0):
    """Total tau over [0,zmax], not the excess over recombination."""
    z = np.linspace(0, zmax, 8001)
    bg = results.get_background_redshift_evolution(z, vars=["x_e", "opacity"])
    h_si = results.hubble_parameter(z) * 1000 / constants.Mpc
    rho_b0 = 3 * (100000 / constants.Mpc) ** 2 / (8 * np.pi * constants.G) * results.Params.ombh2
    n_h0 = (1 - results.Params.YHe) * rho_b0 / constants.m_H
    tau_density = simpson(constants.c * constants.sigma_thomson * n_h0 * (1 + z) ** 2 * bg["x_e"] / h_si, x=z)
    tau_opacity = simpson(bg["opacity"] * (constants.c / 1000) / results.hubble_parameter(z), x=z)
    np.testing.assert_allclose(tau_density, tau_opacity, rtol=2e-6)
    return float(tau_density)


def compare_native(name):
    p = parameters()
    p.set_for_lmax(300, lens_potential_accuracy=0)
    p.Accuracy.AccuracyBoost = 2
    p.Accuracy.BackgroundTimeStepBoost = 2
    p.Accuracy.TimeStepBoost = 2
    if name == "exp":
        p.Reion = ExpReionization()
    p.Reion.set_zrei(8.0)
    calculator = camb.get_results if name == "tanh" else camb.get_background
    reference = calculator(p)
    ref_cl = reference.get_unlensed_scalar_cls(lmax=100, CMB_unit="muK", raw_cl=True) if name == "tanh" else None
    ref_tau = integrated_tau(reference)
    query = np.linspace(0, 80, 4001)
    ref_xe = reference.get_background_redshift_evolution(query, vars=["x_e"])["x_e"]
    records = []

    def native_history(z):
        return reference.get_background_redshift_evolution(z, vars=["x_e"])["x_e"]

    for tolerance in (1e-4, 1e-5, 1e-6):
        test = p.copy()
        test.Reion = TabulatedReionization()
        sampling = test.Reion.set_function(
            native_history,
            0,
            80,
            history_mode="total_xe",
            transition_width=1.0,
            include_helium_fullreion=False,
            atol=tolerance,
            rtol=0,
            max_dz=1.0,
        )
        result = calculator(test)
        xe = result.get_background_redshift_evolution(query, vars=["x_e"])["x_e"]
        records.append(
            {
                "sampling_atol": tolerance,
                "n_nodes": sampling["n_nodes"],
                "max_abs_xe": float(np.max(np.abs(xe - ref_xe))),
                "abs_delta_tau_total_0_80": abs(integrated_tau(result) - ref_tau),
            }
        )
        if ref_cl is not None:
            cl = result.get_unlensed_scalar_cls(lmax=100, CMB_unit="muK", raw_cl=True)
            records[-1].update(spectrum_errors(cl, ref_cl))
    return {"native": name, "tau_total_0_80": ref_tau, "refinements": records}


def spectrum_errors(cl, reference):
    normalization = np.sqrt(reference[2:, 0] * reference[2:, 1])
    return {
        "max_fractional_EE_l2_100": float(np.max(np.abs(cl[2:, 1] / reference[2:, 1] - 1))),
        "max_fractional_TT_l2_100": float(np.max(np.abs(cl[2:, 0] / reference[2:, 0] - 1))),
        "max_delta_TE_over_sqrt_TT_EE": float(np.max(np.abs(cl[2:, 3] - reference[2:, 3]) / normalization)),
    }


class TabulatedReionizationTest(unittest.TestCase):
    def test_zero_fraction_recovers_recombination(self):
        p = parameters()
        p.Reion.Reionization = False
        reference = camb.get_background(p)
        p.Reion = TabulatedReionization().set_table(
            [0, 20, 80], [0, 0, 0], history_mode="ionized_fraction", include_helium_fullreion=False
        )
        result = camb.get_background(p)
        z = np.linspace(0, 100, 101)
        np.testing.assert_allclose(
            result.get_background_redshift_evolution(z, vars=["x_e"])["x_e"],
            reference.get_background_redshift_evolution(z, vars=["x_e"])["x_e"],
            rtol=2e-6,
            atol=1e-9,
        )

    def test_full_fraction_and_second_helium_are_not_double_counted(self):
        p = parameters()
        z = np.array([0.0, 2.0, 3.5, 5.0])
        curves = []
        for helium in (False, True):
            p.Reion = TabulatedReionization().set_table(
                [0, 6, 10, 20],
                [1, 1, 0, 0],
                history_mode="ionized_fraction",
                include_helium_fullreion=helium,
            )
            curves.append(camb.get_background(p).get_background_redshift_evolution(z, vars=["x_e"])["x_e"])
        f_he = p.YHe / (constants.mass_ratio_He_H * (1 - p.YHe))
        np.testing.assert_allclose(curves[0], 1 + f_he, atol=2e-7)
        expected = f_he * (1 + np.tanh((3.5 - z) / 0.4)) / 2
        np.testing.assert_allclose(curves[1] - curves[0], expected, atol=2e-7)

    def test_nonmonotonic_early_peak_survives_above_redshift_50(self):
        p = parameters()
        p.Reion = TabulatedReionization().set_table(
            [0, 5, 8, 12, 40, 50, 55, 60, 70, 80],
            [1, 1, 0.5, 0, 0, 0.05, 0.15, 0.05, 0, 0],
            history_mode="ionized_fraction",
            include_helium_fullreion=False,
        )
        result = camb.get_background(p)
        z = np.array([30.0, 40.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0])
        xe = result.get_background_redshift_evolution(z, vars=["x_e"])["x_e"]
        self.assertTrue(np.all(np.isfinite(xe)))
        self.assertGreater(xe[3], 0.15)
        self.assertLess(xe[1], 1e-3)
        self.assertLess(xe[-1], 1e-3)
        self.assertGreater(integrated_tau(result), 0)

    def test_total_history_can_cross_baseline_without_singular_extra_tau(self):
        p = parameters()
        p.Reion.Reionization = False
        baseline = camb.get_background(p)
        z = np.linspace(0, 80, 321)
        residual = baseline.get_background_redshift_evolution(z, vars=["x_e"])["x_e"]
        values = residual * (1 + 0.25 * np.sin(z))
        p.Reion = TabulatedReionization().set_table(z, values, history_mode="total_xe", transition_width=2)
        result = camb.get_background(p)
        self.assertTrue(np.all(np.isfinite(list(result.get_derived_params().values()))))
        self.assertGreater(integrated_tau(result), 0)

    def test_api_rejects_ambiguous_inversion_and_bad_flags(self):
        model = TabulatedReionization()
        with self.assertRaises(ValueError):
            model.set_tau(0.05)
        with self.assertRaises(ValueError):
            model.set_zrei(8)
        with self.assertRaises(ValueError):
            model.get_zre(parameters())
        with self.assertRaises(TypeError):
            model.set_table([0, 10], [1, 0], include_helium_fullreion="False")
        with self.assertRaises(ValueError):
            model.set_table([0, 10], [1, 0.2])
        with self.assertRaises(ValueError):
            model.set_table([0, 10], [1, 0.2], history_mode="total_xe")
        with self.assertRaises(TypeError):
            model.write_ini(None)

    def test_budget_uses_current_helium_and_invalid_controls_fail_cleanly(self):
        p = parameters()
        p.Reion = TabulatedReionization().set_table([0, 80], [1.2, 0])
        with self.assertRaises(camb.CAMBError):
            camb.get_background(p)
        p.Reion = TabulatedReionization().set_table([0, 80], [1, 0], history_mode="ionized_fraction")
        p.Reion.helium_delta_redshift = float("nan")
        with self.assertRaises(camb.CAMBError):
            camb.get_background(p)
        p.Reion.helium_delta_redshift = 0.4
        p.Reion.source_samples_per_interval = 1
        with self.assertRaises(camb.CAMBError):
            camb.get_background(p)
        p.Reion.source_samples_per_interval = 16
        p.Reion.helium_redshiftstart = 0
        with self.assertRaises(camb.CAMBError):
            camb.get_background(p)
        p.Reion.helium_redshiftstart = 5.5
        p.Reion.include_heating = True
        p.Reion.reion_temperature = float("nan")
        with self.assertRaises(camb.CAMBError):
            camb.get_background(p)

    def test_explicit_connection_preserves_input_and_partial_low_z_endpoint(self):
        p = parameters()
        p.Accuracy.BackgroundTimeStepBoost = 2
        p.Reion.Reionization = False
        baseline = camb.get_background(p)
        p.Reion = TabulatedReionization().set_table(
            [2, 8, 50],
            [0.6, 1, 0.2],
            history_mode="ionized_fraction",
            transition_width=4,
            include_helium_fullreion=False,
        )
        result = camb.get_background(p)
        z = np.array([0.0, 1.0, 2.0, 50.0, 52.0, 54.0, 55.0])
        q = np.array([0.6, 0.6, 0.6, 0.2, 0.1, 0, 0])
        residual = baseline.get_background_redshift_evolution(z, vars=["x_e"])["x_e"]
        actual = result.get_background_redshift_evolution(z, vars=["x_e"])["x_e"]
        f_he = p.YHe / (constants.mass_ratio_He_H * (1 - p.YHe))
        np.testing.assert_allclose(actual, (1 - q) * residual + q * (1 + f_he), atol=1e-5)
        np.testing.assert_array_equal(p.Reion.get_table()["x_e"], [0.6, 1, 0.2])

    def test_function_api_and_mutable_connection_width(self):
        p = parameters()
        p.Reion = TabulatedReionization()
        report = p.Reion.set_function(
            lambda z: 0.5 * (1 - np.tanh((z - 8) / 0.8)),
            0,
            80,
            history_mode="ionized_fraction",
            max_dz=1,
            atol=1e-6,
        )
        self.assertGreater(report["n_nodes"], 80)
        p.Reion.transition_width = 2
        result = camb.get_background(p)
        self.assertEqual(result.Params.Reion.z_start, 82)
        copied = p.copy()
        copied.Reion.set_table([0, 80], [0, 0], history_mode="ionized_fraction")
        self.assertGreater(p.Reion.n_nodes, 2)
        saved = p.Reion.get_table()
        self.assertEqual(len(saved["redshifts"]), p.Reion.n_nodes)
        saved["x_e"][:] = 0
        self.assertGreater(np.max(p.Reion.get_table()["x_e"]), 0.9)

    def test_native_tanh_spectra_and_exponential_background(self):
        for name in ("tanh", "exp"):
            with self.subTest(name=name):
                record = compare_native(name)
                print("TABULATED_GOLDEN", json.dumps(record, sort_keys=True), flush=True)
                coarse, _, fine = record["refinements"]
                self.assertLess(fine["max_abs_xe"], 1e-5)
                self.assertLess(fine["abs_delta_tau_total_0_80"], 1e-6)
                if name == "tanh":
                    self.assertLess(fine["max_fractional_EE_l2_100"], 1e-3)
                    self.assertLess(fine["max_fractional_TT_l2_100"], 1e-3)
                    self.assertLess(fine["max_delta_TE_over_sqrt_TT_EE"], 1e-3)
                    self.assertLess(fine["max_fractional_EE_l2_100"], coarse["max_fractional_EE_l2_100"] + 1e-5)
                # Native Exp cuts a nonzero ~1e-4 tail at its high-z handoff.
                # Its TT is not a converged golden reference: see the validation
                # report, including an independent untouched upstream build.
                # Test its sampled background/total tau, not false TT agreement.

    def test_early_nonmonotonic_spectra_resolve_background_and_sources(self):
        p = parameters()
        p.set_for_lmax(300, lens_potential_accuracy=0)
        p.Accuracy.AccuracyBoost = 2
        p.Accuracy.TimeStepBoost = 2
        p.Reion = TabulatedReionization().set_table(
            [0, 5, 8, 12, 40, 50, 55, 60, 70, 80],
            [1, 1, 0.5, 0, 0, 0.05, 0.15, 0.05, 0, 0],
            history_mode="ionized_fraction",
            include_helium_fullreion=False,
        )
        spectra = []
        taus = []
        for background, samples in ((2, 16), (4, 16), (4, 32)):
            p.Accuracy.BackgroundTimeStepBoost = background
            p.Reion.source_samples_per_interval = samples
            result = camb.get_results(p)
            spectra.append(result.get_unlensed_scalar_cls(lmax=100, CMB_unit="muK", raw_cl=True))
            taus.append(integrated_tau(result))
        report = {
            "background_refinement": spectrum_errors(spectra[1], spectra[0]),
            "source_refinement": spectrum_errors(spectra[2], spectra[1]),
            "tau_total_0_80": taus,
        }
        print("TABULATED_EARLY_CONVERGENCE", json.dumps(report, sort_keys=True), flush=True)
        for refinement in ("background_refinement", "source_refinement"):
            for error in report[refinement].values():
                self.assertLess(error, 1e-3)
        self.assertLess(max(taus) - min(taus), 1e-6)


if __name__ == "__main__":
    unittest.main()
