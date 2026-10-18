"""Accuracy and speed benchmark for TabulatedReionization.

Parts (select with --parts):

nodes   Interpolation error versus node density. Source-time sampling is held fixed and dense
        (large min_steps), so differences come from the table, not from sampling.
speed   The native tanh history tabulated at several node densities: run time and spectra
        relative to the native model, with per-interval sampling capped at the old fixed count
        (source_min_dlnz -> 0) and with the default width-aware sampling.
coarse  A sparse ten-node table: default sampling versus much finer sampling.
tau     get_reionization_optical_depth() versus an independent density integral.

Spectra are lensed totals to --lmax. TT and EE differences are fractional, TE uses
dTE/sqrt(TT*EE), and "max/CV" is the largest difference in units of the per-multipole
cosmic-variance fraction sqrt(2/(2l+1)). Timings are best-of-three after a discarded call.
"""

import argparse
import time

import numpy as np

import camb
from camb import constants
from camb.reionization import TabulatedReionization, TanhReionization


def make_params(reion, lmax):
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.12, mnu=0.06, YHe=0.245)
    p.InitPower.set_params(As=2.1e-9, ns=0.965)
    p.Reion = reion
    p.set_for_lmax(lmax, lens_potential_accuracy=1)
    return p


def spectra(reion, lmax, repeat=1):
    best = np.inf
    for _ in range(repeat):
        start = time.perf_counter()
        results = camb.get_results(make_params(reion, lmax))
        best = min(best, time.perf_counter() - start)
    return results, results.get_cmb_power_spectra(CMB_unit="muK", raw_cl=True)["total"][: lmax + 1], best


def differences(cl, ref):
    ell = np.arange(2, ref.shape[0])
    tt, ee, te = ref[2:, 0], ref[2:, 1], ref[2:, 3]
    cosmic = np.sqrt(2.0 / (2 * ell + 1))
    d_tt = (cl[2:, 0] - tt) / tt
    d_ee = (cl[2:, 1] - ee) / ee
    d_te = (cl[2:, 3] - te) / np.sqrt(tt * ee)
    return {
        "TT": np.abs(d_tt).max(),
        "EE": np.abs(d_ee).max(),
        "TE": np.abs(d_te).max(),
        "max/CV": max(np.abs(d / cosmic).max() for d in (d_tt, d_ee, d_te)),
    }


def text(values):
    return "  ".join(f"{key}={value:.2e}" for key, value in values.items())


def q_tanh(z):
    return 0.5 * (1 - np.tanh((np.asarray(z) - 8.0) / 0.7))


def q_asymmetric(z):
    # Ends at z=6 with a long high-redshift tail; only C1 at z=6, like many simulated histories.
    return np.exp(-np.log(2) * (np.clip(np.asarray(z) - 6.0, 0, None) / 2.5) ** 2)


def log_grid(step, z_max=30.0):
    z = [0.0]
    while z[-1] < z_max:
        z.append((1 + z[-1]) * step - 1)
    z[-1] = z_max
    return np.array(z)


def table(z, q, **options):
    reion = TabulatedReionization().set_table(z, q, history_mode="ionized_fraction")
    for name, value in options.items():
        setattr(reion, name, value)
    return reion


def part_nodes(lmax):
    print("\n== nodes: interpolation error versus node density (fixed dense source sampling)")
    fixed = {"min_steps": 4000, "source_samples_per_interval": 2, "source_min_dlnz": 1.0}
    z_check = np.linspace(0, 30, 6001)
    for name, shape in (("tanh", q_tanh), ("asymmetric", q_asymmetric)):
        z_ref = log_grid(1.002)
        ref_results, ref_cl, _ = spectra(table(z_ref, shape(z_ref), **fixed), lmax)
        ref_tau = ref_results.get_reionization_optical_depth()
        ref_xe = ref_results.get_background_redshift_evolution(z_check, ["x_e"])["x_e"]
        print(f"[{name}] reference: {len(z_ref)} nodes, tau={ref_tau:.6f}")
        for step in (1.01, 1.02, 1.03, 1.05, 1.10, 1.20):
            z = log_grid(step)
            results, cl, _ = spectra(table(z, shape(z), **fixed), lmax)
            xe = results.get_background_redshift_evolution(z_check, ["x_e"])["x_e"]
            summary = {
                "max|dxe|": np.abs(xe - ref_xe).max(),
                "|dtau|": abs(results.get_reionization_optical_depth() - ref_tau),
                **differences(cl, ref_cl),
            }
            print(f"  ln(1+z) step {step - 1:.2f} ({len(z):4d} nodes, dz={9 * (step - 1):.2f} at z=8)  {text(summary)}")


def part_speed(lmax):
    print("\n== speed: native tanh history as a table, old fixed versus width-aware sampling")
    native = TanhReionization()
    native.set_tau(0.0576)
    spectra(native, lmax)  # discard first call (Bessel set-up)
    _, native_cl, native_time = spectra(native, lmax, repeat=3)
    background = camb.get_background(make_params(native, lmax))

    def history(z):
        return background.get_background_redshift_evolution(z, vars=["x_e"])["x_e"]

    print(f"native tanh: {native_time:.2f} s")
    for atol, max_dz in ((1e-5, 0.5), (1e-6, 0.1), (1e-6, 0.04)):
        for label, min_dlnz in (("old: 16 per interval", 1e-12), ("width-aware default", None)):
            reion = TabulatedReionization()
            info = reion.set_function(
                history,
                0,
                30,
                history_mode="total_xe",
                transition_width=1.0,
                include_helium_fullreion=False,
                atol=atol,
                rtol=0,
                max_dz=max_dz,
            )
            if min_dlnz is not None:
                reion.source_min_dlnz = min_dlnz
            _, cl, seconds = spectra(reion, lmax, repeat=3)
            print(
                f"  {info['n_nodes']:4d} nodes, {label:22s} {seconds:5.2f} s ({seconds / native_time:4.1f}x native)  "
                f"vs native: {text(differences(cl, native_cl))}"
            )


def part_coarse(lmax):
    print("\n== coarse: sparse ten-node table, default versus much finer source sampling")
    z = np.array([0, 5, 8, 12, 40, 50, 55, 60, 70, 80.0])
    q = np.array([1, 1, 0.5, 0, 0, 0.05, 0.15, 0.05, 0, 0])
    _, ref, _ = spectra(table(z, q, source_samples_per_interval=128, source_min_dlnz=1e-4), lmax)
    for samples in (4, 8, 16):
        _, cl, _ = spectra(table(z, q, source_samples_per_interval=samples), lmax)
        print(f"  source_samples_per_interval={samples:2d}: {text(differences(cl, ref))}")


def part_tau(lmax):
    print("\n== tau: get_reionization_optical_depth versus independent density integral")
    baseline = TabulatedReionization().set_table(
        [0, 1], [0, 0], history_mode="ionized_fraction", include_helium_fullreion=False
    )
    residual = camb.get_background(make_params(baseline, lmax))
    z = np.linspace(0, 40, 80001)
    for name, shape in (("tanh", q_tanh), ("asymmetric", q_asymmetric)):
        nodes = log_grid(1.01)
        results = camb.get_background(make_params(table(nodes, shape(nodes)), lmax))
        excess = (
            results.get_background_redshift_evolution(z, ["x_e"])["x_e"]
            - residual.get_background_redshift_evolution(z, ["x_e"])["x_e"]
        )
        p = results.Params
        n_h0 = (1 - p.YHe) * 3 * (100000 / constants.Mpc) ** 2 / (8 * np.pi * constants.G) * p.ombh2 / constants.m_H
        weight = (
            constants.c
            * constants.sigma_thomson
            * n_h0
            * (1 + z) ** 2
            / (results.hubble_parameter(z) * 1000 / constants.Mpc)
        )
        independent = np.trapezoid(weight * excess, z)
        tau = results.get_reionization_optical_depth()
        print(
            f"  {name:10s} accessor={tau:.7f} independent={independent:.7f} relative={(tau - independent) / independent:.1e}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--parts", default="tau,nodes,speed,coarse")
    parser.add_argument("--lmax", type=int, default=2500)
    args = parser.parse_args()
    print("camb", camb.__version__, camb.__file__)
    for part in args.parts.split(","):
        {"nodes": part_nodes, "speed": part_speed, "coarse": part_coarse, "tau": part_tau}[part](args.lmax)


if __name__ == "__main__":
    main()
