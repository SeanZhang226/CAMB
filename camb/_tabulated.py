"""Shape-preserving interpolation and controlled sampling of electron histories.

The Bernstein representation is shared with the Fortran evaluator. Its ordered
control points bound each segment by its two data values; shared node derivatives
give C2 continuity, including constant continuation at the two outer endpoints.
"""

import numbers

import numpy as np
from scipy.interpolate import BPoly, CubicSpline


def _node_resolution(left, right):
    return 4 * np.finfo(np.float64).eps * np.maximum(1, np.maximum(np.abs(left), np.abs(right)))


def _merge_sampling_nodes(nodes, priorities):
    """Coalesce only roundoff-scale sampling collisions, preferring markers."""
    order = np.argsort(nodes, kind="stable")
    kept, kept_priorities = [], []
    for node, priority in zip(nodes[order], priorities[order]):
        if kept and node - kept[-1] <= _node_resolution(kept[-1], node):
            if priority > kept_priorities[-1]:
                kept[-1], kept_priorities[-1] = node, priority
        else:
            kept.append(node)
            kept_priorities.append(priority)
    return np.array(kept, dtype=np.float64)


def prepare_table(redshifts, values):
    """Return sorted nodes, values, and C2 quintic Bernstein coefficients.

    Redshifts and electron abundances must be finite, one-dimensional and
    non-negative, with at least two distinct redshifts. There is deliberately no
    imposed upper redshift or abundance bound: their physical meaning is the
    responsibility of the reionization model. No samples are silently removed.
    Neighbours separated by no more than ``4*eps*max(1, abs(z))`` are rejected
    rather than constructing numerically ill-conditioned derivatives.
    """
    if np.iscomplexobj(redshifts) or np.iscomplexobj(values):
        raise ValueError("redshifts and values must be real numeric arrays")
    try:
        z = np.asarray(redshifts, dtype=np.float64)
        y = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("redshifts and values must be real numeric arrays") from exc
    if z.ndim != 1 or y.ndim != 1 or len(z) != len(y) or len(z) < 2:
        raise ValueError("redshifts and values must be equal-length 1D arrays with at least two nodes")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(y)):
        raise ValueError("redshifts and values must be finite")
    if np.any(z < 0) or np.any(y < 0):
        raise ValueError("redshifts and values must be non-negative")
    order = np.argsort(z)
    z, y = z[order].copy(), y[order].copy()
    h = np.diff(z)
    if np.any(h <= 0):
        raise ValueError("redshifts must be distinct")
    if np.any(h <= _node_resolution(z[:-1], z[1:])):
        raise ValueError(
            "adjacent redshifts are too close to floating-point resolution; provide distinct resolved nodes"
        )

    spline = CubicSpline(z, y, bc_type="clamped")
    d, e = spline(z, 1), spline(z, 2)
    # Constant continuation must match both first and second derivatives.
    d[[0, -1]], e[[0, -1]] = 0, 0
    # The zero-curvature choice must remain feasible at every node. Limiting
    # the shared first derivative to 5/4 of each adjacent secant guarantees
    # this. Clipping continuous bounds avoids a discontinuous all-or-nothing
    # rescaling when one Bernstein increment changes sign.
    lower_d = np.full(len(z), -np.inf)
    upper_d = np.full(len(z), np.inf)
    for i, (width, delta) in enumerate(zip(h, np.diff(y))):
        secant_bound = 1.25 * delta / width
        lower_d[i : i + 2] = np.maximum(lower_d[i : i + 2], min(0, secant_bound))
        upper_d[i : i + 2] = np.minimum(upper_d[i : i + 2], max(0, secant_bound))
    d = np.clip(d, lower_d, upper_d)

    lower_e = np.full(len(z), -np.inf)
    upper_e = np.full(len(z), np.inf)
    for i, (width, delta) in enumerate(zip(h, np.diff(y))):
        direction = 1 if delta >= 0 else -1
        curvature_budget = 10 * abs(delta) / width**2
        for node, outgoing in ((i, True), (i + 1, False)):
            directed_slope = direction * d[node] / width
            if outgoing:
                lower = -4 * directed_slope
                upper = max(0, curvature_budget - 8 * directed_slope)
            else:
                lower = min(0, 8 * directed_slope - curvature_budget)
                upper = 4 * directed_slope
            if direction < 0:
                lower, upper = -upper, -lower
            lower_e[node] = max(lower_e[node], lower)
            upper_e[node] = min(upper_e[node], upper)
    e = np.clip(e, lower_e, upper_e)
    controls = np.empty((6, len(h)), dtype=np.float64, order="F")
    controls[0], controls[5] = y[:-1], y[1:]
    controls[1] = y[:-1] + h * d[:-1] / 5
    controls[2] = y[:-1] + 2 * h * d[:-1] / 5 + h**2 * e[:-1] / 20
    controls[4] = y[1:] - h * d[1:] / 5
    controls[3] = y[1:] - 2 * h * d[1:] / 5 + h**2 * e[1:] / 20
    if not np.all(np.isfinite(controls)):
        raise ValueError("table numerical range is too large to construct finite interpolation coefficients")
    return z, y, controls


def _evaluate_function(function, z):
    """Evaluate either an array-aware or scalar-only callable."""
    try:
        result = function(z)
        if np.iscomplexobj(result):
            raise ValueError("function must return real values")
        values = np.asarray(result, dtype=np.float64)
        if values.ndim == 0:
            values = np.full(z.shape, values, dtype=np.float64)
        if values.shape != z.shape:
            raise ValueError("callable returned a shape different from its input")
    except (TypeError, ValueError):
        try:
            result = [function(float(value)) for value in z]
            if np.iscomplexobj(result):
                raise ValueError("function must return real values")
            values = np.array(result, dtype=np.float64)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("function must return real scalar values or an array matching the input shape") from exc
    if values.shape != z.shape:
        raise ValueError("function must return one value per redshift")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("function values must be finite and non-negative")
    return values


def sample_function(
    function,
    z_min,
    z_max,
    *,
    atol=1e-6,
    rtol=1e-4,
    max_dz=1.0,
    breakpoints=(),
    max_nodes=10000,
):
    """Adaptively sample a callable into the same interpolant used by Fortran.

    Each interval is checked at its quarter, half, and three-quarter points
    against ``atol + rtol * abs(function(z))``. This is a numerical probe test,
    not a rigorous bound for an arbitrary black-box function: unresolved narrow
    features can be missed. Set ``max_dz`` small enough for known physical
    scales and provide ``breakpoints`` at known sharp or narrow features.
    Roundoff-close grid/marker collisions are coalesced, preserving explicit
    markers over generated grid points and the exact interval bounds over both.
    Endpoints have zero first and second interpolated derivatives, so callables
    with non-flat endpoints may require additional refinement near the ends.
    """
    if not callable(function):
        raise TypeError("function must be callable")
    try:
        z_min, z_max, atol, rtol, max_dz = map(float, (z_min, z_max, atol, rtol, max_dz))
    except (TypeError, ValueError) as exc:
        raise ValueError("sampling bounds and tolerances must be real numbers") from exc
    if not np.all(np.isfinite([z_min, z_max, atol, rtol, max_dz])):
        raise ValueError("sampling bounds and tolerances must be finite")
    if not 0 <= z_min < z_max:
        raise ValueError("sampling bounds must satisfy 0 <= z_min < z_max")
    if atol < 0 or rtol < 0 or atol + rtol == 0 or max_dz <= 0:
        raise ValueError("tolerances must be non-negative with at least one positive; max_dz must be positive")
    if isinstance(max_nodes, bool) or not isinstance(max_nodes, numbers.Integral) or max_nodes < 2:
        raise ValueError("max_nodes must be an integer >= 2")
    try:
        markers = np.asarray(breakpoints, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError("breakpoints must be a finite 1D sequence within the sampling bounds") from exc
    if markers.ndim != 1 or not np.all(np.isfinite(markers)) or np.any((markers < z_min) | (markers > z_max)):
        raise ValueError("breakpoints must be a finite 1D sequence within the sampling bounds")
    intervals = (z_max - z_min) / max_dz
    if not np.isfinite(intervals) or intervals > max_nodes - 1:
        raise ValueError("max_nodes budget is too small for the requested max_dz")
    grid = np.linspace(z_min, z_max, max(1, int(np.ceil(intervals))) + 1)
    nodes = np.concatenate((grid, markers, [z_min, z_max]))
    # Explicit markers take precedence over generated grid points; the interval
    # bounds in turn take precedence over roundoff-close interior markers.
    priorities = np.concatenate((np.zeros(len(grid)), np.ones(len(markers)), [2, 2]))
    z = _merge_sampling_nodes(nodes, priorities)
    if len(z) < 2:
        raise ValueError("sampling bounds are too close to floating-point resolution")
    if len(z) > max_nodes:
        raise ValueError("max_nodes budget is too small for the initial grid and breakpoints")
    y = _evaluate_function(function, z)
    for iteration in range(1, 101):
        z, y, controls = prepare_table(z, y)
        probes = z[:-1, None] + np.diff(z)[:, None] * np.array([0.25, 0.5, 0.75])
        actual = _evaluate_function(function, probes.ravel()).reshape(probes.shape)
        error = np.abs(BPoly(controls, z)(probes) - actual)
        bad = np.any(error > atol + rtol * np.abs(actual), axis=1)
        if not np.any(bad):
            return (
                z,
                y,
                controls,
                {
                    "n_nodes": len(z),
                    "max_probe_error": float(np.max(error)),
                    "atol": atol,
                    "rtol": rtol,
                    "max_dz": max_dz,
                    "iterations": iteration,
                },
            )
        additions = probes[bad].ravel()
        if len(z) + len(additions) > max_nodes:
            raise ValueError("max_nodes budget exhausted before the function met the requested probe tolerance")
        new_z = np.unique(np.concatenate((z, additions)))
        if len(new_z) == len(z) or np.any(np.diff(new_z) <= _node_resolution(new_z[:-1], new_z[1:])):
            raise ValueError("floating-point redshift resolution exhausted before the requested tolerance was reached")
        # Retain existing evaluations; callables are expected to be deterministic.
        new_y = np.empty(len(new_z))
        old_indices = np.searchsorted(new_z, z)
        new_y[old_indices] = y
        new_mask = np.ones(len(new_z), dtype=bool)
        new_mask[old_indices] = False
        new_y[new_mask] = _evaluate_function(function, new_z[new_mask])
        z, y = new_z, new_y
    raise ValueError("function sampling exceeded 100 refinement iterations; relax tolerance or specify feature scales")
