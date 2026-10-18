from ctypes import POINTER, byref, c_bool, c_double, c_int, c_void_p

import numpy as np

from .baseconfig import F2003Class, f_pointer, fortran_class, numpy_1d


class ReionizationModel(F2003Class):
    """
    Abstract base class for reionization models.

    The optional ``include_heating`` switch adds an approximate model of the photo-heating of the
    intergalactic medium during reionization. It is a very crude instantaneous-Jeans approximation for
    an order of magnitude estimate. When enabled, the baryon temperature is smoothly raised
    towards ``reion_temperature`` (default :math:`10^4` K) following the ionization-fraction shape, and the
    baryon sound speed is mapped towards the corresponding ideal-gas value. This raises the baryon pressure,
    suppressing small-scale baryon (and hence total matter) power. It is off by default and only affects the
    low-redshift matter power spectrum.

    """

    _fields_ = (
        ("Reionization", c_bool, "Is there reionization? (can be off for matter power, which is independent of it)"),
        (
            "include_heating",
            c_bool,
            "Whether to include approximate smooth reionization heating following the x_e shape (default False)",
        ),
        (
            "reion_temperature",
            c_double,
            "Asymptotic gas temperature in K reached during reionization heating (default 1e4 K)",
        ),
    )

    def write_ini(self, state) -> None:
        state.set("reionization", self.Reionization)
        state.set("reion_include_heating", self.include_heating)
        state.set("reion_temperature", self.reion_temperature)


class BaseTauWithHeReionization(ReionizationModel):
    """
    Abstract class for models that map z_re to tau, and include second reionization of Helium
    """

    _fields_ = (
        ("use_optical_depth", c_bool, "Whether to use the optical depth or redshift parameters"),
        ("redshift", c_double, "Reionization redshift (xe=0.5) if use_optical_depth=False"),
        ("optical_depth", c_double, "Optical depth if use_optical_depth=True"),
        (
            "fraction",
            c_double,
            (
                "Reionization fraction when complete, "
                "or -1 for full ionization of hydrogen and first ionization of helium."
            ),
        ),
        ("include_helium_fullreion", c_bool, "Whether to include second reionization of helium"),
        ("helium_redshift", c_double, "Redshift for second reionization of helium"),
        ("helium_delta_redshift", c_double, "Width in redshift for second reionization of helium"),
        ("helium_redshiftstart", c_double, "Include second helium reionization below this redshift"),
        ("tau_solve_accuracy_boost", c_double, "Accuracy boosting parameter for solving for z_re from tau"),
        (
            "timestep_boost",
            c_double,
            "Accuracy boosting parameter for the minimum number of time sampling steps through reionization",
        ),
        ("max_redshift", c_double, "Maximum redshift allowed when mapping tau into reionization redshift"),
        ("__min_redshift", c_double, "Minimum redshift allowed when mapping tau into reionization redshift"),
        ("__fHe", c_double, "Helium fraction"),
        ("__state", f_pointer),
    )

    _fortran_class_module_ = "Reionization"
    _fortran_class_name_ = "TBaseTauWithHeReionization"

    _methods_ = (("GetZreFromTau", [c_void_p, POINTER(c_double)], c_double, {"nopass": True}),)

    def get_zre(self, params, tau=None):
        """
        Get the midpoint redshift of reionization.

        :param params: :class:`.model.CAMBparams` instance with cosmological parameters
        :param tau: if set, calculate the redshift for optical depth tau, otherwise uses currently set parameters
        :return: reionization mid-point redshift
        """
        if self.use_optical_depth or tau is not None:
            from .camb import CAMBparams

            assert isinstance(params, CAMBparams)
            return self.f_GetZreFromTau(byref(params), c_double(self.optical_depth if tau is None else tau))
        else:
            return self.redshift

    def set_zrei(self, zrei):
        """
        Set the mid-point reionization redshift

        :param zrei: mid-point redshift
        :return: self
        """
        self.use_optical_depth = False
        self.redshift = zrei
        return self

    def set_tau(self, tau):
        """
        Set the optical depth

        :param tau: optical depth
        :return: self
        """
        self.use_optical_depth = True
        self.optical_depth = tau
        return self

    def set_extra_params(self, max_zrei=None):
        """
        Set extra parameters (not tau, or zrei)

        :param max_zrei: maximum redshift allowed when mapping tau into reionization redshift
        """
        if max_zrei is not None:
            self.max_redshift = max_zrei

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.set("re_use_optical_depth", self.use_optical_depth)
        if self.use_optical_depth:
            state.set("re_optical_depth", self.optical_depth)
        else:
            state.set("re_redshift", self.redshift)
        state.set("re_ionization_frac", self.fraction)
        state.set("include_helium_fullreion", self.include_helium_fullreion)
        state.set("re_helium_redshift", self.helium_redshift)
        state.set("re_helium_delta_redshift", self.helium_delta_redshift)
        state.set("re_helium_redshiftstart", self.helium_redshiftstart)
        state.set("max_zrei", self.max_redshift)


@fortran_class
class TanhReionization(BaseTauWithHeReionization):
    """
    This default (unphysical) tanh x_e parameterization is described in
    Appendix B of `arXiv:0804.3865 <https://arxiv.org/abs/0804.3865>`_
    """

    _fields_ = (("delta_redshift", c_double, "Duration of reionization"),)

    _fortran_class_name_ = "TTanhReionization"

    def set_extra_params(self, deltazrei=None, max_zrei=None) -> None:
        """
        Set extra parameters (not tau, or zrei)

        :param deltazrei: delta z for reionization
        :param max_zrei: maximum redshift allowed when mapping tau into reionization
        """
        super().set_extra_params(max_zrei)
        if deltazrei is not None:
            self.delta_redshift = deltazrei

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.set("re_delta_redshift", self.delta_redshift)


@fortran_class
class ExpReionization(BaseTauWithHeReionization):
    """
    An ionization fraction that decreases exponentially at high z, saturating to fully ionized at fixed redshift.
    This model has a minimum non-zero tau around 0.04 for reion_redshift_complete=6.1.
    Similar to e.g. arXiv:1509.02785, arXiv:2006.16828, but not attempting to fit shape near x_e~1 at z<6.1
    """

    _fields_ = (
        ("reion_redshift_complete", c_double, "end of reionization"),
        ("reion_exp_smooth_width", c_double, "redshift scale to smooth exponential"),
        ("reion_exp_power", c_double, "power in exponential, exp(-lambda(z-redshift_complete)^exp_power)"),
    )

    _fortran_class_name_ = "TExpReionization"

    def set_extra_params(
        self, reion_redshift_complete=None, reion_exp_power=None, reion_exp_smooth_width=None, max_zrei=None
    ) -> None:
        """
        Set extra parameters (not tau, or zrei)

        :param reion_redshift_complete: redshift at which reionization complete (e.g. around 6)
        :param reion_exp_power: power in exponential decay with redshift
        :param reion_exp_smooth_width: smoothing parameter to keep derivative smooth
        :param max_zrei: maximum redshift allowed when mapping tau into reionization
        """
        super().set_extra_params(max_zrei)
        if reion_redshift_complete is not None:
            self.reion_redshift_complete = reion_redshift_complete
        if reion_exp_power is not None:
            self.reion_exp_power = reion_exp_power
        if reion_exp_smooth_width is not None:
            self.reion_exp_smooth_width = reion_exp_smooth_width

    def write_ini(self, state) -> None:
        super().write_ini(state)
        state.write_fields(self, names=("reion_redshift_complete", "reion_exp_power", "reion_exp_smooth_width"))


@fortran_class
class TabulatedReionization(BaseTauWithHeReionization):
    """External, homogeneous post-recombination ionization histories.

    ``ionization`` tables specify the H/first-He reionization contribution in
    electrons per hydrogen nucleus, excluding the neutral-region residual.
    ``ionized_fraction`` tables specify its normalized effective fraction.
    ``total_xe`` tables specify the total electron abundance, including residuals.
    These are distinct physical inputs, not interchangeable parameter names.

    The interpolation is range-preserving and C2. No global monotonicity or
    redshift-50 cutoff is imposed. This does not solve an ionization/heating
    model, modify recombination, or represent spatially patchy reionization.
    """

    _fields_ = (
        ("z_start", c_double, "High-redshift end including the explicit connection region"),
        ("z_complete", c_double, "Lowest tabulated redshift; not necessarily physical completion"),
        ("transition_width", c_double, "Explicit high-redshift connection width"),
        ("boundary_tolerance", c_double, "Maximum negligible endpoint contribution for automatic connection"),
        ("_history_mode", c_int),
        ("n_nodes", c_int, "Number of input history nodes"),
        ("min_steps", c_int, "Minimum source-time samples; not a thermodynamics accuracy guarantee"),
        ("source_samples_per_interval", c_int, "Maximum additional source samples in each varying history interval"),
        ("source_min_dlnz", c_double, "Smallest ln(1+z) spacing of those samples, before accuracy boosts"),
    )
    _fortran_class_name_ = "TTabulatedReionization"
    _methods_ = (
        ("SetHistory", [POINTER(c_int), numpy_1d, numpy_1d]),
        ("GetHistory", [POINTER(c_int), numpy_1d, numpy_1d]),
    )

    def set_table(
        self,
        redshifts,
        x_e,
        *,
        history_mode="ionization",
        transition_width=None,
        include_helium_fullreion=None,
        boundary_tolerance=1e-6,
    ):
        """Set a table and its physical interpretation; return this object.

        A zero high-z reionization contribution connects to recombination
        directly. A negligible nonzero endpoint receives a connection one last
        node-spacing wide. A significant unfinished endpoint requires an
        explicit positive ``transition_width``; it is never silently cut off.
        ``total_xe`` always requires an explicit positive connection width.
        The table is unchanged inside its specified range. Below it, the first
        value is held constant; partial late ionization is not forced to unity.

        By default the extra He II -> He III contribution is included for the
        two reionization-contribution modes, but excluded for ``total_xe``.
        Physical H/He electron bounds are checked against the actual CAMB YHe
        during initialization, not against a hard-coded electron ceiling.
        """
        from ._tabulated import prepare_table

        modes = {"ionization": 0, "ionized_fraction": 1, "total_xe": 2}
        if history_mode not in modes:
            raise ValueError("history_mode must be ionization, ionized_fraction or total_xe")
        z, values, controls = prepare_table(redshifts, x_e)
        boundary_tolerance = float(boundary_tolerance)
        if not np.isfinite(boundary_tolerance) or boundary_tolerance < 0:
            raise ValueError("boundary_tolerance must be finite and non-negative")
        if history_mode == "ionized_fraction" and np.any(values > 1):
            raise ValueError("ionized_fraction must lie between zero and one")
        if transition_width is None:
            if history_mode == "total_xe" or values[-1] > boundary_tolerance:
                raise ValueError(
                    "The high-redshift history has not been connected to recombination. "
                    "Extend the table to a neutral endpoint, or explicitly supply a positive "
                    "transition_width. total_xe requires an explicit connection width."
                )
            transition_width = z[-1] - z[-2] if values[-1] > 0 else 0.0
        transition_width = float(transition_width)
        if not np.isfinite(transition_width) or transition_width < 0:
            raise ValueError("transition_width must be finite and non-negative")
        if transition_width == 0 and (history_mode == "total_xe" or values[-1] != 0):
            raise ValueError("A nonzero endpoint or total_xe needs a positive transition_width")
        if not np.isfinite(z[-1] + transition_width):
            raise ValueError("The connection endpoint must be finite")
        if include_helium_fullreion is not None and not isinstance(include_helium_fullreion, (bool, np.bool_)):
            raise TypeError("include_helium_fullreion must be a boolean")
        self.transition_width = transition_width
        self.boundary_tolerance = boundary_tolerance
        self._history_mode = modes[history_mode]
        self.use_optical_depth = False
        self.optical_depth = 0
        self.include_helium_fullreion = (
            history_mode != "total_xe" if include_helium_fullreion is None else bool(include_helium_fullreion)
        )
        self.f_SetHistory(byref(c_int(len(z))), z, np.ascontiguousarray(controls.ravel(order="F")))
        return self

    def set_function(
        self,
        function,
        z_min,
        z_max,
        *,
        history_mode="ionization",
        transition_width=None,
        include_helium_fullreion=None,
        boundary_tolerance=1e-6,
        **sampling_options,
    ):
        """Sample a function into the same Fortran history; return diagnostics.

        ``atol``, ``rtol``, ``max_dz``, ``breakpoints`` and ``max_nodes`` control
        sampling. Probe errors refer to function values, not derivatives or C_l.
        Finite sampling cannot guarantee finding arbitrary unannounced narrow
        features. Supply breakpoints and resolve their scales with max_dz.
        No Python callbacks are used in the Fortran calculation.
        """
        from ._tabulated import sample_function

        z, values, _, diagnostics = sample_function(function, z_min, z_max, **sampling_options)
        self.set_table(
            z,
            values,
            history_mode=history_mode,
            transition_width=transition_width,
            include_helium_fullreion=include_helium_fullreion,
            boundary_tolerance=boundary_tolerance,
        )
        diagnostics = dict(diagnostics)
        diagnostics["transition_width"] = self.transition_width
        diagnostics["history_mode"] = history_mode
        return diagnostics

    @property
    def history_mode(self):
        return ("ionization", "ionized_fraction", "total_xe")[self._history_mode]

    def get_table(self):
        """Return independent copies of the sampled redshifts and input values.

        These are in the input mode's units, before residual/helium/connection
        processing. Save the mode and other options alongside these arrays.
        """
        if self.n_nodes < 2:
            raise ValueError("No ionization history has been set")
        z = np.empty(self.n_nodes)
        values = np.empty(self.n_nodes)
        self.f_GetHistory(byref(c_int(self.n_nodes)), z, values)
        return {"redshifts": z, "x_e": values}

    def set_tau(self, tau):
        raise ValueError("tau is derived from a tabulated history; tau-to-redshift inversion is not defined")

    def set_zrei(self, zrei):
        raise ValueError("A general history has no unique zrei; use set_table or set_function")

    def get_zre(self, params, tau=None):
        raise ValueError("A general history has no unique midpoint; inspect the history explicitly")

    def set_extra_params(self, max_zrei=None):
        if max_zrei is not None:
            raise ValueError("max_zrei controls tau inversion, not the range of a tabulated history")

    def write_ini(self, state):
        raise TypeError("TabulatedReionization is Python-only; .ini cannot serialize its history arrays")

    def __getstate__(self):
        raise TypeError("Save history arrays and options explicitly; this Fortran history cannot be pickled")
