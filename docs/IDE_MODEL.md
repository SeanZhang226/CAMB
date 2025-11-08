# IDE (Interacting Dark Energy) Model in CAMB

## Overview

This implementation provides complete support for Interacting Dark Energy (IDE) models in CAMB, ported from [liaocrane/IDECAMB](https://github.com/liaocrane/IDECAMB).

## Features

### Dark Energy Equation of State Forms

1. **CPL (Chevallier-Polarski-Linder)**: w(a) = w₀ + wₐ(1-a)
2. **HDE (Holographic Dark Energy)**: w(a) = -1/3 - 2√(ρ_de/3)/(3cH)
3. **NADE (New Agegraphic Dark Energy)**: w(a) = -1 + 2√(ρ_de/3)/(3anH)

### Coupling Forms

The energy transfer between dark energy and dark matter can take several forms:

1. **Hrde**: Q = β H ρ_de (coupling proportional to Hubble × dark energy density)
2. **Hrc**: Q = β H ρ_c (coupling proportional to Hubble × dark matter density)
3. **H0rde**: Q = β H₀ ρ_de (coupling proportional to H₀ × dark energy density)
4. **H0rc**: Q = β H₀ ρ_c (coupling proportional to H₀ × dark matter density)

### Covariant Forms

- **uc**: Coupling in dark matter rest frame
- **ude**: Coupling in dark energy rest frame

## Usage

### Python Interface

```python
import camb
from camb.dark_energy import DarkEnergyIDE

# Create CAMB parameters
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)

# Set up IDE model
ide = DarkEnergyIDE()
ide.set_params(
    w=-0.95,           # w(z=0)
    wa=0.0,            # -dw/da(0)
    beta=0.02,         # coupling strength
    coupling_form=1,   # 1=Hrde, 2=Hrc, 3=H0rde, 4=H0rc
    eos_form=1,        # 1=CPL, 2=HDE, 3=NADE
    covariant_form=1,  # 1=uc, 2=ude
    cs2=1.0            # sound speed squared
)

# Assign to CAMB parameters
pars.DarkEnergy = ide

# Calculate results
pars.set_for_lmax(2500, lens_potential_accuracy=0)
pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)

results = camb.get_results(pars)
```

### INI File Configuration

```ini
dark_energy_model = ide
w = -0.95
wa = 0.0
beta = 0.02
coupling_form = 1
eos_form = 1
covariant_form = 1
cs2_lam = 1.0
```

## Implementation Details

### Modules

1. **IDEtools.f90**: Utility functions
   - `ConformalH`: Computes H(a) including all energy components
   - `Hermite`: Hermite cubic interpolation for background queries

2. **CoupledFluidModels.f90**: Background evolution solver
   - `Solve_CF`: Main ODE solver using dverk integrator
   - `Eqs_CF`: Background evolution equations
   - `EoS_CF`: Equation of state functions
   - `Coup_CF`: Coupling functions
   - `PerturCoupC_CF`, `PerturCoupD_CF`: Perturbation couplings
   - `Get_grhodea2_grhoca2`: Interpolated background queries

3. **DarkEnergyIDE.f90**: Main IDE model
   - Inherits from `TDarkEnergyEqnOfState`
   - Implements full CAMB dark energy interface
   - Handles initialization and perturbation evolution

### Numerical Methods

- **ODE Integration**: Runge-Kutta-Verner (dverk) method with tolerance 10⁻⁸
- **Interpolation**: Hermite cubic interpolation for smooth background queries
- **Time Stepping**: 2000 logarithmic steps from a = 10⁻¹² to a = 1
- **Initial Conditions**: Broyden iteration for NADE model

### Background Evolution Equations

The coupled evolution equations are:

```
dρ_de/dτ + 3H(1+w)ρ_de = Q
dρ_c/dτ + 3Hρ_c = -Q
```

where τ is conformal time, and Q is the energy transfer rate.

## Physical Interpretation

- **β > 0**: Energy flows from dark energy to dark matter
- **β < 0**: Energy flows from dark matter to dark energy
- **β = 0**: No coupling (standard dark energy)

The coupling affects:
- Background expansion history H(z)
- Dark energy equation of state w(z)
- Growth of structure
- CMB anisotropies
- Matter power spectrum

## Validation

The implementation has been validated to reproduce results from liaocrane/IDECAMB for:
- Background evolution ρ_de(a), ρ_c(a)
- Equation of state w(a)
- Coupling term Q(a)

## References

- Liaocrane IDECAMB: https://github.com/liaocrane/IDECAMB
- CPL parametrization: Chevallier & Polarski (2001), Linder (2003)
- IDE models: Wang et al. (2016), arXiv:1606.07815

## Testing

A Fortran test program is provided in `fortran/test_ide.f90`:

```bash
cd fortran/Release
gfortran -O3 -fopenmp -c ../test_ide.f90 -o test_ide.o
gfortran -O3 -fopenmp test_ide.o libcamb.a -lforutils -o test_ide
./test_ide
```

## Notes

- The IDE solver is initialized during `Init()` call
- Background evolution is precomputed and stored in interpolation tables
- Perturbation equations include coupling terms in both density and velocity
- Dark matter velocity evolution (v_c) is enabled when needed based on coupling form
