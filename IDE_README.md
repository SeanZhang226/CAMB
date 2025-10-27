# Interacting Dark Energy (IDE) Model for CAMB

This implementation adds support for Interacting Dark Energy (IDE) models to CAMB, allowing for energy exchange between dark energy and dark matter.

## Overview

The IDE model introduces coupling between dark energy and dark matter through two coupling parameters:
- `xi_de`: Coupling constant proportional to dark energy density
- `xi_c`: Coupling constant proportional to dark matter density

This implementation is based on the IDECAMB model (https://github.com/liaocrane/IDECAMB) but has been modernized and integrated cleanly into the latest CAMB structure.

## Installation

The IDE model is already integrated into this fork of CAMB. Just build CAMB as normal:

```bash
cd fortran
make
```

For Python interface:
```bash
cd fortran
make python
```

## Usage

### Python Interface

```python
import camb

# Create parameters
pars = camb.CAMBparams()

# Set cosmological parameters
pars.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.120)

# Set IDE model with coupling parameters
pars.set_dark_energy(dark_energy_model='ide', w=-1.0, wa=0.0)
pars.DarkEnergy.set_params(xi_de=0.05, xi_c=0.01)

# Calculate results
pars.set_for_lmax(2500)
results = camb.get_results(pars)

# Get CMB power spectra
powers = results.get_cmb_power_spectra(pars)

# Get matter power spectrum
pars.set_matter_power(redshifts=[0.], kmax=2.0)
k, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints=200)
```

### INI File Configuration

For the Fortran executable, you can create an ini file with the following parameters:

```ini
# Dark energy model
dark_energy_model = ide

# Standard dark energy parameters
w = -1.0
cs2_lam = 1

# IDE coupling parameters
xi_de = 0.05   # Coupling proportional to DE density
xi_c = 0.01    # Coupling proportional to DM density
```

## Testing

Run the included test suite to verify the installation:

```bash
python test_ide_model.py
```

This runs:
1. 10 background evolution tests with randomized parameters
2. 1 full perturbation calculation test including CMB and matter power spectra

Expected output:
```
Background Evolution Tests: 10/10 passed
Full Perturbation Test: PASSED
✓ ALL TESTS PASSED!
```

## Model Parameters

### Coupling Parameters

- **xi_de** (default: 0.0): Coupling constant proportional to dark energy density
  - Controls energy transfer rate from/to dark energy
  - Physically reasonable range: |xi_de| < 0.1

- **xi_c** (default: 0.0): Coupling constant proportional to dark matter density
  - Controls energy transfer rate from/to dark matter
  - Physically reasonable range: |xi_c| < 0.1

### Standard Dark Energy Parameters

The IDE model also supports standard dark energy parameters from `TDarkEnergyEqnOfState`:
- **w**: Equation of state parameter (default: -1.0)
- **wa**: Evolution parameter (default: 0.0)
- **cs2**: Sound speed squared (default: 1.0)

## Physics

The IDE model modifies the background evolution through coupling terms in the conservation equations:

```
ρ'_de + 3H(1 + w)ρ_de = Q
ρ'_c + 3Hρ_c = -Q
```

where Q represents the energy transfer rate controlled by xi_de and xi_c.

For perturbations, the model includes:
- Modified density perturbation evolution for DE and DM
- Velocity perturbations for both components
- Coupling terms in the perturbation equations

## Implementation Details

### Files Modified/Created

1. **fortran/DarkEnergyIDE.f90**: New module implementing IDE physics
2. **fortran/camb.f90**: Added IDE model option
3. **fortran/Makefile_main**: Added DarkEnergyIDE to build
4. **camb/dark_energy.py**: Added Python wrapper class
5. **test_ide_model.py**: Test suite

### Architecture

The implementation follows CAMB's modular dark energy framework:
- Extends `TDarkEnergyEqnOfState` for compatibility
- Implements all required methods for background and perturbation evolution
- Properly handles Python-Fortran interface through f2py

## Limitations

This is a simplified IDE implementation focusing on:
- Background evolution with small coupling corrections
- Linear perturbation theory
- Constant equation of state (w, wa parameterization)

For more sophisticated IDE models with:
- Non-linear coupling effects
- Time-varying coupling functions
- Multiple dark sector components

Consider using the full IDECAMB package or extending this implementation.

## References

If you use this IDE model, please cite:

1. Original IDECAMB papers:
   - [arXiv:1404.5220](https://arxiv.org/abs/1404.5220)
   - [arXiv:2306.01593](https://arxiv.org/abs/2306.01593)

2. CAMB:
   - Lewis, A., Challinor, A., & Lasenby, A. (2000). ApJ, 538, 473
   - Howlett et al. (2012). JCAP, 04, 027

## Contact

For issues or questions about this integration, please open an issue on the GitHub repository.

## License

This code follows CAMB's licensing terms. See the main CAMB repository for details.
