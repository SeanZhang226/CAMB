#!/usr/bin/env python3
"""
Test script for the Interacting Dark Energy (IDE) model integration into CAMB.

This script validates the IDE model implementation by:
1. Testing background evolution with 10 randomized parameter sets
2. Testing full perturbation calculations with a specific parameter set

Usage:
    python test_ide_model.py
"""

import sys
import numpy as np

# Add current directory to path to import local camb
sys.path.insert(0, '.')
import camb

def test_background_evolution():
    """Test background evolution with 10 randomized parameter sets."""
    print("=" * 70)
    print("TEST 1: Background Evolution (10 Parameter Sets)")
    print("=" * 70)
    
    np.random.seed(42)  # For reproducibility
    
    successful_runs = 0
    failed_runs = 0
    
    for i in range(10):
        print(f"\n--- Run {i+1}/10 ---")
        
        try:
            # Randomize parameters within physically reasonable ranges
            H0 = np.random.uniform(60, 80)  # Hubble constant
            ombh2 = np.random.uniform(0.019, 0.026)  # Baryon density
            omch2 = np.random.uniform(0.08, 0.15)  # CDM density
            xi_de = np.random.uniform(-0.1, 0.1)  # DE coupling
            xi_c = np.random.uniform(-0.1, 0.1)  # DM coupling
            
            # Create parameters
            pars = camb.CAMBparams()
            
            # Set cosmological parameters
            pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
            
            # Set dark energy model to IDE
            pars.set_dark_energy(dark_energy_model='ide', w=-1.0, wa=0.0)
            
            # Set IDE coupling parameters
            pars.DarkEnergy.set_params(xi_de=xi_de, xi_c=xi_c)
            
            # Calculate background quantities
            results = camb.get_background(pars)
            
            # Get some derived quantities
            # For age, we can use cosmological calculator
            zre = pars.get_zre()
            
            print(f"Parameters: H0={H0:.2f}, ombh2={ombh2:.4f}, omch2={omch2:.4f}")
            print(f"            xi_de={xi_de:.4f}, xi_c={xi_c:.4f}")
            print(f"Derived: zre={zre:.3f}")
            print("Status: SUCCESS")
            
            successful_runs += 1
            
        except Exception as e:
            print(f"Parameters: H0={H0:.2f}, ombh2={ombh2:.4f}, omch2={omch2:.4f}")
            print(f"            xi_de={xi_de:.4f}, xi_c={xi_c:.4f}")
            print(f"Status: FAILED")
            print(f"Error: {str(e)}")
            failed_runs += 1
    
    print("\n" + "=" * 70)
    print(f"Background Evolution Test Summary:")
    print(f"  Successful: {successful_runs}/10")
    print(f"  Failed: {failed_runs}/10")
    print("=" * 70)
    
    return successful_runs, failed_runs


def test_full_perturbations():
    """Test full perturbation calculation with specific IDE parameters."""
    print("\n" + "=" * 70)
    print("TEST 2: Full Perturbation Calculation")
    print("=" * 70)
    
    try:
        # Set specific parameters
        H0 = 67.5
        ombh2 = 0.0224
        omch2 = 0.120
        xi_de = 0.05
        xi_c = 0.01
        
        print(f"\nParameters:")
        print(f"  H0 = {H0}")
        print(f"  ombh2 = {ombh2}")
        print(f"  omch2 = {omch2}")
        print(f"  xi_de = {xi_de}")
        print(f"  xi_c = {xi_c}")
        
        # Create parameters
        pars = camb.CAMBparams()
        
        # Set cosmological parameters
        pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
        
        # Set dark energy model to IDE
        pars.set_dark_energy(dark_energy_model='ide', w=-1.0, wa=0.0)
        
        # Set IDE coupling parameters
        pars.DarkEnergy.set_params(xi_de=xi_de, xi_c=xi_c)
        
        # Request CMB power spectra
        pars.set_for_lmax(2500, lens_potential_accuracy=0)
        
        # Request matter power spectrum
        pars.set_matter_power(redshifts=[0.], kmax=2.0)
        
        print("\nRunning CAMB with IDE model...")
        
        # Calculate results
        results = camb.get_results(pars)
        
        # Get CMB power spectra
        powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')
        totCL = powers['total']
        
        # Get matter power spectrum
        k, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints=200)
        
        print("\n✓ Full perturbation calculation successful!")
        print(f"\nFirst 5 TT power spectrum values (l*(l+1)*C_l/2π in μK²):")
        for l in range(2, 7):
            print(f"  l={l}: {totCL[l, 0]:.6e}")
        
        print(f"\nMatter power spectrum calculated at {len(k)} k-points")
        print(f"  k range: {k[0]:.6e} to {k[-1]:.6e} h/Mpc")
        print(f"  P(k) at k=0.1 h/Mpc: {pk[0, np.argmin(np.abs(k - 0.1))]:.6e} (Mpc/h)³")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Full perturbation calculation FAILED")
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main test routine."""
    print("\n" + "#" * 70)
    print("# CAMB IDE Model Integration Test Suite")
    print("# Testing Interacting Dark Energy (IDE) implementation")
    print("#" * 70)
    
    # Test 1: Background evolution with multiple parameter sets
    bg_success, bg_failed = test_background_evolution()
    
    # Test 2: Full perturbation calculation
    pert_success = test_full_perturbations()
    
    # Final summary
    print("\n" + "#" * 70)
    print("# FINAL TEST SUMMARY")
    print("#" * 70)
    print(f"Background Evolution Tests: {bg_success}/10 passed")
    print(f"Full Perturbation Test: {'PASSED' if pert_success else 'FAILED'}")
    
    if bg_success == 10 and pert_success:
        print("\n✓ ALL TESTS PASSED!")
        print("#" * 70)
        return 0
    else:
        print("\n✗ SOME TESTS FAILED")
        print("#" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
