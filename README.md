# Monte Carlo Forward Model for Photon Transport

Python implementation of Jacques' `mc321.c` (Chapter 5, Jacques 2011).  
Simulates photon transport from an isotropic point source in an infinite,
homogeneous medium and records absorbed energy in spherical, cylindrical,
and planar geometries.

## Files

- `mc_forward_jacques.py`  
  Core Monte Carlo code (LAUNCH–HOP–DROP–SPIN–ROULETTE) and fluence
  calculation.

- `demo_forward.py`  
  Runs the simulation with Jacques-like parameters and generates
  `figs/Fluence_Slab_Comparison.png` (Jacques-style fluence plot).

- `figs/Fluence_Slab_Comparison.png`  
  Final figure used in the report.

## How to run

```bash
python demo_forward.py
