# Monte Carlo Forward Model (Finite Slab)
Date: November 25, 2025

Version: 1.1

## Project Overview

This project implements a Forward Monte Carlo simulation for photon transport in a finite slab geometry with Fresnel boundaries. It is adapted from the methods described by Steven Jacques [1] and the standard MCML software [2, 3].

This version (v1.1) implements the Finite Slab Geometry, representing a "Pencil Beam" experiment.

## 📁 Folder Structure
```
mc-forward-jacques/
├── mc_forward_jacques.py     # The core simulation logic (Finite Slab).
├── demo_forward.py           # A script to run the simulation and generate fluence plots.
├── figs/                     # Contains output plots.
├── report/                   # Project report (LaTeX source)
└── README.md
```

## 🚀 Key Features

- Geometry: Finite tissue slab of thickness d with Air/Tissue interfaces.
- Source: Collimated "Pencil Beam" incident perpendicular to the surface (Impulse Response).
- Physics:
  - Fresnel Reflection & Refraction: Implemented at top and bottom boundaries using Snell's Law.
  - Henyey-Greenstein Scattering: Models anisotropic scattering (g).
- Verification: Tracks Total Energy (Reflectance + Transmittance + Absorption) to ensure conservation = 1.0.


## ▶️ How to Run

Run the forward Monte Carlo simulation and produce the Jacques-style fluence plot:

```bash
python demo_forward.py
```
This generates:
```bash
figs/Fluence_Slab_Comparison.png
```

## 📦 Requirements
- Python 3.10+
- NumPy
- Matplotlib

Install dependencies:
```bash
pip install numpy matplotlib
```

## 📚 Reference
1. S. L. Jacques, Optical-Thermal Response of Laser-Irradiated Tissue, 2011. [PDF Chapter](https://omlc.org/software/mc/Jacques2011_MonteCarlo_Welch&VanGemert.pdf)
2. L. Wang and S. L. Jacques, "Monte Carlo modeling of light transport in multi-layered tissues in standard C," Computer Methods and Programs in Biomedicine, 1995. [PDF Paper](https://omlc.org/software/mc/mcpubs/1995LWCMPBMcml.pdf)
3. Oregon Medical Laser Center (OMLC), Monte Carlo Software Resources. [Website](https://omlc.org/software/mc/)
