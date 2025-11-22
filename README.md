# Monte Carlo Forward Model for Photon Transport

Python implementation of Jacques’ classic **`mc321.c`**  (Chapter 5 of *Steven L. Jacques, 2011*).

This project simulates **photon transport from an isotropic point source** in an infinite, homogeneous turbid medium using the standard  
**HOP → DROP → SPIN → ROULETTE** Monte Carlo workflow.

Absorbed energy is recorded in **spherical**, **cylindrical**, and **planar** geometries and converted into fluence using Jacques’ normalization formulas (Appendix 5.7).

## 📁 Folder Structure
```
mc-forward-jacques/
│
├── mc_forward_jacques.py     # Core Monte Carlo engine
├── demo_forward.py           # Runs simulation and generates Fig. 1.7
│
├── figs/
│ └── Fsph_Fcyl_Fpla.png      # Jacques-style fluence plot
│
├── report/
│ └── Research_Project_Monte_Carlo_Forward_Model.pdf
│
└── README.md
```

## 🚀 Key Features

- Faithful Python translation of **Jacques' Monte Carlo model**
- Exact implementation of:
  - **HOP** – exponential free-path sampling  
  - **DROP** – absorption + scoring  
  - **SPIN** – Henyey–Greenstein scattering  
  - **ROULETTE** – termination of low-weight photons  
- Correct HG sampling formula (fixes the known book typo)
- Fluence output:
  - `T_sph(r)` — spherical shells  
  - `T_cyl(r)` — cylindrical shells  
  - `T_pla(r)` — planar slabs  
- Output fluence curves reproduce **Fig. 5.10** from Jacques (2011)
- Includes full **LaTeX report** with derivations and figures


## ▶️ How to Run

Run the forward Monte Carlo simulation and produce the Jacques-style fluence plot:

```bash
python demo_forward.py
```
This generates:
```bash
figs/Fsph_Fcyl_Fpla.png
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
Jacques, S. L. (2011). Optical Properties of Biological Tissues: A Review.<br>
Chapter 5 provides the original mc321.c Monte Carlo model used here.
