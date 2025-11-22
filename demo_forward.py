"""
demo_forward.py

Author: Binger Yu
Date: 2025-11-22
version: 0.1

Generate Jacques-style fluence plots (spherical, cylindrical, planar)
using the forward Monte Carlo simulator implemented in `mc_forward_jacques.py`.

This script:
1. Runs the MC simulation with Jacques (2011) optical parameters.
2. Extracts fluence T(r) for r ≤ 1 cm.
3. Produces the Jacques-style log-scale fluence plot (Fig. 5.10).
4. Saves the figure to figs/Fsph_Fcyl_Fpla.png.
"""

import os
os.environ.pop("MPLBACKEND", None)   # avoid backend override issues

import matplotlib
matplotlib.use("MacOSX")             # ensures correct display on macOS

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext

from pathlib import Path

from mc_forward_jacques import run_mc


# -------------------------------------------------------------
# 1. Run Monte Carlo with Jacques-like optical parameters
# -------------------------------------------------------------
r, Fsph, Fcyl, Fpla = run_mc(
    Nphotons=10_000,
    mua=1.673,        # absorption [cm^-1]
    mus=312.0,        # scattering [cm^-1]
    g=0.90,           # anisotropy
    radial_size=2.0,  # scoring range [cm]
    NR=500,
)

# -------------------------------------------------------------
# 2. Remove overflow bin and restrict to r ≤ 1 cm
# -------------------------------------------------------------
r_plot    = r[:-1]
Fsph_plot = Fsph[:-1]
Fcyl_plot = Fcyl[:-1]
Fpla_plot = Fpla[:-1]

mask = (r_plot >= 0.0) & (r_plot <= 1.0)
r_plot    = r_plot[mask]
Fsph_plot = Fsph_plot[mask]
Fcyl_plot = Fcyl_plot[mask]
Fpla_plot = Fpla_plot[mask]

# Ensure output folder exists
out_dir = Path("figs")
out_dir.mkdir(exist_ok=True)


# -------------------------------------------------------------
# 3. Plotting — Jacques-style fluence curves (Fig. 5.10 mimic)
# -------------------------------------------------------------
fig, ax = plt.subplots()

# Three geometry curves
ax.semilogy(r_plot, Fsph_plot, '-',  label="spherical (Fsph)")
ax.semilogy(r_plot, Fcyl_plot, '--', label="cylindrical (Fcyl)")
ax.semilogy(r_plot, Fpla_plot, ':',  label="planar (Fpla)")

# Axis limits (approximately matching Jacques Fig. 5.10)
ax.set_xlim(0.0, 1.0)
ax.set_ylim(1e-5, 1e4)
ax.margins(x=0)

# Tick placement similar to the book
ax.set_xticks(np.linspace(0.0, 1.0, 6))   # 0, 0.2, ..., 1.0
ax.set_yticks([1e-4, 1e-2, 1.0, 1e2])

# y-axis ticks and labels: 10^3, 10^2, …, 10^-4
ax.yaxis.set_major_locator(LogLocator(base=10.0,
                                      subs=(1.0,)))  # only decades
ax.yaxis.set_major_formatter(LogFormatterMathtext())

# light horizontal rules (like the book)
ax.grid(True, which="both", axis="y", linestyle="-", linewidth=0.5)

# Labels and title
ax.set_xlabel(r"$r$ [cm]")
ax.set_ylabel(r"$T(r)$ [1/cm$^2$]")
ax.set_title("Fluence from Isotropic Point Source (Jacques-style)")
ax.legend()

# Layout and save
fig.tight_layout()
fig.savefig(out_dir / "Fsph_Fcyl_Fpla.png", dpi=300)
plt.close(fig)

print("Saved plot to", out_dir / "Fsph_Fcyl_Fpla.png")
