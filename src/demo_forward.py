"""
demo_forward.py

Author: Binger Yu
Date: 2025-11-28
version: 1.1

Generate Jacques-style fluence plots for a Finite Slab using
the forward Monte Carlo simulator `mc_forward_jacques.py`.

Key Fixes in v1.1:
- Handles "Data has no positive values" error by masking zeros before log-plotting.
- Prints debug statistics (Max Fluence) to verify simulation data.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from mc_forward_jacques import run_mc

# -------------------------------------------------------------
# 1. Setup Parameters
# -------------------------------------------------------------
THICKNESS = 1.0     # Slab thickness [cm]
RADIAL_SIZE = 1.0   # Max scoring distance [cm]
NR_BINS = 100       # Number of bins

# -------------------------------------------------------------
# 2. Run Monte Carlo (Finite Slab)
# -------------------------------------------------------------
print(f"Running MC Simulation for a {THICKNESS} cm slab...")

r_centers, Fsph, Fcyl, Fpla, R_tot, T_tot = run_mc(
    Nphotons=10_000,
    mua=1.0,           # absorption [cm^-1]
    mus=100.0,         # scattering [cm^-1]
    g=0.90,            # anisotropy
    nt=1.33,           # Refractive index (Tissue)
    n_env=1.0,         # Refractive index (Air)
    thickness=THICKNESS,
    radial_size=RADIAL_SIZE,
    NR=NR_BINS,
)

# Verify Energy Conservation
Absorbed_tot = 1.0 - (R_tot + T_tot)
print("-" * 30)
print(f"MACROSCOPIC RESULTS:")
print(f"  Diffuse Reflectance (R): {R_tot:.4f}")
print(f"  Total Transmittance (T): {T_tot:.4f}")
print(f"  Total Absorbed (A):      {Absorbed_tot:.4f}")
print(f"  Conservation (R+T+A):    {R_tot + T_tot + Absorbed_tot:.4f}")
print("-" * 30)

# -------------------------------------------------------------
# 3. Prepare Plotting Data (Masking Zeros)
# -------------------------------------------------------------
# Remove overflow bin
r_plot = r_centers[:-1]
Fsph_plot = Fsph[:-1]
Fcyl_plot = Fcyl[:-1]
Fpla_plot = Fpla[:-1]

# Debug: Check for data presence
print(f"DATA DEBUG:")
print(f"  Max Fsph: {np.max(Fsph_plot):.4e}")
print(f"  Max Fcyl: {np.max(Fcyl_plot):.4e}")
print(f"  Max Fpla: {np.max(Fpla_plot):.4e}")

# Helper to mask zeros for log plots
def mask_zeros(arr):
    arr_masked = arr.copy()
    # Replace <= 0 with NaN so matplotlib ignores them in log plots
    arr_masked[arr_masked <= 0] = np.nan
    return arr_masked

Fsph_log = mask_zeros(Fsph_plot)
Fcyl_log = mask_zeros(Fcyl_plot)
Fpla_log = mask_zeros(Fpla_plot)

# Ensure output folder exists
out_dir = Path("../figs")
out_dir.mkdir(exist_ok=True)

# -------------------------------------------------------------
# 4. Plotting
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the three fluence profiles
# We use the masked arrays to avoid "UserWarning: Data has no positive values"
ax.semilogy(r_plot, Fsph_log, '-',  linewidth=1.5, label=r"$F_{sph}$ (vs distance $r$)")
ax.semilogy(r_plot, Fcyl_log, '--', linewidth=1.5, label=r"$F_{cyl}$ (vs radius $\rho$)")
ax.semilogy(r_plot, Fpla_log, ':',  linewidth=2.0, label=r"$F_{pla}$ (vs depth $z$)")

# Limits & Formatting
ax.set_xlim(0.0, RADIAL_SIZE)

# Dynamic Y-limits: Find min/max of valid positive data
all_valid_data = np.concatenate([
    Fsph_plot[Fsph_plot > 0],
    Fcyl_plot[Fcyl_plot > 0],
    Fpla_plot[Fpla_plot > 0]
])

if len(all_valid_data) > 0:
    y_min = np.min(all_valid_data) * 0.5
    y_max = np.max(all_valid_data) * 2.0
    ax.set_ylim(y_min, y_max)
else:
    print("\nWARNING: No positive fluence data found to plot!")

ax.grid(True, which="both", axis="y", linestyle="-", linewidth=0.5, alpha=0.5)
ax.grid(True, which="major", axis="x", linestyle="-", linewidth=0.5, alpha=0.5)

# Labels
ax.set_xlabel(r"Distance / Depth [cm]")
ax.set_ylabel(r"Fluence Rate [1/cm$^2$]")
ax.set_title(f"Fluence Profiles (Slab d={THICKNESS}cm)")
ax.legend()

# Layout and save
fig.tight_layout()
output_path = out_dir / "Fluence_Slab_Comparison.png"
fig.savefig(output_path, dpi=300)
plt.close(fig)

print(f"\nSaved plot to {output_path}")