"""
vis_trajectories.py

Author: Binger Yu
Date: 2025-11-29
version: 1.1

Visualizes random photon trajectories in a finite slab.

Generates TWO 3D plots:
1. Left-to-Right Flow (Side View): Ideal for schematic comparison.
   - Save as: figs/trajectories_side_view.png
2. Top-to-Bottom Flow (Top-Down View): Standard depth visualization.
   - Save as: figs/trajectories_top_down.png
"""

import math
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Global plotting style
# ----------------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "figure.autolayout": True,  # Helps prevent label cutoff
})

# ----------------------------------------------------------------------
# Physical / simulation parameters
# ----------------------------------------------------------------------
N_PHOTONS = 100
THICKNESS = 0.1  # [cm]
MUA = 10.0       # [1/cm]
MUS = 90.0       # [1/cm]
G = 0.9
NT = 1.33
N_ENV = 1.0

STEP_EPS = 1.0e-9  # safety for log(0)

# List of (x_arr, y_arr, z_arr, status_string)
paths: list[tuple[list[float], list[float], list[float], str]] = []


def run_trajectory_sim(seed: int = 42) -> None:
    """
    Run a simple Monte Carlo propagation for N_PHOTONS photons.
    """
    rng = np.random.default_rng(seed)
    mut = MUA + MUS

    # Specular reflection at z = 0 (normal incidence)
    r_sp = (N_ENV - NT) / (N_ENV + NT)
    R_sp = r_sp ** 2

    for _ in range(N_PHOTONS):
        # Initial position and direction (collimated beam at z = 0)
        x, y, z = 0.0, 0.0, 0.0
        ux, uy, uz = 0.0, 0.0, 1.0
        W = 1.0 - R_sp  # subtract specular component

        history_x = [x]
        history_y = [y]
        history_z = [z]

        status = "absorbed"

        while W > 1.0e-4:
            # 1. HOP
            rnd = STEP_EPS + (1.0 - STEP_EPS) * rng.random()
            s = -math.log(rnd) / mut

            # Distance to nearest slab boundary
            if uz > 0.0:
                d_bound = (THICKNESS - z) / uz
            elif uz < 0.0:
                d_bound = -z / uz
            else:
                d_bound = float("inf")

            if s > d_bound:
                # Move to boundary
                x += d_bound * ux
                y += d_bound * uy
                z += d_bound * uz
                s -= d_bound

                history_x.append(x)
                history_y.append(y)
                history_z.append(z)

                # Fresnel reflection / transmission
                ca1 = abs(uz)
                sa1_sq = 1.0 - ca1 ** 2
                val = (NT / N_ENV) ** 2 * sa1_sq

                if val >= 1.0:
                    R = 1.0  # total internal reflection
                else:
                    ca2 = math.sqrt(1.0 - val)
                    rs = (NT * ca1 - N_ENV * ca2) / (NT * ca1 + N_ENV * ca2)
                    rp = (NT * ca2 - N_ENV * ca1) / (NT * ca2 + N_ENV * ca1)
                    R = 0.5 * (rs ** 2 + rp ** 2)

                if rng.random() > R:
                    if uz < 0.0:
                        status = "reflected"  # leaves through top surface
                    else:
                        status = "transmitted"  # leaves through bottom
                    break
                else:
                    uz = -uz

            else:
                # Full step
                x += s * ux
                y += s * uy
                z += s * uz

                history_x.append(x)
                history_y.append(y)
                history_z.append(z)

                # 2. DROP
                dW = W * (MUA / mut)
                W -= dW

                # 3. SPIN
                rnd = rng.random()
                if G == 0.0:
                    costheta = 2.0 * rnd - 1.0
                else:
                    temp = (1.0 - G ** 2) / (1.0 - G + 2.0 * G * rnd)
                    costheta = (1.0 + G ** 2 - temp ** 2) / (2.0 * G)

                sintheta = math.sqrt(max(0.0, 1.0 - costheta ** 2))
                psi = 2.0 * math.pi * rng.random()
                cpsi, spsi = math.cos(psi), math.sin(psi)

                if 1.0 - abs(uz) <= 1.0e-12:
                    uxx = sintheta * cpsi
                    uyy = sintheta * spsi
                    uzz = costheta * math.copysign(1.0, uz)
                else:
                    temp = math.sqrt(1.0 - uz * uz)
                    uxx = (
                            sintheta * (ux * uz * cpsi - uy * spsi) / temp
                            + ux * costheta
                    )
                    uyy = (
                            sintheta * (uy * uz * cpsi + ux * spsi) / temp
                            + uy * costheta
                    )
                    uzz = -sintheta * cpsi * temp + uz * costheta

                ux, uy, uz = uxx, uyy, uzz

                # 4. ROULETTE
                if W < 1.0e-3:
                    if rng.random() <= 0.1:
                        W /= 0.1
                    else:
                        break

        paths.append((history_x, history_y, history_z, status))


def _set_equal_aspect_3d(ax, xs, ys, zs) -> None:
    """Set equal aspect ratio for a 3D axes based on data ranges."""
    x_min, x_max = xs.min(), xs.max()
    y_min, y_max = ys.min(), ys.max()
    z_min, z_max = zs.min(), zs.max()

    max_range = max(x_max - x_min, y_max - y_min, z_max - z_min)

    x_mid = 0.5 * (x_max + x_min)
    y_mid = 0.5 * (y_max + y_min)
    z_mid = 0.5 * (z_max + z_min)

    half = 0.5 * max_range
    ax.set_xlim(x_mid - half, x_mid + half)
    ax.set_ylim(y_mid - half, y_mid + half)
    ax.set_zlim(z_mid - half, z_mid + half)


def _common_plot_elements(ax, title):
    """Helper to set common plot styles and legend."""
    color_map = {"reflected": "tab:red", "transmitted": "tab:green", "absorbed": "tab:blue"}
    ax.plot([], [], [], color=color_map["reflected"], label="Reflected")
    ax.plot([], [], [], color=color_map["transmitted"], label="Transmitted")
    ax.plot([], [], [], color=color_map["absorbed"], label="Absorbed")
    ax.legend(loc="upper right", frameon=False)
    ax.set_title(title)
    # Label padding to prevent overlap - Increased to 12
    ax.xaxis.labelpad = 12
    ax.yaxis.labelpad = 12
    ax.zaxis.labelpad = 12


def plot_left_to_right(output_path: str = "figs/trajectories_side_view.png") -> None:
    """
    Plot 1: Side View (Left-to-Right).
    Simulation Depth (z) -> Plot X-axis.
    """
    all_px, all_py, all_pz = [], [], []
    plot_paths = []

    for sim_x, sim_y, sim_z, status in paths:
        px = np.array(sim_z)  # Depth -> X
        py = np.array(sim_y)  # Lateral -> Y
        pz = np.array(sim_x)  # Lateral -> Z
        plot_paths.append((px, py, pz, status))
        all_px.append(px);
        all_py.append(py);
        all_pz.append(pz)

    all_px = np.concatenate(all_px)
    all_py = np.concatenate(all_py)
    all_pz = np.concatenate(all_pz)

    fig = plt.figure(figsize=(7, 6))  # Increased height to 6
    ax = fig.add_subplot(111, projection="3d")

    # Boundaries (Vertical Planes at X=0 and X=d)
    y_span = np.linspace(all_py.min(), all_py.max(), 10)
    z_span = np.linspace(all_pz.min(), all_pz.max(), 10)
    YY, ZZ = np.meshgrid(y_span, z_span)

    ax.plot_surface(np.zeros_like(YY), YY, ZZ, alpha=0.15, color="grey", linewidth=0)
    ax.plot_surface(np.full_like(YY, THICKNESS), YY, ZZ, alpha=0.15, color="grey", linewidth=0)

    color_map = {"reflected": "tab:red", "transmitted": "tab:green", "absorbed": "tab:blue"}
    for px, py, pz, status in plot_paths:
        ax.plot(px, py, pz, color=color_map.get(status), linewidth=0.6, alpha=0.6)

    ax.set_xlabel(r"Depth $z$ [cm]")
    ax.set_ylabel(r"Lateral $y$ [cm]")
    ax.set_zlabel(r"Lateral $x$ [cm]")

    _set_equal_aspect_3d(ax, all_px, all_py, all_pz)
    ax.view_init(elev=10, azim=-80)

    # Increase padding for side view labels specifically
    ax.xaxis.labelpad = 15
    ax.yaxis.labelpad = 15
    ax.zaxis.labelpad = 15

    _common_plot_elements(ax, rf"Photon Trajectories in a Finite Slab (Side View, $N={N_PHOTONS}$)")

    # Use generous padding to ensure labels aren't cut off
    fig.tight_layout(pad=3.0)
    fig.savefig(output_path, bbox_inches='tight')
    print(f"Saved {output_path}")


def plot_top_to_bottom(output_path: str = "figs/trajectories_top_down.png") -> None:
    """
    Plot 2: Top-Down View.
    Simulation Depth (z) -> Plot Z-axis (Inverted).
    """
    all_px, all_py, all_pz = [], [], []
    plot_paths = []

    for sim_x, sim_y, sim_z, status in paths:
        px = np.array(sim_x)  # Lateral -> X
        py = np.array(sim_y)  # Lateral -> Y
        pz = np.array(sim_z)  # Depth -> Z (Vertical)
        plot_paths.append((px, py, pz, status))
        all_px.append(px);
        all_py.append(py);
        all_pz.append(pz)

    all_px = np.concatenate(all_px)
    all_py = np.concatenate(all_py)
    all_pz = np.concatenate(all_pz)

    # Slightly wider figure to accommodate Z labels
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")

    x_span = np.linspace(all_px.min(), all_px.max(), 10)
    y_span = np.linspace(all_py.min(), all_py.max(), 10)
    XX, YY = np.meshgrid(x_span, y_span)

    ax.plot_surface(XX, YY, np.zeros_like(XX), alpha=0.15, color="grey", linewidth=0)
    ax.plot_surface(XX, YY, np.full_like(XX, THICKNESS), alpha=0.15, color="grey", linewidth=0)

    color_map = {"reflected": "tab:red", "transmitted": "tab:green", "absorbed": "tab:blue"}
    for px, py, pz, status in plot_paths:
        ax.plot(px, py, pz, color=color_map.get(status), linewidth=0.6, alpha=0.6)

    ax.set_xlabel(r"Lateral $x$ [cm]")
    ax.set_ylabel(r"Lateral $y$ [cm]")
    ax.set_zlabel(r"Depth $z$ [cm]")

    _set_equal_aspect_3d(ax, all_px, all_py, all_pz)
    ax.invert_zaxis()  # 0 at top, 0.1 at bottom
    ax.view_init(elev=20, azim=-60)
    _common_plot_elements(ax, rf"Photon Trajectories in Finite Slab (Top-Down View, $N={N_PHOTONS}$)")

    # Increase padding to fix label cutoff on Z-axis
    fig.tight_layout(pad=3.0)
    fig.savefig(output_path, bbox_inches='tight')
    print(f"Saved {output_path}")


if __name__ == "__main__":
    run_trajectory_sim()
    plot_left_to_right()
    plot_top_to_bottom()