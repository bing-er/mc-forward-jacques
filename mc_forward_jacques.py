"""
mc_forward_jacques.py

Author: Binger Yu
Date: 2025-11-28
version: 1.1

Python translation of Jacques' Monte Carlo model (Chapter 5, Welch & van Gemert, 2011)
<https://omlc.org/software/mc/Jacques2011_MonteCarlo_Welch&VanGemert.pdf>
and the MCML standard (Wang et al., 1995)<https://omlc.org/software/mc/man_mcml.pdf>.

REFERENCES:
1. Jacques, S. L. "Monte Carlo Modeling..." in Optical-Thermal Response of Laser-Irradiated Tissue.
   - Eq 5.34-5.35: Boundary interaction (Hop/Step size).
   - Eq 5.36: Fresnel Reflection.
2. Wang, L., Jacques, S. L., & Prahl, S. A. (1995). MCML Manual.
   - Chapter 3: Rules for Photon Propagation.

KEY FEATURES:
- Geometry: Finite slab (0 <= z <= d).
- Source: Collimated pencil beam (Green's function impulse response).
- Physics: Fresnel reflection/refraction at boundaries.
"""

import math
from typing import Tuple
import numpy as np

# ---- global constants ----
PI = math.pi
ALIVE = 1
DEAD = 0
THRESHOLD = 0.01  # roulette threshold
CHANCE = 0.1  # roulette survival probability
ONE_MINUS_COSZERO = 1.0e-12
STEP_EPS = 1.0e-9  # safety for log(0)


def sign(x: float) -> int:
    """Return the sign of x as +1 or -1."""
    return 1 if x >= 0.0 else -1


def fresnel_reflectance(n1: float, n2: float, ca1: float) -> float:
    """
    Calculate Fresnel reflectance R for unpolarized light.
    Matches Eq. 5.36 (Jacques 2011) / Eq. 3.26 (MCML Manual),
    but uses efficient cosine formulation to avoid expensive trig functions.

    Parameters
    ----------
    n1, n2 : Refractive indices of current and entered medium
    ca1    : Cosine of the angle of incidence (0 <= ca1 <= 1)
    """
    if n1 == n2: return 0.0

    if ca1 > 1.0: ca1 = 1.0
    if ca1 < 0.0: ca1 = 0.0

    sa1_sq = 1.0 - ca1 * ca1

    # Snell's Law: n1 * sin(theta1) = n2 * sin(theta2)
    # Check for Total Internal Reflection (TIR)
    val = (n1 / n2) * (n1 / n2) * sa1_sq
    if val >= 1.0:
        return 1.0  # TIR

    ca2 = math.sqrt(1.0 - val)  # Cosine of transmission angle

    # Fresnel formulas (Amplitude Coefficients)
    # This is mathematically equivalent to the tan/sin formula in Eq 5.36
    # but faster computationally.
    rs = (n1 * ca1 - n2 * ca2) / (n1 * ca1 + n2 * ca2)
    rp = (n1 * ca2 - n2 * ca1) / (n1 * ca2 + n2 * ca1)

    # Unpolarized light is average of s and p polarizations
    return 0.5 * (rs ** 2 + rp ** 2)


def run_mc(
        mua: float = 1.673,      # absorption coeff [cm^-1]
        mus: float = 312.0,      # scattering coeff [cm^-1]
        g: float = 0.90,         # anisotropy
        nt: float = 1.33,        # refractive index of tissue (slab)
        n_env: float = 1.0,      # refractive index of environment (air)
        thickness: float = 0.1,  # thickness of the slab [cm]
        Nphotons: int = 10_000,
        radial_size: float = 2.0,
        NR: int = 500,
        seed: int | None = 1234,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """
    Run MC simulation for a finite slab (Pencil Beam).
    """

    rng = np.random.default_rng(seed)

    # Optical Properties
    mut = mua + mus
    if mut <= 0: raise ValueError("mua + mus must be > 0")
    absorption_rate = mua / mut
    dr = radial_size / NR

    # Scoring Arrays
    Csph = np.zeros(NR + 1, dtype=float)
    Ccyl = np.zeros(NR + 1, dtype=float)
    Cpla = np.zeros(NR + 1, dtype=float)

    # Counters for macroscopic results
    total_reflection = 0.0
    total_transmission = 0.0

    # Specular Reflection at surface (normal incidence)
    # See Eq. 5.36 for theta=0 (Normal incidence)
    r_sp = (n_env - nt) / (n_env + nt)
    specular_R = r_sp ** 2

    # Main Loop
    i_photon = 0
    while i_photon < Nphotons:
        i_photon += 1

        # LAUNCH: Collimated Beam (straight down)
        x, y, z = 0.0, 0.0, 0.0
        ux, uy, uz = 0.0, 0.0, 1.0

        # Subtract specular reflection immediately
        W = 1.0 - specular_R

        if W <= 0: continue

        photon_status = ALIVE

        while photon_status == ALIVE:
            # 1. HOP: Sample step size
            # Eq 5.31 (Jacques 2011) / Eq 3.13 (MCML)
            rnd = STEP_EPS + (1.0 - STEP_EPS) * rng.random()
            s = -math.log(rnd) / mut

            # 2. CHECK BOUNDARIES (Eq 5.33 - 5.35)
            # We must loop because a reflected photon might hit the other
            # boundary before the step 's' is finished.
            while s > 0:
                d_bound = 0.0

                # Distance to boundary (Eq 3.23 MCML / Eq 5.34 Jacques)
                if uz > 0:
                    d_bound = (thickness - z) / uz  # To bottom (z=d)
                elif uz < 0:
                    d_bound = -z / uz  # To top (z=0)
                else:
                    d_bound = float('inf')  # Horizontal

                if s > d_bound:
                    # Hit boundary! Move to boundary first.
                    x += d_bound * ux
                    y += d_bound * uy
                    z += d_bound * uz
                    s -= d_bound  # Remaining step

                    # Interaction: Calculate Fresnel Reflection (Eq 5.36)
                    Refl = fresnel_reflectance(nt, n_env, abs(uz))

                    if rng.random() > Refl:
                        # Transmit (Escape) - Eq 5.37
                        photon_status = DEAD
                        if uz < 0:
                            total_reflection += W  # Escaped Top
                        else:
                            total_transmission += W  # Escaped Bottom
                        break
                    else:
                        # Reflect (Bounce back)
                        uz = -uz
                        # Photon continues remaining 's' inside this loop
                else:
                    # No boundary hit. Move full step.
                    x += s * ux
                    y += s * uy
                    z += s * uz
                    s = 0.0

            if photon_status == DEAD:
                break

            # 3. DROP (Absorb)
            absorbed = W * absorption_rate
            W -= absorbed

            # Scoring
            r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
            ir = int(r / dr)
            if ir > NR: ir = NR
            Csph[ir] += absorbed

            r_cyl = math.sqrt(x ** 2 + y ** 2)
            ir_cyl = int(r_cyl / dr)
            if ir_cyl > NR: ir_cyl = NR
            Ccyl[ir_cyl] += absorbed

            r_pla = abs(z)
            ir_pla = int(r_pla / dr)
            if ir_pla > NR: ir_pla = NR
            Cpla[ir_pla] += absorbed

            # 4. SPIN (Scatter)
            # Sample Henyey-Greenstein (Eq 3.19 in Jacques 2011)
            rnd = rng.random()
            if g == 0.0:
                costheta = 2.0 * rnd - 1.0
            else:
                temp = (1.0 - g ** 2) / (1.0 - g + 2.0 * g * rnd)
                costheta = (1.0 + g ** 2 - temp ** 2) / (2.0 * g)

            costheta = max(min(costheta, 1.0), -1.0)
            sintheta = math.sqrt(max(0.0, 1.0 - costheta ** 2))

            psi = 2.0 * PI * rng.random()
            cospsi = math.cos(psi)
            if psi < PI:
                sinpsi = math.sqrt(max(0.0, 1.0 - cospsi ** 2))
            else:
                sinpsi = -math.sqrt(max(0.0, 1.0 - cospsi ** 2))

            if 1.0 - abs(uz) <= ONE_MINUS_COSZERO:
                uxx = sintheta * cospsi
                uyy = sintheta * sinpsi
                uzz = costheta * sign(uz)
            else:
                temp = math.sqrt(1.0 - uz ** 2)
                uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta
                uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta
                uzz = -sintheta * cospsi * temp + uz * costheta

            ux, uy, uz = uxx, uyy, uzz

            # 5. ROULETTE
            if W < THRESHOLD:
                if rng.random() <= CHANCE:
                    W /= CHANCE
                else:
                    photon_status = DEAD

    # Normalization
    r_centers = np.arange(NR + 1) * dr + (dr / 2.0)

    vol_sph = 4.0 * PI * (r_centers ** 2) * dr
    Fsph = Csph / Nphotons / vol_sph / mua

    vol_cyl = 2.0 * PI * r_centers * dr
    Fcyl = Ccyl / Nphotons / vol_cyl / mua

    vol_pla = dr
    Fpla = Cpla / Nphotons / vol_pla / mua

    R_total = total_reflection / Nphotons
    T_total = total_transmission / Nphotons

    return r_centers, Fsph, Fcyl, Fpla, R_total, T_total


if __name__ == "__main__":
    # Test run
    r, Fsph, Fcyl, Fpla, R_tot, T_tot = run_mc(
        Nphotons=10_000,
        mua=10.0,
        mus=90.0,
        g=0.9,
        nt=1.33,
        n_env=1.0,
        thickness=0.02,
        radial_size=0.05,
        NR=100
    )
    print(f"Simulation Complete. R={R_tot:.4f}, T={T_tot:.4f}")