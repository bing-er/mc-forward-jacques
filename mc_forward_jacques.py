"""
mc_forward_jacques.py

Author: Binger Yu
Date: 2025-11-22
version: 0.1

Python translation of Jacques' mc321.c (Chapter 5, Jacques 2011).

This implements a forward Monte Carlo simulation of photon transport
from an isotropic point source in an infinite, homogeneous medium.

Key features
-----------
- isotropic point source at the origin
- photon interactions: absorption, scattering, roulette termination
- no boundaries (photons never escape)
- absorbed energy scored in:
    * spherical shells  (r = sqrt(x^2 + y^2 + z^2))
    * cylindrical shells (r = sqrt(x^2 + y^2))
    * planar slabs      (r = |z|)
"""

import math
from typing import Tuple

import numpy as np

# ---- global constants (mirroring the C code) ----
PI = math.pi

ALIVE = 1
DEAD = 0

THRESHOLD = 0.01          # roulette threshold for photon weight
CHANCE = 0.1              # roulette survival probability

# numerical guard used when |uz| is very close to 1
ONE_MINUS_COSZERO = 1.0e-12


def sign(x: float) -> int:
    """Return the sign of x as +1 or -1."""
    return 1 if x >= 0.0 else -1


def run_mc(
    mua: float = 1.673,         # absorption coefficient [cm^-1]
    mus: float = 312.0,         # scattering coefficient [cm^-1]
    g: float = 0.90,            # anisotropy factor (Henyey–Greenstein)
    nt: float = 1.33,           # refractive index (not used here, kept for completeness)
    Nphotons: int = 10_000,     # number of launched photons
    radial_size: float = 2.0,   # maximum radius for scoring bins [cm]
    NR: int = 500,              # number of radial bins
    seed: int | None = 1234,    # RNG seed (None -> unpredictable)
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Run the forward Monte Carlo simulation (Jacques-style).

    Parameters
    ----------
    mua, mus : float
        Absorption and scattering coefficients [cm^-1].
    g : float
        Anisotropy parameter of the Henyey–Greenstein phase function (-1 <= g <= 1).
    nt : float
        Refractive index (unused in this infinite, boundary-free model).
    Nphotons : int
        Number of photons to simulate.
    radial_size : float
        Maximum radius covered by scoring bins [cm].
    NR : int
        Number of radial bins; an additional overflow bin NR is used.
    seed : int or None
        Seed for NumPy's default RNG.

    Returns
    -------
    r_centers : (NR+1,) ndarray
        Radial bin center positions [cm].
    Fsph, Fcyl, Fpla : (NR+1,) ndarray
        Fluence T(r) [1/cm^2] for spherical, cylindrical, and planar geometries.
        The last element (index NR) is the overflow bin.
    """
    rng = np.random.default_rng(seed)

    # ---- derived optical properties ----
    dr = radial_size / NR                    # radial bin width [cm]
    mut = mua + mus                          # total attenuation coefficient
    if mut <= 0:
        raise ValueError("mua + mus must be > 0")
    albedo = mus / mut                       # single-scattering albedo

    # absorbed energy accumulators (last bin is overflow)
    Csph = np.zeros(NR + 1, dtype=float)
    Ccyl = np.zeros(NR + 1, dtype=float)
    Cpla = np.zeros(NR + 1, dtype=float)

    # ------------------------------------------------------------------
    # main photon loop
    # ------------------------------------------------------------------
    i_photon = 0
    while i_photon < Nphotons:
        i_photon += 1

        # ---- LAUNCH: start at origin with isotropic direction ----
        W = 1.0
        photon_status = ALIVE

        x = y = z = 0.0

        # isotropic direction: cos(theta) and psi ~ U(0, 2π)
        costheta = 2.0 * rng.random() - 1.0
        sintheta = math.sqrt(max(0.0, 1.0 - costheta * costheta))
        psi = 2.0 * PI * rng.random()
        ux = sintheta * math.cos(psi)
        uy = sintheta * math.sin(psi)
        uz = costheta

        # ---- HOP–DROP–SPIN–ROULETTE loop for a single photon ----
        while photon_status == ALIVE:
            # ---------------- HOP: sample step length ----------------
            # Draw a strictly positive random number 0 < rnd <= 1
            rnd = 0.0
            while rnd <= 0.0:
                rnd = rng.random()
            s = -math.log(rnd) / mut   # step size [cm] (exponential free path)

            # update position
            x += s * ux
            y += s * uy
            z += s * uz

            # ---------------- DROP: absorb weight -------------------
            absorb = W * (1.0 - albedo)   # fraction absorbed in this step
            W -= absorb

            # Deposit absorbed energy into the three scoring geometries.

            # (1) Spherical shells: r = sqrt(x^2 + y^2 + z^2)
            r = math.sqrt(x * x + y * y + z * z)
            ir = int(r / dr)
            if ir >= NR:
                ir = NR
            Csph[ir] += absorb

            # (2) Cylindrical shells: r = sqrt(x^2 + y^2)
            r = math.sqrt(x * x + y * y)
            ir = int(r / dr)
            if ir >= NR:
                ir = NR
            Ccyl[ir] += absorb

            # (3) Planar slabs: r = |z|
            r = abs(z)
            ir = int(r / dr)
            if ir >= NR:
                ir = NR
            Cpla[ir] += absorb

            # ---------------- SPIN: scatter photon ------------------
            # Sample scattering angle from Henyey–Greenstein (HG).
            rnd = rng.random()
            if g == 0.0:
                # isotropic scattering
                costheta = 2.0 * rnd - 1.0
            else:
                # HG sampling consistent with mc321.c (no 3/2 exponent)
                temp = (1.0 - g * g) / (1.0 - g + 2.0 * g * rnd)
                costheta = (1.0 + g * g - temp * temp) / (2.0 * g)

            # Clamp for numerical safety
            costheta = max(min(costheta, 1.0), -1.0)
            sintheta = math.sqrt(max(0.0, 1.0 - costheta * costheta))

            # Azimuthal angle uniformly distributed
            psi = 2.0 * PI * rng.random()
            cospsi = math.cos(psi)

            # Compute sin(psi) with sign, like the original C code
            if psi < PI:
                sinpsi = math.sqrt(max(0.0, 1.0 - cospsi * cospsi))
            else:
                sinpsi = -math.sqrt(max(0.0, 1.0 - cospsi * cospsi))

            # Update direction cosines (ux, uy, uz)
            if 1.0 - abs(uz) <= ONE_MINUS_COSZERO:
                # Nearly parallel to z-axis: use simplified formulas
                uxx = sintheta * cospsi
                uyy = sintheta * sinpsi
                uzz = costheta * sign(uz)
            else:
                # General case: rotate around current direction
                temp = math.sqrt(max(0.0, 1.0 - uz * uz))
                uxx = (
                    sintheta * (ux * uz * cospsi - uy * sinpsi) / temp
                    + ux * costheta
                )
                uyy = (
                    sintheta * (uy * uz * cospsi + ux * sinpsi) / temp
                    + uy * costheta
                )
                uzz = -sintheta * cospsi * temp + uz * costheta

            ux, uy, uz = uxx, uyy, uzz

            # ------------- ROULETTE: terminate photon ---------------
            if W < THRESHOLD:
                # Survive with probability CHANCE; otherwise die.
                if rng.random() <= CHANCE:
                    W /= CHANCE   # boost weight to keep expectation unbiased
                else:
                    photon_status = DEAD

    # ------------------------------------------------------------------
    # POSTPROCESS: convert absorbed energy Csph/Ccyl/Cpla to fluence T(r)
    # using Jacques' normalization (Appendix 5.7).
    # ------------------------------------------------------------------
    r_centers = np.zeros(NR + 1, dtype=float)
    Fsph = np.zeros(NR + 1, dtype=float)
    Fcyl = np.zeros(NR + 1, dtype=float)
    Fpla = np.zeros(NR + 1, dtype=float)

    for ir in range(NR + 1):
        # bin center radius
        r = (ir + 0.5) * dr
        r_centers[ir] = r

        # (1) spherical shells: volume = 4π r^2 Δr
        shellvolume = 4.0 * PI * r * r * dr
        Fsph[ir] = Csph[ir] / Nphotons / shellvolume / mua

        # (2) cylindrical shells (per cm length): volume = 2π r Δr
        shellvolume = 2.0 * PI * r * dr
        Fcyl[ir] = Ccyl[ir] / Nphotons / shellvolume / mua

        # (3) planar slabs (per cm^2): thickness = Δr
        shellvolume = dr
        Fpla[ir] = Cpla[ir] / Nphotons / shellvolume / mua

    return r_centers, Fsph, Fcyl, Fpla


if __name__ == "__main__":
    # Example run with fewer photons for a quick sanity check
    r, Fsph, Fcyl, Fpla = run_mc(
        Nphotons=5_000,
        mua=1.0,
        mus=10.0,
        g=0.9,
    )

    # Print the first few spherical fluence values
    for i in range(5):
        print(f"r = {r[i]:.3f} cm,  Fsph = {Fsph[i]:.3e}")
