"""Independent regression test for the general-N transversality theorem.

The test uses only the general polynomial, not the rod or beam implementation.
For conservative real local models with N=2,...,5 it checks that the finite
root through c=s=(-1)^n has the slope of Eq. (6.5),

    c_s'(w_B) = d_1'(w_B) / sum_{j=2}^N d_j(w_B)/(C_j(w_B)-s),

where d_j=i*a_j*S_j and d_1' is given by Eq. (6.4).  It also verifies that the
root enters (-1,1) on the side predicted by that slope, repeats the slope test
on hosts with sum_j d_j = 0 but T_B != 0 (Eq. (6.3)), and checks the N=1
coincidence residue c(w_B) = s(1 + Gamma_B) of Eq. (6.7).

This is an independent algebraic regression of the general-N result of Section 6
(Section 6.2).  It is intentionally separate from the physical Timoshenko
implementation.

Produces the host counts and worst errors quoted in Section 6.2 (1000, 300 and 500 hosts).
Run from python/nchannel/:  python nchan_transversality.py
(NumPy; prints RESULT: PASS; output in nchan_transversality_output.txt)
"""

from __future__ import annotations

import numpy as np


RNG = np.random.default_rng(20260828)


def product_polynomial(roots: np.ndarray) -> np.ndarray:
    """Coefficients of prod_j (C_j-c), descending in c."""
    return (-1) ** len(roots) * np.poly(roots)


def normalized_coefficients(
    omega: float,
    omega_b: float,
    n: int,
    slope_k_l: float,
    other_phases: np.ndarray,
    alpha: np.ndarray,
    mass: float,
) -> np.ndarray:
    """Coefficients of F=K_T^{-1} prod(C_j-c)+sum d_j prod_{l!=j}(C_l-c), Eq. (6.2)."""
    delta = omega - omega_b
    phase_1 = n * np.pi + slope_k_l * delta
    phases = np.concatenate(([phase_1], other_phases))
    C = np.cos(phases)
    S = np.sin(phases)
    d = alpha * S  # choose a_j=-i*alpha_j, so d_j=i*a_j*S_j is real

    epsilon = (omega**2 - omega_b**2) / (mass * omega_b**2 * omega**2)
    coeff = epsilon * product_polynomial(C)
    for j in range(len(C)):
        others = np.delete(C, j)
        term = d[j] * product_polynomial(others)
        coeff = np.polyadd(coeff, term)
    return np.real_if_close(coeff, tol=1000).astype(float)


def finite_root(coeff: np.ndarray, target: float) -> float:
    roots = np.roots(coeff)
    root = roots[np.argmin(np.abs(roots - target))]
    if abs(root.imag) > 2e-7:
        raise AssertionError(f"expected a real finite root, obtained {root}")
    return float(root.real)


def check_transversality() -> tuple[float, float, int]:
    worst_slope_error = 0.0
    worst_imaginary_part = 0.0
    tested = 0
    omega_b = 5.0
    mass = 0.7
    h = 2.0e-6

    for N in range(2, 6):
        accepted = 0
        while accepted < 250:
            n = int(RNG.integers(1, 5))
            s = float((-1) ** n)
            slope_k_l = float(RNG.uniform(0.25, 2.0))

            # Keep the other waves away from a simultaneous Bragg degeneracy.
            other_phases = RNG.uniform(0.2, np.pi - 0.2, N - 1)
            if s < 0:
                other_phases = RNG.uniform(-1.25, 1.25, N - 1)
            C_other = np.cos(other_phases)
            if np.min(np.abs(C_other - s)) < 0.12:
                continue

            alpha = RNG.uniform(0.3, 2.0, N)
            d_other = alpha[1:] * np.sin(other_phases)
            trans_den = np.sum(d_other / (C_other - s))
            d1_prime = alpha[0] * s * slope_k_l
            predicted = d1_prime / trans_den
            if abs(trans_den) < 0.2 or abs(predicted) > 25:
                continue

            coeff_minus = normalized_coefficients(
                omega_b - h, omega_b, n, slope_k_l, other_phases, alpha, mass
            )
            coeff_plus = normalized_coefficients(
                omega_b + h, omega_b, n, slope_k_l, other_phases, alpha, mass
            )
            c_minus = finite_root(coeff_minus, s - predicted * h)
            c_plus = finite_root(coeff_plus, s + predicted * h)
            measured = (c_plus - c_minus) / (2 * h)
            rel_error = abs(measured - predicted) / max(1.0, abs(predicted))
            worst_slope_error = max(worst_slope_error, rel_error)

            # Select the side on which s*(c-s)<0; the root must be propagating.
            delta_pass = -np.sign(s * predicted) * h
            coeff_pass = normalized_coefficients(
                omega_b + delta_pass,
                omega_b,
                n,
                slope_k_l,
                other_phases,
                alpha,
                mass,
            )
            c_pass = finite_root(coeff_pass, s + predicted * delta_pass)
            if not (-1.0 < c_pass < 1.0):
                raise AssertionError(
                    f"N={N}: predicted pass-side root {c_pass} is outside (-1,1)"
                )

            roots_pass = np.roots(coeff_pass)
            worst_imaginary_part = max(
                worst_imaginary_part,
                float(np.min(np.abs(roots_pass.imag))),
            )
            accepted += 1
            tested += 1

    return worst_slope_error, worst_imaginary_part, tested


def check_without_branch_count() -> tuple[float, int]:
    """Force sum_j d_j=0 while retaining T_B!=0 of Eq. (6.3) (possible for N>=3)."""
    worst_slope_error = 0.0
    tested = 0
    omega_b = 5.0
    mass = 0.7
    h = 2.0e-6

    for N in range(3, 6):
        accepted = 0
        while accepted < 100:
            n = int(RNG.integers(1, 5))
            s = float((-1) ** n)
            slope_k_l = float(RNG.uniform(0.25, 2.0))
            other_phases = RNG.uniform(0.25, np.pi - 0.25, N - 1)
            C_other = np.cos(other_phases)
            if np.min(np.abs(C_other - s)) < 0.12:
                continue

            alpha = np.empty(N)
            alpha[0] = RNG.uniform(0.3, 2.0)
            alpha[1:-1] = RNG.uniform(-2.0, 2.0, N - 2)
            alpha[-1] = -np.sum(
                alpha[1:-1] * np.sin(other_phases[:-1])
            ) / np.sin(other_phases[-1])
            if abs(alpha[-1]) > 8.0:
                continue

            d_other = alpha[1:] * np.sin(other_phases)
            if abs(np.sum(d_other)) > 2.0e-14:
                raise AssertionError("failed to construct the branch-count degeneracy")
            trans_den = np.sum(d_other / (C_other - s))
            d1_prime = alpha[0] * s * slope_k_l
            predicted = d1_prime / trans_den
            if abs(trans_den) < 0.2 or abs(predicted) > 25:
                continue

            c_minus = finite_root(
                normalized_coefficients(
                    omega_b - h, omega_b, n, slope_k_l, other_phases, alpha, mass
                ),
                s - predicted * h,
            )
            c_plus = finite_root(
                normalized_coefficients(
                    omega_b + h, omega_b, n, slope_k_l, other_phases, alpha, mass
                ),
                s + predicted * h,
            )
            measured = (c_plus - c_minus) / (2 * h)
            rel_error = abs(measured - predicted) / max(1.0, abs(predicted))
            worst_slope_error = max(worst_slope_error, rel_error)

            delta_pass = -np.sign(s * predicted) * h
            c_pass = finite_root(
                normalized_coefficients(
                    omega_b + delta_pass,
                    omega_b,
                    n,
                    slope_k_l,
                    other_phases,
                    alpha,
                    mass,
                ),
                s + predicted * delta_pass,
            )
            if not (-1.0 < c_pass < 1.0):
                raise AssertionError(
                    f"N={N}: branch-count-degenerate root {c_pass} missed (-1,1)"
                )

            accepted += 1
            tested += 1

    return worst_slope_error, tested


def check_single_wave_residue() -> tuple[float, float, int]:
    """N=1: compare the pole-zero residue of c = C_1 + i K_T a_1 S_1 with Eq. (6.7)."""
    worst_limit_error = 0.0
    smallest_margin = np.inf
    tested = 0
    omega_b = 5.0
    h = 1.0e-7

    for _ in range(500):
        n = int(RNG.integers(1, 5))
        s = float((-1) ** n)
        slope_k_l = float(RNG.uniform(0.2, 2.5))
        alpha = float(RNG.uniform(0.2, 2.0))  # i*a_1=alpha>0: passive
        mass = float(RNG.uniform(0.1, 1.5))
        gamma = alpha * mass * omega_b**3 * slope_k_l / 2.0
        predicted = s * (1.0 + gamma)

        values = []
        for delta in (-h, h):
            omega = omega_b + delta
            phase = n * np.pi + slope_k_l * delta
            C1 = np.cos(phase)
            S1 = np.sin(phase)
            Kt = mass * omega**2 * omega_b**2 / (omega**2 - omega_b**2)
            values.append(C1 + Kt * alpha * S1)
        measured = 0.5 * (values[0] + values[1])
        rel_error = abs(measured - predicted) / max(1.0, abs(predicted))
        worst_limit_error = max(worst_limit_error, rel_error)
        smallest_margin = min(smallest_margin, abs(predicted) - 1.0)
        if abs(predicted) <= 1.0:
            raise AssertionError("passive forward N=1 residue did not enter a stop band")
        tested += 1

    return worst_limit_error, float(smallest_margin), tested


def main() -> None:
    slope_error, imaginary, n_multi = check_transversality()
    degenerate_error, n_degenerate = check_without_branch_count()
    residue_error, margin, n_single = check_single_wave_residue()

    print("GENERAL-N TRANSVERSALITY REGRESSION")
    print(f"  N=2..5 conservative hosts tested : {n_multi}")
    print(f"  worst relative slope error       : {slope_error:.3e}")
    print(f"  worst selected-root imaginary    : {imaginary:.3e}")
    print("  pass-side test                   : every predicted root entered (-1,1)")
    print(f"  sum(d_j)=0, T_B!=0 cases tested  : {n_degenerate}")
    print(f"  worst slope error in those cases : {degenerate_error:.3e}")
    print(f"  passive forward N=1 cases        : {n_single}")
    print(f"  worst coincidence-residue error  : {residue_error:.3e}")
    print(f"  smallest N=1 stop-band margin    : {margin:.3e}")

    if slope_error > 2e-5 or degenerate_error > 2e-5 or residue_error > 2e-7:
        raise SystemExit("FAIL: theorem regression exceeded tolerance")
    print("  RESULT                           : PASS")


if __name__ == "__main__":
    main()
