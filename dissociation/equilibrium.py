"""Core acid-base equilibrium calculations."""

from __future__ import annotations

import math
from typing import Sequence

KW = 1.0e-14


def pk(k: float) -> float:
    """Return pKa, pKb, or pH from a constant or concentration."""
    if k <= 0:
        raise ValueError("the constant must be > 0")
    return -math.log10(k)


def k_from_pk(value: float) -> float:
    """Return a constant from its pK."""
    return 10.0 ** (-value)


def _pkw(kw: float) -> float:
    return -math.log10(kw)


def monoprotic_h(ka: float, c: float) -> float:
    """Return exact [H+] for a weak monoprotic acid, ignoring water."""
    if ka <= 0 or c < 0:
        raise ValueError("Ka must be > 0 and C must be >= 0")
    if c == 0:
        return 0.0
    return 2.0 * ka * c / (ka + math.sqrt(ka * ka + 4.0 * ka * c))


def _classic_monoprotic_h(ka: float, c: float) -> float:
    """Evaluate the same formula directly from the quadratic."""
    return (ka / 2.0) * (math.sqrt(1.0 + 4.0 * c / ka) - 1.0)


def _betas(kas: Sequence[float]) -> list[float]:
    """Return cumulative constants beta_j, including beta_0 = 1."""
    products = [1.0]
    for ka in kas:
        if ka <= 0:
            raise ValueError("every Ka must be > 0")
        products.append(products[-1] * ka)
    return products


def fractions(kas: Sequence[float], h: float) -> list[float]:
    """Return the fraction alpha_j of each species at a given [H+]."""
    return _fractions_from_betas(_betas(kas), h)


def _fractions_from_betas(betas: Sequence[float], h: float) -> list[float]:
    n = len(betas) - 1
    terms = [betas[j] * h ** (n - j) for j in range(n + 1)]
    total = sum(terms)
    return [term / total for term in terms]


def _released_protons(
    betas: Sequence[float],
    c: float,
    h: float,
) -> float:
    alphas = _fractions_from_betas(betas, h)
    return c * sum(j * alpha for j, alpha in enumerate(alphas))


def _solve_charge_balance(
    betas: Sequence[float],
    c: float,
    kw: float,
    spectator: float = 0.0,
) -> float:
    """Close the charge balance by bisection on a logarithmic scale."""

    def residual(h: float) -> float:
        released = _released_protons(
            betas=betas,
            c=c,
            h=h,
        )
        return h + spectator - released - kw / h

    lo = 1e-20
    hi = max(1.0, c * (len(betas) - 1) + abs(spectator) + 1.0)
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if residual(mid) > 0:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def polyprotic_h(
    kas: Sequence[float],
    c: float,
    kw: float = KW,
) -> float:
    """Return exact [H+] for one or several constants, including water."""
    constants = list(kas)
    if not constants:
        raise ValueError("at least one Ka is required")
    if c < 0:
        raise ValueError("C must be >= 0")
    return _solve_charge_balance(
        betas=_betas(constants),
        c=c,
        kw=kw,
    )


def acid_ph(
    ka: float | Sequence[float],
    c: float,
    kw: float = KW,
    exact: bool = False,
) -> float:
    """Return the pH of a weak acid."""
    if isinstance(ka, (int, float)):
        if not exact:
            return pk(monoprotic_h(float(ka), c))
        kas: Sequence[float] = [float(ka)]
    else:
        kas = list(ka)
    h = polyprotic_h(
        kas=kas,
        c=c,
        kw=kw,
    )
    return pk(h)


def base_ph(
    kb: float | Sequence[float],
    c: float,
    kw: float = KW,
    exact: bool = False,
) -> float:
    """Return the pH of a weak base."""
    poh = acid_ph(
        ka=kb,
        c=c,
        kw=kw,
        exact=exact,
    )
    return _pkw(kw) - poh
