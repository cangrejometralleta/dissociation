"""Approximate and exact buffer calculations."""

from __future__ import annotations

import math

from .equilibrium import KW, _betas, _solve_charge_balance, pk


def henderson_hasselbalch(
    ka: float,
    acid_c: float,
    base_c: float,
) -> float:
    """Return the Henderson-Hasselbalch approximation."""
    if acid_c <= 0 or base_c <= 0:
        raise ValueError("a buffer needs both members of the pair")
    return pk(ka) + math.log10(base_c / acid_c)


def buffer_h(
    ka: float,
    acid_c: float,
    base_c: float,
    kw: float = KW,
) -> float:
    """Return exact [H+] for a weak acid and its conjugate base."""
    if ka <= 0:
        raise ValueError("Ka must be > 0")
    if acid_c < 0 or base_c < 0:
        raise ValueError("concentrations must be >= 0")
    return _solve_charge_balance(
        betas=_betas([ka]),
        c=acid_c + base_c,
        kw=kw,
        spectator=base_c,
    )


def buffer_ph(
    ka: float,
    acid_c: float,
    base_c: float,
    kw: float = KW,
) -> float:
    """Return the exact pH of a buffer."""
    h = buffer_h(
        ka=ka,
        acid_c=acid_c,
        base_c=base_c,
        kw=kw,
    )
    return pk(h)
