"""Titration calculations for a weak acid and strong base."""

from __future__ import annotations

from typing import Sequence

from .equilibrium import KW, _betas, _solve_charge_balance, pk


def equivalence_volume(
    analyte_c: float,
    analyte_volume: float,
    titrant_c: float,
    protons: int = 1,
) -> float:
    """Return the titrant volume that neutralizes the analyte."""
    if titrant_c <= 0:
        raise ValueError("the titrant must have a concentration > 0")
    return protons * analyte_c * analyte_volume / titrant_c


def titration_ph(
    kas: Sequence[float],
    analyte_c: float,
    analyte_volume: float,
    titrant_c: float,
    titrant_volume: float,
    kw: float = KW,
) -> float:
    """Return the pH at one point of a strong-base titration."""
    total_volume = analyte_volume + titrant_volume
    if total_volume <= 0:
        raise ValueError("the total volume must be > 0")
    diluted_analyte = analyte_c * analyte_volume / total_volume
    spectator = titrant_c * titrant_volume / total_volume
    h = _solve_charge_balance(
        betas=_betas(kas),
        c=diluted_analyte,
        kw=kw,
        spectator=spectator,
    )
    return pk(h)


def titration_curve(
    kas: Sequence[float],
    analyte_c: float,
    analyte_volume: float,
    titrant_c: float,
    up_to: float | None = None,
    points: int = 101,
    kw: float = KW,
) -> list[tuple[float, float]]:
    """Return the whole curve as titrant-volume and pH pairs."""
    if points < 2:
        raise ValueError("a curve needs at least two points")
    if up_to is None:
        final_equivalence = equivalence_volume(
            analyte_c=analyte_c,
            analyte_volume=analyte_volume,
            titrant_c=titrant_c,
            protons=len(kas),
        )
        up_to = 2.0 * final_equivalence
    step = up_to / (points - 1)
    return [
        (
            i * step,
            titration_ph(
                kas=kas,
                analyte_c=analyte_c,
                analyte_volume=analyte_volume,
                titrant_c=titrant_c,
                titrant_volume=i * step,
                kw=kw,
            ),
        )
        for i in range(points)
    ]
