"""Public operations that cross equilibrium and catalog boundaries."""

from __future__ import annotations

from .equilibrium import KW, acid_ph, base_ph
from .species import Species, find_species


def ph(
    species: Species | str,
    c: float,
    kw: float = KW,
    exact: bool = False,
) -> float:
    """Return the pH of a tabulated species."""
    selected = find_species(species) if isinstance(species, str) else species
    if selected.is_base:
        return base_ph(
            kb=selected.k,
            c=c,
            kw=kw,
            exact=exact,
        )
    return acid_ph(
        ka=selected.k,
        c=c,
        kw=kw,
        exact=exact,
    )
