"""Worked examples used by the console entry point."""

from __future__ import annotations

from .api import ph
from .buffers import buffer_ph, henderson_hasselbalch
from .equilibrium import fractions
from .species import find_species
from .titration import titration_ph


def _demo_species() -> None:
    acetic_ph = ph(
        species="CH3COOH",
        c=0.1,
    )
    ammonia_ph = ph(
        species="NH3",
        c=0.1,
    )
    phosphoric_ph = ph(
        species="H3PO4",
        c=0.1,
    )
    hcn_ph = ph(
        species="HCN",
        c=1e-6,
        exact=True,
    )
    print("acetic acid 0.1 M        -> pH", round(acetic_ph, 3))
    print("ammonia 0.1 M            -> pH", round(ammonia_ph, 3))
    print("phosphoric acid 0.1 M    -> pH", round(phosphoric_ph, 3))
    print("HCN 1e-6 M (exact)       -> pH", round(hcn_ph, 3))
    print("H3PO4 speciation at pH 7 ->",
          [round(value, 4) for value in fractions(find_species("H3PO4").k, 1e-7)])


def _demo_buffer() -> None:
    acetic = find_species("CH3COOH").k[0]
    concentrated_ph = buffer_ph(
        ka=acetic,
        acid_c=0.1,
        base_c=0.1,
    )
    concentrated_approximation = henderson_hasselbalch(
        ka=acetic,
        acid_c=0.1,
        base_c=0.1,
    )
    dilute_ph = buffer_ph(
        ka=acetic,
        acid_c=1e-5,
        base_c=1e-5,
    )
    dilute_approximation = henderson_hasselbalch(
        ka=acetic,
        acid_c=1e-5,
        base_c=1e-5,
    )
    print("acetate buffer 1:1       -> pH", round(concentrated_ph, 3),
          "(Henderson-Hasselbalch:", round(concentrated_approximation, 3), ")")
    print("the same pair at 1e-5 M  -> pH", round(dilute_ph, 3),
          "(Henderson-Hasselbalch:", round(dilute_approximation, 3), ")")


def _demo_titration() -> None:
    acetic = find_species("CH3COOH").k[0]
    print("0.1 M acetic acid, 50 mL, titrated with 0.1 M NaOH:")
    for volume in (0.0, 25.0, 49.9, 50.0, 60.0):
        value = titration_ph(
            kas=[acetic],
            analyte_c=0.1,
            analyte_volume=50.0,
            titrant_c=0.1,
            titrant_volume=volume,
        )
        print(f"  {volume:5.1f} mL -> pH", round(value, 3))


def _demo_names() -> None:
    rows = (
        ("hydrofluoric acid", "ácido fluorhídrico", "ácido fluorídrico"),
        ("hydrocyanic acid", "ácido cianhídrico", "ácido cianídrico"),
        ("ammonia", "amoníaco", "amônia"),
    )
    for english, spanish, portuguese in rows:
        found = {find_species(key).formula for key in (english, spanish, portuguese)}
        print(f"  {english:18s} {spanish:20s} {portuguese:20s} -> {found.pop()}")


def demo() -> None:
    """Run the worked examples."""
    print("the same entry, three languages:")
    _demo_names()
    _demo_species()
    _demo_buffer()
    _demo_titration()
