"""Common acid and base constants, names, and aliases."""

from __future__ import annotations

import unicodedata
from dataclasses import dataclass, field
from importlib import resources
from typing import Any

import yaml

from .equilibrium import pk


@dataclass(frozen=True)
class Species:
    name: str
    formula: str
    k: tuple[float, ...]
    is_base: bool = False
    alias: tuple[str, ...] = field(default=())

    @property
    def pk(self) -> tuple[float, ...]:
        return tuple(pk(value) for value in self.k)

    def __repr__(self) -> str:
        label = "Kb" if self.is_base else "Ka"
        values = ", ".join(f"{value:.3g}" for value in self.k)
        return f"<{self.name} ({self.formula}) {label}=[{values}]>"


def _load_catalog() -> dict[str, Any]:
    catalog = resources.files("dissociation").joinpath("constants/common-species.yaml")
    with catalog.open(encoding="utf-8") as stream:
        loaded = yaml.safe_load(stream)
    if not isinstance(loaded, dict):
        raise ValueError("common-species.yaml must contain a mapping")
    return loaded


def _build_table(
    catalog: dict[str, Any],
    group: str,
    is_base: bool,
) -> dict[str, Species]:
    records = catalog.get(group)
    if not isinstance(records, list):
        raise ValueError(f"common-species.yaml must contain a {group} list")

    table: dict[str, Species] = {}
    for record in records:
        if not isinstance(record, dict):
            raise ValueError(f"every {group} entry must be a mapping")
        constants = record.get("constants")
        aliases = record.get("aliases", {})
        if not isinstance(constants, list) or not isinstance(aliases, dict):
            raise ValueError(f"every {group} entry needs a constants list and aliases mapping")

        flattened_aliases: list[str] = []
        for language, translated_names in aliases.items():
            if language not in {"en", "es", "pt"}:
                raise ValueError(f"unsupported alias language: {language}")
            if not isinstance(translated_names, list):
                raise ValueError(f"aliases.{language} must be a list")
            for translated_name in translated_names:
                name = str(translated_name)
                if name not in flattened_aliases:
                    flattened_aliases.append(name)

        species = Species(
            name=str(record["name"]),
            formula=str(record["formula"]),
            k=tuple(float(value) for value in constants),
            is_base=is_base,
            alias=tuple(flattened_aliases),
        )
        if species.formula in table:
            raise ValueError(f"duplicate species formula: {species.formula}")
        table[species.formula] = species
    return table


_CATALOG = _load_catalog()
ACIDS = _build_table(
    catalog=_CATALOG,
    group="acids",
    is_base=False,
)
BASES = _build_table(
    catalog=_CATALOG,
    group="bases",
    is_base=True,
)


def _normalize(value: str) -> str:
    decomposed = unicodedata.normalize("NFD", value)
    stripped = "".join(
        character for character in decomposed
        if not unicodedata.combining(character)
    )
    return stripped.strip().lower()


def find_species(key: str) -> Species:
    """Find a species by formula, name, or alias."""
    target = _normalize(key)
    for table in (ACIDS, BASES):
        for species in table.values():
            candidates = {
                _normalize(species.formula),
                _normalize(species.name),
                *(_normalize(alias) for alias in species.alias),
            }
            if target in candidates:
                return species
    raise KeyError(f"no constants available for {key!r}")
