"""Locked v1 constants shared across drivers.

Single source of truth — duplicating these in drivers loses the DECISIONS.md
provenance. Import from here, do not redefine.
"""

from __future__ import annotations

BOOTSTRAP_N: int = 1000
SEED: int = 42
MIN_CELLS_PER_TYPE: int = 50

FDR_PRIMARY: float = 0.05
FDR_SATURATION: float = 0.01  # used when >=80% saturation at FDR_PRIMARY in both atlases
KAPPA_SATURATION_THRESHOLD: float = 0.80

TOP_K_BROAD: tuple[int, ...] = (5, 10)
TOP_K_FINE: tuple[int, ...] = (5, 10, 20)

SCDRS_N_CTRL: int = 1000
SEISMIC_N_PERMUTATIONS: int = 1000

UC_ATLASES: tuple[str, ...] = ("smillie", "garrido_trigo", "taurus")
UC_GWAS: tuple[str, ...] = ("delange", "liu")
METHODS: tuple[str, ...] = ("scdrs", "seismic")
TIERS: tuple[str, ...] = ("broad", "fine")
