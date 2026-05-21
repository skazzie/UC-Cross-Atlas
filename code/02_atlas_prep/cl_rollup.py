"""Cell-type fine-tier roll-up.

High-cardinality author annotations (e.g. HCA Gut's ~120
``author_cell_type`` values) exceed the v1 "fine" tier target of ~30-50
categories. This module collapses those labels using:

  1. Strip marker-gene parentheticals: "Stromal 1 (ADAMDEC1+)" -> "Stromal 1"
  2. Strip trailing marker-suffix tokens after the base label
     (e.g. "M0 Ribhi" -> "M0"; "BEST4+ epithelial" preserved because the
     marker is the lineage identifier, not a sub-cluster suffix).
  3. Strip numeric / sub-cluster suffixes: "Neutrophil 1/2/3" -> "Neutrophil",
     "S1.2"/"S2a" -> "S1"/"S2".
  4. Per-atlas overrides loaded from ``cl_rollup_maps.yaml`` for edge
     cases the heuristics cannot handle.

The original ``cell_type_fine`` column is preserved. The new column is
added as ``cell_type_fine_rolled``.

Spec: DECISIONS.md correction 2026-05-20 (6/7); ``atlas_schemas.md``.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

from anndata import AnnData
import pandas as pd

logger = logging.getLogger(__name__)

_PAREN_RE = re.compile(r"\s*\([^()]*\)\s*$")
_TRAILING_SUFFIX_RE = re.compile(
    r"\s+(?:hi|lo|low|high|pos|neg|\+|-|[A-Z][a-z]*hi|[A-Z][a-z]*lo)$"
)
_NUMERIC_SUFFIX_RE = re.compile(r"\s*[0-9]+(?:\.[0-9]+)?[a-z]?$")
_LETTER_SUFFIX_RE = re.compile(r"([A-Z]+[0-9]+)[a-z]+$")  # S2a -> S2


def _strip_parenthetical(label: str) -> str:
    while True:
        new = _PAREN_RE.sub("", label).rstrip()
        if new == label:
            return label
        label = new


def _strip_letter_subcluster(label: str) -> str:
    m = _LETTER_SUFFIX_RE.fullmatch(label)
    if m:
        return m.group(1)
    return label


def _strip_numeric_suffix(label: str) -> str:
    # Apply to whole label first (e.g. "S1.2" -> "S1").
    label = _strip_letter_subcluster(label)
    new = _NUMERIC_SUFFIX_RE.sub("", label).rstrip()
    return new if new else label


def _strip_trailing_marker(label: str) -> str:
    """Remove a trailing "Ribhi"/"hi"/"lo"/"+"/"-" style suffix token.

    Only applied when there is at least one preceding whitespace-separated
    token (so single-token labels like "BEST4+" are preserved by callers
    that choose not to apply this step).
    """
    new = _TRAILING_SUFFIX_RE.sub("", label).rstrip()
    return new if new else label


def _heuristic_collapse(label: str) -> str:
    s = str(label).strip()
    s = _strip_parenthetical(s)
    s = _strip_numeric_suffix(s)
    if " " in s:
        s = _strip_trailing_marker(s)
    return s.strip()


def _load_overrides(atlas_slug: str) -> dict[str, str]:
    """Read per-atlas override mapping from ``cl_rollup_maps.yaml``.

    YAML structure:

        garrido_trigo: {}
        pangi:
          "Original Label": "Rolled Label"
        hca_gut:
          ...

    YAML is parsed with a very small custom reader so the package has no
    hard dependency on PyYAML. Missing file or missing section returns {}.
    """
    path = Path(__file__).with_name("cl_rollup_maps.yaml")
    if not path.exists():
        return {}
    try:
        import yaml  # type: ignore

        with path.open("r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
    except ImportError:
        data = _parse_simple_yaml(path)
    section = data.get(atlas_slug) or {}
    if not isinstance(section, dict):
        return {}
    return {str(k): str(v) for k, v in section.items()}


def _parse_simple_yaml(path: Path) -> dict[str, dict[str, str]]:
    """Tiny YAML reader for the limited cl_rollup_maps.yaml shape.

    Supports:
      - Top-level keys ``slug:`` opening a section.
      - Section entries ``  "Key": "Value"`` or ``  Key: Value``.
      - Empty sections (``slug: {}``).
      - ``#`` comments.
    """
    result: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    current_key: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line:
            continue
        if not line.startswith(" "):
            key, _, rest = line.partition(":")
            current_key = key.strip()
            rest = rest.strip()
            if rest in ("", "{}", "null"):
                current = {}
                result[current_key] = current
            continue
        if current is None or current_key is None:
            continue
        body = line.strip()
        if ":" not in body:
            continue
        k, _, v = body.partition(":")
        k = k.strip().strip('"').strip("'")
        v = v.strip().strip('"').strip("'")
        current[k] = v
    return result


def rollup_fine_tier(
    adata: AnnData,
    atlas_slug: str,
    target_n: int = 50,
) -> AnnData:
    """Add ``cell_type_fine_rolled`` to ``adata.obs``.

    Reads from ``cell_type_fine``. Heuristics first, then applies the
    per-atlas override map. Raises ``AssertionError`` if the resulting
    cardinality exceeds ``target_n``.
    """
    if "cell_type_fine" not in adata.obs.columns:
        raise KeyError(
            "rollup_fine_tier requires 'cell_type_fine' in obs; "
            "run the atlas loader first."
        )

    overrides = _load_overrides(atlas_slug)
    src = adata.obs["cell_type_fine"]

    def collapse(label: object) -> object:
        if pd.isna(label):
            return label
        s = str(label)
        if s in overrides:
            return overrides[s]
        rolled = _heuristic_collapse(s)
        return overrides.get(rolled, rolled)

    rolled = src.map(collapse)
    n_in = src.dropna().nunique()
    n_out = rolled.dropna().nunique()
    logger.info(
        "cl_rollup[%s]: %d -> %d categories (target <= %d)",
        atlas_slug,
        n_in,
        n_out,
        target_n,
    )

    adata = adata.copy()
    adata.obs["cell_type_fine_rolled"] = rolled.astype("category")
    assert n_out <= target_n, (
        f"rollup_fine_tier[{atlas_slug}]: post-rollup cardinality {n_out} "
        f"exceeds target_n={target_n}. Add overrides to cl_rollup_maps.yaml."
    )
    return adata
