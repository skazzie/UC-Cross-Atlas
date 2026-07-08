"""Parser for LDSC `.log` output — extract intercept, ratio, h2.

LDSC's `--h2` mode writes a multi-section text log with values like:

    Total Observed scale h2: 0.0712 (0.0034)
    Lambda GC: 1.1721
    Mean Chi^2: 1.2746
    Intercept: 1.0231 (0.0091)
    Ratio: 0.084 (0.033)

The parser is regex-based and tolerant of minor whitespace + format
variation across LDSC versions. Returns a dict with parsed numeric
fields and a `status` flag: "ok" if intercept + ratio both present,
"partial" if only one, "missing" if neither (LDSC crashed mid-run).

Importable for tests in `test_parse_log.py`.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional


_RE_VAL_SE = re.compile(
    r"^\s*(?P<label>[A-Za-z][\w\s/^.]*?)\s*:\s*"
    r"(?P<value>-?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)"
    r"(?:\s*\(\s*(?P<se>-?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s*\))?\s*$"
)


_LABEL_TO_KEY = {
    "Total Observed scale h2": "h2",
    "Total Liability scale h2": "h2",
    "Intercept": "intercept",
    "Ratio": "ratio",
    "Lambda GC": "lambda_gc",
    "Mean Chi^2": "mean_chisq",
    "Max Chi^2": "max_chisq",
}


def parse(log_path: Path) -> dict:
    """Parse one LDSC `.log` file. Returns dict with parsed fields + status.

    Always returns a dict (never raises on missing fields); inspect
    `status` to see whether the parse was complete.
    """
    out: dict[str, Optional[float]] = {
        "h2": None, "h2_se": None,
        "intercept": None, "intercept_se": None,
        "ratio": None, "ratio_se": None,
        "lambda_gc": None, "mean_chisq": None, "max_chisq": None,
    }
    if not Path(log_path).exists():
        out["status"] = "missing-file"
        return out

    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        m = _RE_VAL_SE.match(line)
        if not m:
            continue
        label = m.group("label").strip()
        key = _LABEL_TO_KEY.get(label)
        if key is None:
            continue
        try:
            out[key] = float(m.group("value"))
            se = m.group("se")
            if se is not None and f"{key}_se" in out:
                out[f"{key}_se"] = float(se)
        except (TypeError, ValueError):
            continue

    have_int = out["intercept"] is not None
    have_rat = out["ratio"] is not None
    if have_int and have_rat:
        out["status"] = "ok"
    elif have_int or have_rat:
        out["status"] = "partial"
    else:
        out["status"] = "missing"
    return out


def in_band_for_f10(intercept: Optional[float], ratio: Optional[float]) -> bool:
    """F10 acceptance: intercept <= 1.10 OR ratio <= 0.20.

    Returns False if both are None — F10 cannot be claimed without
    evidence.
    """
    if intercept is not None and intercept <= 1.10:
        return True
    if ratio is not None and ratio <= 0.20:
        return True
    return False
