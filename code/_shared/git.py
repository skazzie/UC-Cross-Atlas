"""Tiny git helper for reproducibility columns."""

from __future__ import annotations

import subprocess
from pathlib import Path


def git_sha() -> str:
    """Short git SHA of HEAD, or 'unknown' if not in a git tree."""
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=Path(__file__).resolve().parent,
            text=True,
        ).strip()
    except Exception:
        return "unknown"
