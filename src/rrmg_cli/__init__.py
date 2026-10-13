"""Packaging/dispatch only; scientific analysis stays in the original scripts."""
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

try:
    __version__ = version("rrmg-fungal-spikein")
except PackageNotFoundError:
    __version__ = (Path(__file__).resolve().parents[2] / "VERSION").read_text().strip()
