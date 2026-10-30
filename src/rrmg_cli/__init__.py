"""Packaging/dispatch only; scientific analysis stays in the original scripts."""
from importlib.metadata import version
from pathlib import Path

_source_version = Path(__file__).resolve().parents[2] / "VERSION"
if _source_version.is_file():
    # A checked-out source must not report an older installed distribution.
    __version__ = _source_version.read_text().strip()
else:
    __version__ = version("rrmg-fungal-spikein")
