"""Build the unchanged flat-layout workflow into wheel package resources.

Root scripts remain the only maintained source. No second hand-edited copy is
kept under src/. pip/build use this hook through the PEP 517 backend.
"""
from pathlib import Path
import shutil
from setuptools import setup
from setuptools.command.build_py import build_py

ROOT = Path(__file__).resolve().parent


class BuildWorkflow(build_py):
    def run(self):
        super().run()
        # Editable installs use the root sources via core_root(), not a stale copy.
        if self.editable_mode:
            return
        build_root = Path(self.build_lib).resolve()
        destination = build_root / "rrmg_cli" / "_workflow"
        # This directory contains generated resources only. A previous release's
        # build cache must not put removed workflow entries into a new wheel.
        if destination.is_symlink() or not destination.resolve().is_relative_to(build_root):
            raise RuntimeError("Unsafe workflow build destination: " + str(destination))
        if destination.exists():
            shutil.rmtree(destination)
        names = (ROOT / "src/rrmg_cli/core_files.txt").read_text(encoding="utf-8").splitlines()
        for name in names:
            if not name or name.startswith("#"):
                continue
            source = ROOT / name
            if not source.is_file() or source.is_symlink():
                raise RuntimeError("Missing/unsafe core source: " + name)
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            self.copy_file(str(source), str(target))


setup(cmdclass={"build_py": BuildWorkflow})
