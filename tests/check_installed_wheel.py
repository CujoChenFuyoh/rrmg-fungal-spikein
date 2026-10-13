"""Install a built wheel in an isolated environment and exercise its real CLI.

Uses temporary test data only. Run from repository root:
  python tests/check_installed_wheel.py dist/rrmg_fungal_spikein-0.2.0-py3-none-any.whl
Linux additionally exercises the original runner with mock stage wrappers and
the real supplementary collectors on an empty fixture. No scheduler is called.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import venv
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    args = parser.parse_args()
    wheel = args.wheel.resolve()
    with zipfile.ZipFile(wheel) as z:
        names = z.namelist()
        assert any(n.endswith("_workflow/rrmg_run_local.sh") for n in names)
        assert any(n.endswith("_workflow/docs/INSTALL.zh-CN.md") for n in names)
        assert not any(n.endswith("/rrmg_config.sh") for n in names)
        assert not any("/_workflow/temp/" in n or "/_workflow/input/" in n for n in names)
        expected = z.read("rrmg_cli/_workflow/VERSION").decode().strip()
    with tempfile.TemporaryDirectory(prefix="rrmg-wheel-test-") as tmp:
        root = Path(tmp).resolve()
        envdir = root / "env"
        venv.EnvBuilder(with_pip=True).create(envdir)
        bindir = envdir / ("Scripts" if os.name == "nt" else "bin")
        python = bindir / ("python.exe" if os.name == "nt" else "python")
        env = dict(os.environ, PYTHONPATH="", PYTHONNOUSERSITE="1")
        env.pop("RRMG_CONFIG", None)
        subprocess.run([str(python), "-m", "pip", "install", "--no-index", "--no-deps", str(wheel)], check=True, env=env, cwd=root)
        exe = bindir / ("rrmg.exe" if os.name == "nt" else "rrmg")
        def run(*parts, ok=True):
            result = subprocess.run([str(exe), *map(str, parts)], cwd=root, env=env, text=True, capture_output=True)
            if ok:
                assert result.returncode == 0, result.stdout + result.stderr
            else:
                assert result.returncode != 0, "Expected refusal"
            return result.stdout
        assert run("--version").strip() == "rrmg " + expected
        installed = Path(run("path").strip())
        assert envdir in installed.parents, installed
        export = root / "exported core"
        run("export", export)
        manifest = json.loads((export / "rrmg_core_manifest.json").read_text())
        for name, digest in manifest["sha256"].items():
            assert hashlib.sha256((export / name).read_bytes()).hexdigest() == digest
        run("export", export, ok=False)
        wd = root / "existing data"
        wd.mkdir()
        marker = wd / "original_result.tsv"
        marker.write_text("preserve\n")
        run("init", "--workdir", wd, "--pipeline-dir", export)
        before = (wd / "rrmg_config.sh").read_bytes()
        run("init", "--workdir", wd, ok=False)
        assert (wd / "rrmg_config.sh").read_bytes() == before
        assert marker.read_text() == "preserve\n"
        run("init", "--workdir", root / "snapshot project")
        assert (root / "snapshot project/.rrmg" / ("core-" + expected) / "rrmg_core_manifest.json").is_file()
        assert "usage" in run("script", "rrmg_make_plan.py", "--", "--help").lower()
        if os.name != "nt":
            pipeline = root / "mock pipeline"
            pipeline.mkdir()
            (pipeline / "rrmg_run_local.sh").write_bytes((export / "rrmg_run_local.sh").read_bytes())
            wrapper = pipeline / "rrmg_step2_target_merge_array.pbs"
            wrapper.write_text('#!/bin/bash\nprintf "%s\\n" "$RRMG_ROW_INDEX" >> "$TEST_CALLS"\n')
            plan = wd / "plan.tsv"
            plan.write_text("run_index\trun_id\n1\trun1\n2\trun2\n3\trun3\n")
            config = root / "mock.sh"
            config.write_text("PIPELINE_DIR=" + shlex.quote(str(pipeline)) + "\nWD=" + shlex.quote(str(wd)) +
                              "\nPLAN=" + shlex.quote(str(plan)) + "\nLOCAL_JOBS=1\n")
            env["TEST_CALLS"] = str(root / "calls.txt")
            run("run", "step2", "2-3", "--config", config, "--jobs", "1")
            assert (root / "calls.txt").read_text().splitlines() == ["2", "3"]
            run("run", "step2", "4", "--config", config, ok=False)
            wrapper.write_text("#!/bin/bash\nexit 9\n")
            run("run", "step2", "1", "--config", config, ok=False)
            env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")
            run("script", "rrmg_collect_supplementary_downstream.sh", "--", wd)
            assert (wd / "result/rrmg_reads_to_gene_depth_summary.tsv").is_file()
            assert (wd / "result/rrmg_afum_target_reads_genome_coverage_summary.tsv").is_file()
            print("[OK] Linux original-runner row/failure dispatch and separate-directory collectors")
        else:
            print("[NOTE] Linux execution tests skipped on Windows; covered by Linux CI configuration")
        print("[OK] Installed wheel CLI, bundled assets, export hashes, init, no overwrite, Python utility")


if __name__ == "__main__":
    main()
