"""Install a built wheel in an isolated environment and exercise its real CLI.

Uses temporary test data only. Run from repository root:
  python tests/check_installed_wheel.py dist/rrmg_fungal_spikein-*.whl
All platforms exercise the real Python collectors with missing inputs and
two-row synthetic fixtures. Linux also checks the original runner with mock
stage wrappers and the combined shell collector. No scheduler is called.
Requires pip >=22.3 on the interpreter running this test.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
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
        figure_assets = ("runbook.sh", "RUNBOOK_FIGURES_README.txt", "rrmg_windows_r.sh",
                         "rrmg_plot_common.R", "rrmg_plot_preview.R", "rrmg_plot_alignment.R",
                         "rrmg_prepare_alignment_inputs.py")
        for asset in figure_assets:
            assert "rrmg_cli/_workflow/" + asset in names, "Missing figure asset: " + asset
        assert "rrmg_cli/_workflow/rrmg_mix_paired_fastq.py" in names
        assert not any(n.endswith("/rrmg_config.sh") for n in names)
        assert not any("/_workflow/temp/" in n or "/_workflow/input/" in n for n in names)
        expected = z.read("rrmg_cli/_workflow/VERSION").decode().strip()
    with tempfile.TemporaryDirectory(prefix="rrmg-wheel-test-") as tmp:
        root = Path(tmp).resolve()
        envdir = root / "env"
        # Install through the caller's pip so Debian/Ubuntu Python without
        # ensurepip can also test a genuinely isolated, pip-free environment.
        venv.EnvBuilder(with_pip=False).create(envdir)
        bindir = envdir / ("Scripts" if os.name == "nt" else "bin")
        python = bindir / ("python.exe" if os.name == "nt" else "python")
        env = dict(os.environ, PYTHONPATH="", PYTHONNOUSERSITE="1")
        env.pop("RRMG_CONFIG", None)
        subprocess.run([sys.executable, "-m", "pip", "--python", str(python),
                        "install", "--no-index", "--no-deps", str(wheel)], check=True, env=env, cwd=root)
        exe = bindir / ("rrmg.exe" if os.name == "nt" else "rrmg")
        def run(*parts, ok=True, error=None):
            result = subprocess.run([str(exe), *map(str, parts)], cwd=root, env=env, text=True, capture_output=True)
            if ok:
                assert result.returncode == 0, result.stdout + result.stderr
            else:
                assert result.returncode != 0, "Expected refusal"
            if error is not None:
                assert error in result.stdout + result.stderr, result.stdout + result.stderr
            return result.stdout
        assert run("--version").strip() == "rrmg " + expected
        installed = Path(run("path").strip())
        assert envdir in installed.parents, installed
        export = root / "exported core"
        run("export", export)
        for asset in figure_assets:
            assert (export / asset).is_file(), "Missing exported figure asset: " + asset
        runbook = (export / "runbook.sh").read_text(encoding="utf-8")
        assert "RUNBOOK_DIR" not in runbook
        assert 'FIG_OUT="$WD/figure_output/local_preview"' in runbook
        assert 'ALIGN_INPUT="$WD/figure_data/CYP51A_alignment"' in runbook
        assert 'ALIGN_FIG_OUT="$WD/figure_output/alignment_sites"' in runbook
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
        # Collectors must reject missing results, not produce an empty success.
        # Test this on Windows as well, rather than hiding it in the Linux branch.
        run_root = wd / "temp/rrmg_permutation"
        collectors = [
            ("rrmg_collect_reads_to_gene_depth.py",
             "readmap/reads_to_GENE.coverage_summary.tsv",
             "rrmg_reads_to_gene_depth_summary.tsv", "No readmap coverage summaries",
             [{"RunID": "run1", "RowIndex": "1", "MeanDepth": "2.500000"},
              {"RunID": "run2", "RowIndex": "2", "MeanDepth": "0.000000"}]),
            ("rrmg_collect_afum_genome_cov.py",
             "target_readmap/afum_target_reads_to_genome.summary.tsv",
             "rrmg_afum_target_reads_genome_coverage_summary.tsv", "No target_readmap summaries",
             [{"RunID": "run1", "RowIndex": "1", "ObservedMeanDepth": "1.250000"},
              {"RunID": "run2", "RowIndex": "2", "ObservedMeanDepth": "0.000000"}]),
        ]
        def assert_table(path, expected_rows):
            with path.open(newline="") as fh:
                assert list(csv.DictReader(fh, delimiter="\t")) == expected_rows, path
        for script, relative, output_name, error, records in collectors:
            out = root / "python collector outputs" / output_name
            run("script", script, "--", "--root", run_root, "--out", out, ok=False, error=error)
            assert not out.exists(), "Missing-input test must not create a success output"
            # Only tiny synthetic TSVs, never project/HPC data.
            for record in records:
                path = run_root / record["RunID"] / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", newline="") as fh:
                    writer = csv.DictWriter(fh, fieldnames=list(record), delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerow(record)
            run("script", script, "--", "--root", run_root, "--out", out)
            assert_table(out, records)
        print("[OK] Installed Python collectors reject missing inputs and preserve two fixture rows")
        if os.name != "nt":
            pipeline = root / "mock pipeline"
            pipeline.mkdir()
            (pipeline / "rrmg_run_local.sh").write_bytes((export / "rrmg_run_local.sh").read_bytes())
            wrapper = pipeline / "rrmg_run_one.sh"
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
            empty_wd = root / "empty project"
            empty_wd.mkdir()
            run("script", "rrmg_collect_supplementary_downstream.sh", "--", empty_wd,
                ok=False, error="No readmap coverage summaries")
            assert not (empty_wd / "result/rrmg_reads_to_gene_depth_summary.tsv").exists()
            run("script", "rrmg_collect_supplementary_downstream.sh", "--", wd)
            for _, _, output_name, _, records in collectors:
                assert_table(wd / "result" / output_name, records)
            # Also exercise the explicit absolute ROOT argument, with spaces.
            run("script", "rrmg_collect_supplementary_downstream.sh", "--", wd, run_root)
            for _, _, output_name, _, records in collectors:
                assert_table(wd / "result" / output_name, records)
            print("[OK] Linux original-runner row/failure dispatch and separate-directory collectors")
        else:
            print("[NOTE] Linux execution tests skipped on Windows; covered by Linux CI configuration")
        print("[OK] Installed wheel CLI, bundled assets, export hashes, init, no overwrite, Python utility")


if __name__ == "__main__":
    main()
