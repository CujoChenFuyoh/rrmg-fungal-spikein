#!/usr/bin/env python3
"""Build a complete portable RRMG figure-data bundle without rerunning mapping."""

import argparse
from contextlib import contextmanager
import csv
import hashlib
import json
import shutil
import tempfile
import uuid
from pathlib import Path

from rrmg_prepare_alignment_inputs import build_alignment_bundle, input_file


OUTPUT_MARKER = ".rrmg_figure_data.json"


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def inventory(path):
    result = {}
    for item in sorted(Path(path).rglob("*")):
        if item.is_symlink():
            raise ValueError("Output symlinks are not supported: " + str(item))
        if item.is_file() and item.name != OUTPUT_MARKER:
            result[item.relative_to(path).as_posix()] = digest(item)
    return result


def check_output(path):
    path = Path(path)
    if path.is_symlink():
        raise ValueError("Output directory cannot be a symlink: " + str(path))
    if not path.exists():
        return
    marker = path / OUTPUT_MARKER
    if not marker.is_file():
        raise ValueError("Existing figure_data is unmanaged; preserve/move it aside first: " + str(path))
    saved = json.loads(marker.read_text(encoding="utf-8"))
    if saved.get("kind") != "rrmg_figure_data" or saved.get("sha256") != inventory(path):
        raise ValueError("figure_data contains added, changed or missing files; preserve/move it aside first: " + str(path))


@contextmanager
def managed_output(path):
    path = Path(path).absolute()
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = path.with_name(path.name + ".lock")
    try:
        lock.mkdir()
    except FileExistsError:
        raise ValueError("Output is locked: " + str(lock)) from None
    try:
        check_output(path)
        with tempfile.TemporaryDirectory(prefix="." + path.name + "-stage-", dir=path.parent) as tmp:
            stage = Path(tmp) / "output"
            stage.mkdir()
            yield stage
            (stage / OUTPUT_MARKER).write_text(
                json.dumps({"kind": "rrmg_figure_data", "sha256": inventory(stage)}, indent=2) + "\n",
                encoding="utf-8",
            )
            check_output(path)
            backup = path.with_name("." + path.name + "-previous-" + uuid.uuid4().hex)
            had_old = path.exists()
            if had_old:
                path.rename(backup)
            try:
                stage.rename(path)
            except OSError:
                if had_old:
                    backup.rename(path)
                raise
            if had_old:
                shutil.rmtree(backup)
    finally:
        lock.rmdir()


def read_tsv(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Missing/empty table: " + str(path))
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def exact_summary(path, expected_ids):
    rows = read_tsv(path)
    if not rows or "RunID" not in rows[0]:
        raise ValueError("Summary lacks RunID: " + str(path))
    ids = [row["RunID"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate RunID in summary: " + str(path))
    missing = sorted(set(expected_ids) - set(ids))
    extra = sorted(set(ids) - set(expected_ids))
    if missing or extra or len(ids) != len(expected_ids):
        raise ValueError(f"Incomplete summary {path}: missing={missing}; extra={extra}")
    order = {value: index for index, value in enumerate(expected_ids)}
    return sorted(rows, key=lambda row: order[row["RunID"]])


def build_manifest(root):
    rows = []
    for path in sorted(Path(root).rglob("*")):
        if path.is_file() and path.name not in {OUTPUT_MARKER, "FILE_MANIFEST.tsv", "COMPLETE.json"}:
            rows.append({
                "relative_path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": digest(path),
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wd", required=True, type=Path)
    parser.add_argument("--plan", default="input/rrmg/rrmg_permutation_plan.tsv")
    parser.add_argument("--root", default="temp/rrmg_permutation")
    parser.add_argument("--gene", required=True, type=input_file)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--sites", type=input_file)
    parser.add_argument("--catalogue", type=input_file)
    args = parser.parse_args()
    if bool(args.sites) != bool(args.catalogue):
        parser.error("Use --sites and --catalogue together")

    wd = args.wd.resolve()
    plan_path = Path(args.plan)
    if not plan_path.is_absolute():
        plan_path = wd / plan_path
    root = Path(args.root)
    if not root.is_absolute():
        root = wd / root
    out = args.out.resolve() if args.out else wd / "figure_data"
    plan = read_tsv(plan_path)
    expected_ids = [row.get("run_id", "") for row in plan]
    indices = [row.get("run_index", "") for row in plan]
    if not plan or any(not value for value in expected_ids + indices):
        parser.error("Plan requires non-empty run_id and run_index")
    if len(expected_ids) != len(set(expected_ids)) or len(indices) != len(set(indices)):
        parser.error("Duplicate plan run_id/run_index")

    summaries = {
        "contig_mapping/rrmg_target_gene_map_merged.tsv": wd / "result/rrmg_target_gene_map_merged.tsv",
        "read_to_gene/rrmg_reads_to_gene_depth_summary.tsv": wd / "result/rrmg_reads_to_gene_depth_summary.tsv",
        "read_to_gene/rrmg_target_only_reads_to_gene_depth_summary.tsv": wd / "result/rrmg_target_only_reads_to_gene_depth_summary.tsv",
        "target_genome/rrmg_afum_target_reads_genome_coverage_summary.tsv": wd / "result/rrmg_afum_target_reads_genome_coverage_summary.tsv",
    }
    checked = {relative: exact_summary(path, expected_ids) for relative, path in summaries.items()}

    try:
        check_output(out)
        with managed_output(out) as stage:
            (stage / "metadata").mkdir()
            shutil.copy2(plan_path, stage / "metadata/selected_plan.tsv")
            version_file = Path(__file__).resolve().with_name("VERSION")
            if version_file.is_file():
                shutil.copy2(version_file, stage / "metadata/pipeline_VERSION")
            for relative, source in summaries.items():
                destination = stage / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                write_tsv(destination, checked[relative], list(checked[relative][0]))

            alignment_out = stage / "CYP51A_alignment"
            run_count, site_count = build_alignment_bundle(
                wd, args.gene, alignment_out, "all", args.sites, args.catalogue, plan_path
            )

            for row in plan:
                plot_id = "row_" + str(int(row["run_index"])).zfill(5)
                source = root / row["run_id"] / "target_readmap/afum_target_reads_to_genome.bins.tsv"
                if not source.is_file() or source.stat().st_size == 0:
                    raise ValueError("Missing/empty target-genome bins: " + str(source))
                destination = stage / "target_genome" / plot_id / source.name
                destination.parent.mkdir(parents=True)
                shutil.copy2(source, destination)

            files = build_manifest(stage)
            write_tsv(stage / "FILE_MANIFEST.tsv", files, ["relative_path", "bytes", "sha256"])
            manifest_sha256 = digest(stage / "FILE_MANIFEST.tsv")
            complete = {
                "kind": "rrmg_figure_data",
                "sample_count": len(plan),
                "alignment_sample_count": run_count,
                "validated_site_count": site_count,
                "source_workdir": str(wd),
                "plan": str(plan_path),
                "gene": str(Path(args.gene).resolve()),
                "file_manifest_sha256": manifest_sha256,
            }
            (stage / "COMPLETE.json").write_text(json.dumps(complete, indent=2) + "\n", encoding="utf-8")
    except (ValueError, OSError) as exc:
        parser.error(str(exc))

    print(f"[DONE] complete portable figure data: {out}")
    print(f"[DONE] samples={len(plan)} files={len(files)}")


if __name__ == "__main__":
    main()
