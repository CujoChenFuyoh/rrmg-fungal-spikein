"""Portable figure-data collector tests using two tiny synthetic runs."""

import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "rrmg_collect_figure_data.py"
SHELL_COLLECTOR = REPO / "rrmg_collect_figure_data.sh"


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class FigureDataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.wd = Path(self.temp.name)
        (self.wd / "ref").mkdir()
        self.gene = self.wd / "ref/gene.fasta"
        self.gene.write_text(">toy\nATGAAACCC\n")
        plan_rows = [
            {"run_index": str(index), "run_id": f"run{index}", "seed": str(index),
             "target_abundance": "0.01", "total_pairs": "1000"}
            for index in (1, 2)
        ]
        write_tsv(self.wd / "input/rrmg/rrmg_permutation_plan.tsv", plan_rows)
        summary_rows = [{"RunID": f"run{index}", "RowIndex": str(index)} for index in (1, 2)]
        write_tsv(self.wd / "result/rrmg_target_gene_map_merged.tsv", summary_rows)
        write_tsv(self.wd / "result/rrmg_reads_to_gene_depth_summary.tsv", summary_rows)
        write_tsv(self.wd / "result/rrmg_target_only_reads_to_gene_depth_summary.tsv", summary_rows)
        write_tsv(self.wd / "result/rrmg_afum_target_reads_genome_coverage_summary.tsv", summary_rows)
        paf = "q\t9\t0\t9\t+\ttoy\t9\t0\t9\t9\t9\t60\n"
        depth = "".join(f"toy\t{position}\t1\n" for position in range(1, 10))
        for index in (1, 2):
            run = self.wd / f"temp/rrmg_permutation/run{index}"
            (run / "asm").mkdir(parents=True)
            (run / "readmap").mkdir()
            (run / "target_readmap").mkdir()
            (run / "map").mkdir()
            (run / "asm/contigs.fasta.to_GENE.paf").write_text(paf)
            for prefix in ("reads_to_GENE", "target_only_reads_to_GENE"):
                (run / f"readmap/{prefix}.depth.tsv").write_text(depth)
                with gzip.open(run / f"readmap/{prefix}.paf.gz", "wt") as handle:
                    handle.write(paf)
                write_tsv(
                    run / f"readmap/{prefix}.coverage_summary.tsv",
                    [{"RunID": f"run{index}", "RowIndex": str(index), "MeanDepth": "1.0"}],
                )
            write_tsv(
                run / "map/target_call_summary.tsv",
                [{"SampleID": f"sample{index}", "TargetType": "GENE", "TargetFA": str(self.gene),
                  "ContigsFA": str(run / "asm/contigs.fasta"), "HasHit": "Y", "BestAlnLen": "9",
                  "BestIdentity": "100", "TargetBreadthPct": "100", "TargetCoveredBp": "9",
                  "VCF": "", "VarN": "0"}],
            )
            write_tsv(
                run / "target_readmap/afum_target_reads_to_genome.summary.tsv",
                [{"RunID": f"run{index}", "RowIndex": str(index), "ObservedMeanDepth": "1.0"}],
            )
            (run / "target_readmap/afum_target_reads_to_genome.bins.tsv").write_text(
                "contig\tstart0\tend\tbin_len\tdepth_sum\tmean_depth\nchr1\t0\t9\t9\t9\t1\n"
            )
        write_tsv(
            self.wd / "result/rrmg_split_pipeline_summary.tsv",
            [{"RunDir": str(self.wd / f"temp/rrmg_permutation/run{index}"),
              "Seed": str(index), "Sample": f"sample{index}",
              "TargetAbundance": "0.01", "TotalPairs": "1000"}
             for index in (1, 2)],
        )

    def run_collector(self, ok=True):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--wd", str(self.wd), "--gene", str(self.gene)],
            capture_output=True, text=True,
        )
        if ok:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def test_complete_bundle_and_managed_rerun(self):
        self.run_collector()
        out = self.wd / "figure_data"
        complete = json.loads((out / "COMPLETE.json").read_text())
        self.assertEqual(complete["sample_count"], 2)
        self.assertEqual(complete["alignment_sample_count"], 2)
        self.assertTrue((out / "metadata/selected_plan.tsv").is_file())
        self.assertTrue((out / "CYP51A_alignment/row_00002/target_only_reads_to_GENE.paf.gz").is_file())
        self.assertTrue((out / "target_genome/row_00001/afum_target_reads_to_genome.bins.tsv").is_file())
        self.run_collector()
        (out / "manual.txt").write_text("preserve\n")
        self.run_collector(ok=False)
        self.assertEqual((out / "manual.txt").read_text(), "preserve\n")

    def test_missing_sample_fails_without_partial_output(self):
        (self.wd / "temp/rrmg_permutation/run2/readmap/target_only_reads_to_GENE.depth.tsv").unlink()
        self.run_collector(ok=False)
        self.assertFalse((self.wd / "figure_data").exists())

    @unittest.skipUnless(os.name != "nt", "Bash collector orchestration is tested on Linux")
    def test_shell_collector_rebuilds_all_four_summaries_then_bundles(self):
        for name in (
            "rrmg_target_gene_map_merged.tsv",
            "rrmg_reads_to_gene_depth_summary.tsv",
            "rrmg_target_only_reads_to_gene_depth_summary.tsv",
            "rrmg_afum_target_reads_genome_coverage_summary.tsv",
        ):
            (self.wd / "result" / name).unlink()
        result = subprocess.run(
            ["bash", str(SHELL_COLLECTOR), "--wd", str(self.wd), "--gene", str(self.gene)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.wd / "figure_data/COMPLETE.json").is_file())
        target_summary = self.wd / "figure_data/read_to_gene/rrmg_target_only_reads_to_gene_depth_summary.tsv"
        with target_summary.open(newline="") as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 2)


if __name__ == "__main__":
    unittest.main()
