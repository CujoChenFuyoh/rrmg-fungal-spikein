"""Reference-coordinate and input-safety tests; all generated data are toy fixtures."""
import importlib.util
import argparse
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / "rrmg_prepare_alignment_inputs.py"
spec = importlib.util.spec_from_file_location("alignment", MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class AlignmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gene = self.root / "gene.fasta"
        self.gene.write_text(">toy NC_000001.1:c108-100\nATGAAACCC\n")
        self.catalogue = self.root / "catalogue.tsv"
        self.sites = self.root / "sites.tsv"
        self.cat = dict(marker_id="test_M1V", gene="CYP51A", aa_change="M1V", chrom="NC_000001.1",
                        strand="-", codon_positions="108|107|106", ref_codon="ATG", status="ok")
        self.site = dict(label="M1V", plot_start=1, plot_end=3, plot_mid=2, is_key_site="Y")
        self.write()

    def write(self):
        mod.write_table(self.catalogue, [self.cat], list(self.cat))
        mod.write_table(self.sites, [self.site], list(self.site))

    def test_reverse_reference_exact_codon(self):
        sites, audit = mod.validated_sites(self.gene, self.sites, self.catalogue)
        self.assertEqual((sites[0]["plot_start"], sites[0]["plot_end"]), (1, 3))
        self.assertEqual(sites[0]["ref_codon"], "ATG")
        self.assertEqual(audit[0]["status"], "included")

    def test_forward_reference(self):
        self.gene.write_text(">toy NC_000001.1:100-108\nATGAAACCC\n")
        self.cat.update(strand="+", codon_positions="100|101|102")
        self.write()
        self.assertEqual(mod.validated_sites(self.gene, self.sites, self.catalogue)[0][0]["plot_mid"], 2)

    def test_wrong_ref_or_coordinates_rejected(self):
        for field, value in (("status", "reference_aa_mismatch"), ("ref_codon", "AAA"),
                             ("codon_positions", "10810710600000"), ("chrom", "NC_000002.1")):
            original = self.cat[field]
            self.cat[field] = value
            self.write()
            with self.assertRaises(ValueError):
                mod.validated_sites(self.gene, self.sites, self.catalogue)
            self.cat[field] = original
        self.site["plot_start"] = 2
        self.write()
        with self.assertRaises(ValueError):
            mod.validated_sites(self.gene, self.sites, self.catalogue)

    def test_depth_and_gzip_paf_validation(self):
        depth = self.root / "depth.tsv"
        depth.write_text("".join(f"toy\t{i}\t{i % 2}\n" for i in range(1, 10)))
        self.assertAlmostEqual(mod.depth_stats(depth, "toy", 9), 100 * 5 / 9)
        with gzip.open(self.root / "reads.paf.gz", "wt") as handle:
            handle.write("q\t9\t0\t9\t+\ttoy\t9\t0\t9\t9\t9\t60\n")
        self.assertEqual(mod.paf_stats(self.root / "reads.paf.gz", "toy", 9), (1, 100.0))

    def test_empty_paf_is_valid_zero_hit(self):
        paf = self.root / "test.paf"
        paf.touch()
        self.assertEqual(mod.paf_stats(paf, "toy", 9), (0, 0))
        paf.write_text("malformed\n")
        with self.assertRaises(ValueError):
            mod.paf_stats(paf, "toy", 9)

    def test_empty_directory_and_missing_file_arguments(self):
        for value in ("", " ", str(self.root), str(self.root / "missing.tsv")):
            with self.assertRaises(argparse.ArgumentTypeError):
                mod.input_file(value)
        self.assertEqual(mod.input_file(str(self.catalogue)), self.catalogue)
        result = subprocess.run([sys.executable, str(MODULE), "--wd", str(self.root),
                                 "--gene", str(self.gene), "--out", str(self.root / "out"),
                                 "--sites", str(self.sites), "--catalogue", ""],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--catalogue: Empty file path", result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        self.assertFalse((self.root / "out").exists())

    def test_collects_existing_step3_step4_outputs_without_mapping(self):
        plan = self.root / "input/rrmg/rrmg_permutation_plan.tsv"
        plan.parent.mkdir(parents=True)
        plan.write_text(
            "run_index\trun_id\tseed\ttarget_abundance\ttotal_pairs\n"
            "1\trun1\t1\t0.01\t1000\n"
        )
        run = self.root / "temp/rrmg_permutation/run1"
        (run / "asm").mkdir(parents=True)
        (run / "readmap").mkdir()
        (run / "asm/contigs.fasta.to_GENE.paf").touch()
        depth_text = "".join(f"toy\t{i}\t1\n" for i in range(1, 10))
        for prefix in ("reads_to_GENE", "target_only_reads_to_GENE"):
            (run / f"readmap/{prefix}.depth.tsv").write_text(depth_text)
            with gzip.open(run / f"readmap/{prefix}.paf.gz", "wt") as handle:
                handle.write("q\t9\t0\t9\t+\ttoy\t9\t0\t9\t9\t9\t60\n")
        out = self.root / "collected"
        count, sites = mod.build_alignment_bundle(self.root, self.gene, out)
        self.assertEqual((count, sites), (1, 0))
        self.assertTrue((out / "row_00001/reads_to_GENE.paf.gz").is_file())
        self.assertTrue((out / "row_00001/target_only_reads_to_GENE.depth.tsv").is_file())
        self.assertTrue((out / "COMPLETE.json").is_file())

    def test_fixed_output_rerun_and_failure(self):
        out = self.root / "fixed"
        with mod.managed_output(out) as stage:
            (stage / "old.tsv").write_text("previous\n")
        with self.assertRaises(RuntimeError):
            with mod.managed_output(out) as stage:
                (stage / "partial.tsv").write_text("partial\n")
                raise RuntimeError("mapping failed")
        self.assertEqual((out / "old.tsv").read_text(), "previous\n")
        self.assertFalse((out / "partial.tsv").exists())
        with mod.managed_output(out) as stage:
            (stage / "new.tsv").write_text("complete\n")
        self.assertFalse((out / "old.tsv").exists())
        self.assertEqual((out / "new.tsv").read_text(), "complete\n")
        self.assertFalse((self.root / "fixed.lock").exists())

    def test_protect_unmanaged_modified_added_and_locked_outputs(self):
        out = self.root / "fixed"
        out.mkdir()
        (out / "notes.txt").write_text("user\n")
        with self.assertRaises(ValueError), mod.managed_output(out):
            pass
        self.assertEqual((out / "notes.txt").read_text(), "user\n")
        managed = self.root / "managed"
        with mod.managed_output(managed) as stage:
            (stage / "result.tsv").write_text("tool\n")
        (managed / "extra.txt").write_text("user\n")
        with self.assertRaises(ValueError), mod.managed_output(managed):
            pass
        (managed / "extra.txt").unlink()
        (managed / "result.tsv").write_text("edited\n")
        with self.assertRaises(ValueError), mod.managed_output(managed):
            pass
        self.assertEqual((managed / "result.tsv").read_text(), "edited\n")
        (self.root / "locked.lock").mkdir()
        with self.assertRaisesRegex(ValueError, "Output is locked"), mod.managed_output(self.root / "locked"):
            pass


if __name__ == "__main__":
    unittest.main()
