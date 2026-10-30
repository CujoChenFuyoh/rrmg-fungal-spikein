"""Toy sequence evidence tests; no project data or biological performance claims."""
import csv
import gzip
import importlib.util
import json
import os
import random
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import rrmg_gene_evidence as evidence
import rrmg_pileup_base_counts as counts
import rrmg_prepare_alignment_inputs as collector


class BaseEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.gene = self.root / "gene.fa"
        self.gene.write_text(">toy\n" + "A" * 200 + "\n")

    def test_base_and_strand_counts(self):
        r = counts.parse_read_bases("A", ".,TtNn", 6)
        self.assertEqual((r["A"], r["T"], r["N"]), (2, 2, 2))
        self.assertEqual((r["T_F"], r["T_R"]), (1, 1))
        self.assertEqual(r["ObservationStatus"], "NONREF_OBSERVED")

    def test_indel_start_end_and_deletion_placeholders(self):
        r = counts.parse_read_bases("C", "^].+2tt,-1g$*#<>", 6)
        self.assertEqual(r["InsertionAlleles"], "TT:1")
        self.assertEqual(r["DeletionAlleles"], "G:1")
        self.assertEqual((r["DeletionPlaceholders"], r["ReferenceSkips"]), (2, 2))
        self.assertEqual(r["SNVNonRefCount"], 0)

    def test_malformed_and_truncated_pileup_rejected(self):
        for text, depth in (("^", 1), (".+AC", 1), (".+3A", 1), (".", 2), ("?", 1)):
            with self.subTest(text=text), self.assertRaises(ValueError):
                counts.parse_read_bases("A", text, depth)

    def test_zero_design_gzip_and_resume_fingerprint(self):
        prefix = self.root / "target_only_reads_to_GENE"
        self.assertEqual(evidence.ensure_reads(self.gene, prefix, "synthetic_target", "gen1", empty=True), "OK/OK")
        output = Path(str(prefix) + ".base_counts.q20.tsv.gz")
        with gzip.open(output, "rt") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(len(rows), 200)
        self.assertTrue(all(r["ObservationStatus"] == "NO_COVERAGE" for r in rows))
        before = output.stat().st_mtime_ns
        self.assertEqual(evidence.ensure_reads(self.gene, prefix, "synthetic_target", "gen1", empty=True), "SKIP_CURRENT/SKIP_CURRENT")
        self.assertEqual(output.stat().st_mtime_ns, before)
        self.assertEqual(evidence.ensure_reads(self.gene, prefix, "synthetic_target", "gen2", empty=True), "OK/OK")

    def test_missing_bam_not_zero(self):
        with self.assertRaisesRegex(ValueError, "cannot infer zero"):
            evidence.ensure_reads(self.gene, self.root / "missing", "metagenomic")

    def test_reference_whitespace_does_not_change_sha(self):
        a = counts.reference_inventory(self.gene)[1]
        self.gene.write_bytes(b">toy\r\n" + b"A" * 100 + b"\r\n" + b"A" * 100 + b"\r\n")
        self.assertEqual(counts.reference_inventory(self.gene)[1], a)

    def paf(self, lines):
        path = self.root / "contigs.cs.paf"
        path.write_text("".join(lines))
        return path

    def record(self, name, start, end, cs=None, strand="+"):
        n = end - start
        return f"{name}\t{n}\t0\t{n}\t{strand}\ttoy\t200\t{start}\t{end}\t{n}\t{n}\t60\tcg:Z:{n}M\tcs:Z:{cs or ':'+str(n)}\n"

    def test_full_is_per_contig_not_combined_breadth(self):
        # 200-bp toy reference; records must pass the 150-bp recovery filter.
        p = self.paf([self.record("c1", 0, 160), self.record("c2", 40, 200)])
        r, _ = evidence.contig_details(p, "A" * 200)
        self.assertEqual(r["RecoveryClass"], "FRAGMENTED_COMPLETE")
        self.assertEqual(float(r["MaxSingleContigBreadth"]), 80)
        r, _ = evidence.contig_details(self.paf([self.record("c1", 0, 190)]), "A" * 200)
        self.assertEqual(r["RecoveryClass"], "FULL_SINGLE")

    def test_reverse_strand_snv_is_reference_oriented(self):
        r, changes = evidence.contig_details(self.paf([self.record("c1", 0, 200, ":99*at:100", "-")]), "A" * 200)
        self.assertEqual((changes[0]["GenePosition"], changes[0]["RefAllele"], changes[0]["QueryAllele"]), (100, "A", "T"))
        self.assertEqual(r["DifferenceEventN"], 1)

    def test_indel_cs_coordinates_and_bad_reference(self):
        p = self.paf(["c\t201\t0\t201\t+\ttoy\t200\t0\t200\t199\t202\t60\tcg:Z:99M2I1D100M\tcs:Z::99+tt-a:100\n"])
        _, events = evidence.contig_details(p, "A" * 200)
        self.assertEqual([(e["EventType"], e["GenePosition"]) for e in events], [("INS", 99), ("DEL", 100)])
        with self.assertRaises(ValueError):
            evidence.contig_details(self.paf([self.record("c", 0, 200, ":99*ct:100")]), "A" * 200)

    def test_plan_holes_and_uniqueness(self):
        plan = self.root / "plan.tsv"
        plan.write_text("run_index\trun_id\n1\tr1\n300\tr300\n")
        self.assertEqual(len(evidence.select_plan(plan, "all")), 2)
        with self.assertRaises(ValueError): evidence.select_plan(plan, "1-300")
        plan.write_text("run_index\trun_id\n1\tr1\n2\tr1\n")
        with self.assertRaises(ValueError): evidence.select_plan(plan, "all")

    def test_generation_unknown_not_match(self):
        self.assertEqual(evidence.generation_status("gen1", ""), "UNVERIFIED")
        self.assertEqual(evidence.generation_status("gen1", "gen2"), "MISMATCH")
        self.assertEqual(evidence.generation_status("gen1", "gen1"), "MATCH")

    def test_collector_includes_validated_evidence_and_recovery(self):
        plan = self.root / "input/rrmg/rrmg_permutation_plan.tsv"
        plan.parent.mkdir(parents=True)
        plan.write_text("run_index\trun_id\tseed\ttarget_abundance\ttotal_pairs\n1\tr1\t1\t0\t10\n")
        run = self.root / "temp/rrmg_permutation/r1"
        (run / "asm").mkdir(parents=True)
        (run / "readmap").mkdir()
        (run / "design").mkdir()
        (run / "design/read_mix_audit.tsv").write_text("key\tvalue\nmix_generation_id\tgen1\n")
        paf = run / "asm/contigs.fasta.to_GENE.paf"
        paf.write_text(self.record("c1", 0, 200))
        for prefix in ("reads_to_GENE", "target_only_reads_to_GENE"):
            (run / f"readmap/{prefix}.depth.tsv").write_text("".join(f"toy\t{i}\t0\n" for i in range(1,201)))
            with gzip.open(run / f"readmap/{prefix}.paf.gz", "wt") as stream: stream.write("")
        out = self.root / "strict-missing"
        with self.assertRaisesRegex(ValueError, "run rrmg_backfill"):
            collector.build_alignment_bundle(self.root, self.gene, out, require_base_evidence=True)
        cs = run / "asm/contigs.fasta.to_GENE.cs.paf"
        cs.write_text(self.record("c1", 0, 200))
        refsha = counts.reference_inventory(self.gene)[1]
        evidence.mark(cs, dict(schema=1, reference_sha256=refsha, reference_name="toy", generation="gen1"))
        for prefix, label in (("reads_to_GENE", "metagenomic"), ("target_only_reads_to_GENE", "synthetic_target")):
            evidence.ensure_reads(self.gene, run / f"readmap/{prefix}", label, "gen1", empty=True)
        out = self.root / "strict-complete"
        collector.build_alignment_bundle(self.root, self.gene, out, require_base_evidence=True)
        self.assertTrue(json.loads((out / "COMPLETE.json").read_text())["base_evidence_complete"])
        with (out / "contig_recovery.tsv").open() as stream:
            row = next(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(row["RecoveryClass"], "FULL_SINGLE")
        self.assertEqual(row["GenerationStatus"], "MATCH")
        # Corruption must not silently become a zero-evidence sample.
        cs.write_text("bad\n")
        with self.assertRaisesRegex(ValueError, "Missing/invalid"):
            collector.build_alignment_bundle(self.root, self.gene, self.root / "bad", require_base_evidence=True)


@unittest.skipUnless(sys.platform.startswith("linux") and shutil.which("minimap2") and shutil.which("samtools"),
                     "Optional real minimap2/samtools evidence integration")
class RealEvidenceTests(unittest.TestCase):
    def test_twelve_runs_four_toy_groups_backfill_collect_and_preserve(self):
        with tempfile.TemporaryDirectory(prefix="rrmg-evidence-real-") as tmp:
            wd = Path(tmp)
            (wd / "ref").mkdir()
            rng = random.Random(934)
            seq = "".join(rng.choice("ACGT") for _ in range(1200))
            gene = wd / "ref/CYP51A_wt.fasta"
            gene.write_text(">toy\n" + seq + "\n")
            subprocess.run(["samtools", "faidx", str(gene)], check=True)
            alt = next(b for b in "ACGT" if b != seq[300])
            altered = seq[:300] + alt + seq[301:]
            unknown = "".join(rng.choice("ACGT") for _ in range(1200))
            plans = []
            protected = {}
            for i in range(1, 13):
                group = (i - 1) // 3
                rid = f"group{group}_seed{i}"
                plans.append(dict(run_index=i, run_id=rid, seed=i,
                                  target_abundance=0.1 + group * 0.1, total_pairs=10))
                run = wd / "temp/rrmg_permutation" / rid
                for folder in ("asm", "readmap", "design"):
                    (run / folder).mkdir(parents=True)
                contigs = ([('whole', seq)] if group == 0 else [('whole_alt', altered)] if group == 1
                           else [('left', seq[:650]), ('right', seq[550:])] if group == 2 else [('unrelated', unknown)])
                (run / "asm/contigs.fasta").write_text("".join(f">{name}\n{s}\n" for name, s in contigs))
                (run / "asm/RRMG_MIX_GENERATION").write_text("gen1\n")
                (run / "design/read_mix_audit.tsv").write_text("key\tvalue\nmix_generation_id\tgen1\n")
                # Legacy files are never rewritten by evidence backfill.
                paf = run / "asm/contigs.fasta.to_GENE.paf"
                paf.write_text(f"legacy\t1200\t0\t1200\t+\ttoy\t1200\t0\t1200\t1200\t1200\t60\n")
                for prefix, synthetic in (("reads_to_GENE", False), ("target_only_reads_to_GENE", True)):
                    records = [("ref", 0, 60, seq[250:400], "I" * 150)]
                    if not synthetic:
                        records += [("alt", 16, 60, altered[250:400], "I" * 150),
                                    ("lowmq", 0, 5, altered[250:400], "I" * 150),
                                    ("lowbq", 0, 60, altered[250:400], "!" * 150)]
                    sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:toy\tLN:1200\n" + "".join(
                        f"{name}\t{flag}\ttoy\t251\t{mq}\t150M\t*\t0\t0\t{s}\t{qual}\n"
                        for name, flag, mq, s, qual in records)
                    bam = run / f"readmap/{prefix}.bam"
                    subprocess.run(["samtools", "view", "-b", "-o", str(bam), "-"], input=sam, text=True, check=True)
                    (run / f"readmap/{prefix}.mix_generation.txt").write_text("gen1\n")
                    (run / f"readmap/{prefix}.depth.tsv").write_text("".join(f"toy\t{n}\t{1 if 251<=n<=400 else 0}\n" for n in range(1,1201)))
                    with gzip.open(run / f"readmap/{prefix}.paf.gz", "wt") as stream:
                        stream.write("q\t150\t0\t150\t+\ttoy\t1200\t250\t400\t150\t150\t60\n")
                for path in run.rglob("*"):
                    if path.is_file():
                        protected[path] = (evidence.sha256(path), path.stat().st_mtime_ns)
            plan = wd / "input/rrmg/rrmg_permutation_plan.tsv"
            plan.parent.mkdir(parents=True)
            with plan.open("w") as stream:
                writer = csv.DictWriter(stream, list(plans[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader(); writer.writerows(plans)
            command = ["bash", str(REPO / "rrmg_backfill_cyp51a_evidence.sh"), "-w", str(wd), "-r", "all", "-t", "2"]
            first = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
            out = wd / "alignment"
            collector.build_alignment_bundle(wd, gene, out, require_base_evidence=True)
            with (out / "contig_recovery.tsv").open() as stream:
                recoveries = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual([r["RecoveryClass"] for r in recoveries],
                             ["FULL_SINGLE"] * 6 + ["FRAGMENTED_COMPLETE"] * 3 + ["NOT_DETECTED"] * 3)
            self.assertTrue(all(r["GenerationStatus"] == "MATCH" for r in recoveries))
            with (out / "contig_differences.tsv").open() as stream:
                differences = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(differences), 3)
            self.assertTrue(all(r["GenePosition"] == "301" and r["QueryAllele"] == alt for r in differences))
            for level, expected_dp, expected_alt in (("raw", 4, 3), ("q20", 2, 1)):
                with gzip.open(out / f"row_00001/reads_to_GENE.base_counts.{level}.tsv.gz", "rt") as stream:
                    rows = list(csv.DictReader(stream, delimiter="\t"))
                self.assertEqual(len(rows), 1200)
                self.assertEqual(int(rows[300]["BaseCountDepth"]), expected_dp)
                self.assertEqual(int(rows[300][alt]), expected_alt)
                self.assertEqual(rows[0]["ObservationStatus"], "NO_COVERAGE")
            derived = {p: p.stat().st_mtime_ns for p in (wd / "temp").rglob("*.evidence.json")}
            second = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
            self.assertEqual(derived, {p: p.stat().st_mtime_ns for p in derived})
            self.assertEqual(protected, {p: (evidence.sha256(p), p.stat().st_mtime_ns) for p in protected})
            print("[OK] real evidence: 12 runs / 4 toy groups; exact SNV/base counts, resume and legacy preservation")


if __name__ == "__main__":
    unittest.main()
