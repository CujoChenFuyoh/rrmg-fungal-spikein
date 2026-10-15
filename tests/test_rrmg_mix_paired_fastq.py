import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zlib


REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("rrmg_mix_paired_fastq", REPO / "rrmg_mix_paired_fastq.py")
MIXER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MIXER)


def record(name, mate):
    sequence = ("ACGT" if mate == 1 else "TGCA") * 3
    return f"@{name}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode()


def append_gzip_member(path, names, mate):
    with path.open("ab") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
            for name in names:
                out.write(record(name, mate))


def read_names(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [lines[i][1:].split("/")[0] for i in range(0, len(lines), 4)]


def gzip_members(path):
    remaining = path.read_bytes()
    count = 0
    while remaining:
        stream = zlib.decompressobj(16 + zlib.MAX_WBITS)
        stream.decompress(remaining)
        stream.flush()
        count += 1
        if not stream.unused_data:
            break
        remaining = stream.unused_data
    return count


class ExactMixTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="rrmg-mix-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bg1 = self.root / "background_1.fastq.gz"
        self.bg2 = self.root / "background_2.fastq.gz"
        self.target1 = self.root / "target_1.fastq.gz"
        self.target2 = self.root / "target_2.fastq.gz"
        # Deliberately concatenate two gzip members to reproduce BBMap BGZF-like input.
        for path, mate in ((self.bg1, 1), (self.bg2, 2)):
            append_gzip_member(path, ["bg0", "bg1"], mate)
            append_gzip_member(path, ["bg2", "bg3", "bg4"], mate)
        for path, mate in ((self.target1, 1), (self.target2, 2)):
            append_gzip_member(path, ["target0", "target1", "target2"], mate)

    def mix(self, prefix):
        out1 = self.root / f"{prefix}_1.fastq.gz"
        out2 = self.root / f"{prefix}_2.fastq.gz"
        audit = self.root / f"{prefix}.audit.tsv"
        result = MIXER.main([
            "mix", "--background-r1", str(self.bg1), "--background-r2", str(self.bg2),
            "--background-source-pairs", "5", "--background-pairs", "4",
            "--target-r1", str(self.target1), "--target-r2", str(self.target2),
            "--target-source-pairs", "3", "--target-pairs", "3",
            "--out-r1", str(out1), "--out-r2", str(out2), "--audit", str(audit),
            "--seed", "99", "--compressor", "python",
        ])
        self.assertEqual(result, 0)
        return out1, out2, audit

    def test_exact_pairing_random_interleave_and_single_member(self):
        self.assertEqual(MIXER.count_pairs(self.bg1, self.bg2), 5)
        out1, out2, audit = self.mix("first")
        names1, names2 = read_names(out1), read_names(out2)
        self.assertEqual(names1, names2)
        self.assertEqual(len(names1), 7)
        self.assertEqual(sum(name.startswith("bg") for name in names1), 4)
        self.assertEqual(sum(name.startswith("target") for name in names1), 3)
        self.assertNotEqual(names1[-3:], ["target0", "target1", "target2"])
        self.assertEqual(gzip_members(out1), 1)
        self.assertEqual(gzip_members(out2), 1)
        with audit.open(newline="") as handle:
            values = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(values["total_pairs_written"], "7")
        self.assertEqual(values["gzip_members_per_output"], "1")
        self.assertEqual(len(values["mix_generation_id"]), 32)
        self.assertGreater(int(values["source_switches"]), 0)

        again1, again2, again_audit = self.mix("again")
        self.assertEqual(gzip.open(out1, "rb").read(), gzip.open(again1, "rb").read())
        self.assertEqual(gzip.open(out2, "rb").read(), gzip.open(again2, "rb").read())
        with again_audit.open(newline="") as handle:
            again_values = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertNotEqual(values["mix_generation_id"], again_values["mix_generation_id"])

    def test_declared_source_count_mismatch_is_rejected(self):
        with self.assertRaises(SystemExit):
            MIXER.main([
                "mix", "--background-r1", str(self.bg1), "--background-r2", str(self.bg2),
                "--background-source-pairs", "4", "--background-pairs", "4",
                "--target-pairs", "0", "--out-r1", str(self.root / "bad1.gz"),
                "--out-r2", str(self.root / "bad2.gz"), "--audit", str(self.root / "bad.tsv"),
                "--seed", "1", "--compressor", "python",
            ])
        self.assertFalse((self.root / "bad1.gz").exists())
        self.assertFalse((self.root / "bad2.gz").exists())

    def test_fastp_and_mix_audits_fail_closed(self):
        out1, out2, audit = self.mix("audit")
        MIXER.validate_mix_audit(audit, 7)
        MIXER.validate_mix_audit(audit, 7, out1, out2, 4, 3)
        with self.assertRaises(MIXER.FastqError):
            MIXER.validate_mix_audit(audit, 8)
        with self.assertRaises(MIXER.FastqError):
            MIXER.validate_mix_audit(audit, 7, expected_target_pairs=2)
        with out1.open("ab") as handle:
            handle.write(b"unexpected")
        with self.assertRaises(MIXER.FastqError):
            MIXER.validate_mix_audit(audit, 7, out1, out2)
        report = self.root / "fastp.json"
        report.write_text(json.dumps({"summary": {"before_filtering": {"total_reads": 14}}}))
        MIXER.validate_fastp(report, 7)
        with self.assertRaises(MIXER.FastqError):
            MIXER.validate_fastp(report, 8)

    def test_mate_name_mismatch_is_rejected(self):
        wrong = self.root / "wrong_2.fastq.gz"
        append_gzip_member(wrong, ["different0", "bg1", "bg2", "bg3", "bg4"], 2)
        with self.assertRaises(MIXER.FastqError):
            MIXER.count_pairs(self.bg1, wrong)


if __name__ == "__main__":
    unittest.main()
