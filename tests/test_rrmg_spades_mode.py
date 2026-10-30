"""Offline native Step3 mode/provenance tests; only tiny arbitrary DNA fixtures."""
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


@unittest.skipUnless(sys.platform.startswith('linux'), 'Linux shell utilities required')
class SpadesModeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='rrmg mode ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.run_dir = self.root / 'temp/rrmg_permutation/Toy__Seed1'
        reads = self.run_dir / 'reads'
        reads.mkdir(parents=True)
        for mate in (1, 2):
            with gzip.open(self.root / f'background_{mate}.fq.gz', 'wt') as handle:
                for i in range(2):
                    handle.write(f'@toy{i}/{mate}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n')
        self.invoke([
            sys.executable, REPO / 'rrmg_mix_paired_fastq.py', 'mix',
            '--background-r1', self.root / 'background_1.fq.gz',
            '--background-r2', self.root / 'background_2.fq.gz',
            '--background-source-pairs', 2, '--background-pairs', 2,
            '--out-r1', reads / 'Toy.Seed1_1.fastq.gz',
            '--out-r2', reads / 'Toy.Seed1_2.fastq.gz',
            '--audit', self.run_dir / 'design/read_mix_audit.tsv',
            '--seed', 1, '--compressor', 'python',
        ])
        self.plan = self.root / 'plan.tsv'
        self.plan.write_text('run_index\trun_id\tsample\tseed\ttotal_pairs\ttarget_abundance\n'
                             '1\tToy__Seed1\tToy\t1\t2\t0\n')
        self.gene = self.root / 'gene.fasta'
        self.gene.write_text('>toy\nACGTACGTACGT\n')
        bindir = self.root / 'bin'
        bindir.mkdir()
        tools = {
            'metaspades.py':
                'import json, os, pathlib, sys\n'
                'with open(os.environ["SPADES_CALLS"], "a") as f: f.write(json.dumps(sys.argv[1:]) + "\\n")\n'
                'if os.environ.get("SPADES_FAIL"): sys.exit(13)\n'
                'out = pathlib.Path(sys.argv[sys.argv.index("-o") + 1])\n'
                'out.mkdir(parents=True, exist_ok=True)\n'
                '(out / "contigs.fasta").write_text(">toy\\nACGTACGTACGT\\n")\n',
            'minimap2': 'import sys\nsys.exit(0)\n',
            'fastp':
                'import json, pathlib, shutil, sys\n'
                'args = sys.argv\n'
                'def value(flag): return args[args.index(flag) + 1]\n'
                'shutil.copyfile(value("-i"), value("-o"))\n'
                'shutil.copyfile(value("-I"), value("-O"))\n'
                'pathlib.Path(value("-j")).write_text(json.dumps({"summary":{"before_filtering":{"total_reads":4}}}))\n'
                'pathlib.Path(value("-h")).write_text("<html>fixture</html>")\n',
        }
        for name, content in tools.items():
            path = bindir / name
            path.write_text('#!' + sys.executable + '\n' + content)
            path.chmod(0o755)
        self.env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'],
                        SPADES_CALLS=str(self.root / 'spades_calls.jsonl'))
        self.env.pop('FASTQ_MIXER', None)
        self.env.pop('SPADES_FAIL', None)

    def invoke(self, command, env=None):
        return subprocess.run(list(map(str, command)), env=env, cwd=self.root,
                              capture_output=True, text=True, check=True)

    def stage3(self, *args, ok=True):
        command = ['bash', REPO / 'rrmg_step3_downstream_analysis.sh', self.plan, 1,
                   '-w', self.root, '-gene', self.gene, '--skip-kraken', '-t', 2, '-mem', 1, *args]
        result = subprocess.run(list(map(str, command)), env=self.env, cwd=self.root,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, ok, result.stdout + result.stderr)
        return result

    def calls(self):
        return [json.loads(line) for line in (self.root / 'spades_calls.jsonl').read_text().splitlines()]

    def test_default_skips_bayeshammer_but_still_runs_fastp(self):
        self.stage3()
        self.assertIn('--only-assembler', self.calls()[0])
        self.assertEqual((self.run_dir / 'asm/RRMG_ASSEMBLY_MODE').read_text().strip(), 'only-assembler')
        self.assertIn('assembly_mode=only-assembler\n', (self.run_dir / 'STEP3_SUCCESS').read_text())
        self.assertTrue((self.run_dir / 'qc/Toy.Seed1.fastp.json').is_file())
        # Same input/mode uses the completed row without launching SPAdes again.
        self.assertIn('[SKIP] Completed Step 3', self.stage3().stdout)
        self.assertEqual(len(self.calls()), 1)

    def test_completed_row_only_backfills_cs_without_touching_old_results(self):
        self.stage3()
        files = [self.run_dir / "asm/contigs.fasta", self.run_dir / "asm/contigs.fasta.to_GENE.paf",
                 self.run_dir / "map/target_call_summary.tsv", self.run_dir / "STEP3_SUCCESS",
                 self.run_dir / "qc/Toy.Seed1_1.fastq.gz"]
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in files}
        (self.run_dir / "asm/contigs.fasta.to_GENE.cs.paf").unlink()
        self.stage3()
        self.assertTrue((self.run_dir / "asm/contigs.fasta.to_GENE.cs.paf").is_file())
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in files})
        self.assertEqual(len(self.calls()), 1)

    def test_explicit_error_correction_is_recorded(self):
        self.stage3('--with-error-correction')
        self.assertNotIn('--only-assembler', self.calls()[0])
        self.assertEqual((self.run_dir / 'asm/RRMG_ASSEMBLY_MODE').read_text().strip(),
                         'error-correction-and-assembly')

    def test_mode_change_requires_explicit_rebuild(self):
        self.stage3()
        before = (self.run_dir / 'asm/contigs.fasta').read_bytes()
        result = self.stage3('--with-error-correction', ok=False)
        self.assertIn('Existing assembly mode is different or unrecorded', result.stderr)
        self.assertEqual((self.run_dir / 'asm/contigs.fasta').read_bytes(), before)
        self.assertEqual(len(self.calls()), 1)
        self.stage3('--with-error-correction', '--force')
        self.assertEqual(len(self.calls()), 2)
        self.assertNotIn('--only-assembler', self.calls()[1])
        self.assertIn('assembly_mode=error-correction-and-assembly\n',
                      (self.run_dir / 'STEP3_SUCCESS').read_text())

    def test_unrecorded_mode_is_not_silently_relabelled(self):
        self.stage3()
        marker = self.run_dir / 'asm/RRMG_ASSEMBLY_MODE'
        marker.unlink()
        self.assertIn('unrecorded', self.stage3(ok=False).stderr)
        self.assertFalse(marker.exists())
        self.assertEqual(len(self.calls()), 1)

    def test_skip_assemble_also_checks_mode(self):
        self.stage3()
        (self.run_dir / 'STEP3_SUCCESS').unlink()
        result = self.stage3('--with-error-correction', '--skip-assemble', ok=False)
        self.assertIn('Existing assembly mode is different or unrecorded', result.stderr)
        self.assertEqual(len(self.calls()), 1)

    def test_failure_never_marks_assembly_complete(self):
        self.env['SPADES_FAIL'] = '1'
        self.stage3(ok=False)
        self.assertFalse((self.run_dir / 'STEP3_SUCCESS').exists())
        self.assertFalse((self.run_dir / 'asm/RRMG_ASSEMBLY_MODE').exists())

    def test_fasta_line_endings_do_not_change_contig_breadth(self):
        # Nine covered bases in a 12-base wrapped reference must be 75%,
        # regardless of whether it was saved on Linux or Windows.
        minimap = self.root / 'bin/minimap2'
        minimap.write_text('#!' + sys.executable + '\n'
                           'print("contig\\t12\\t0\\t9\\t+\\ttoy\\t12\\t0\\t9\\t9\\t9\\t60\\tcg:Z:9M\\tcs:Z::9")\n')
        self.gene.write_bytes(b'>toy\nACGTAC\nGTACGT\n')
        self.stage3()
        with (self.run_dir / 'map/target_call_summary.tsv').open() as f:
            linux = next(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(linux['HasHit'], 'Y')
        self.assertEqual(linux['TargetCoveredBp'], '9')
        self.assertEqual(linux['TargetBreadthPct'], '75.000000')
        protected = [self.run_dir / 'asm/contigs.fasta',
                     self.run_dir / 'qc/Toy.Seed1.fastp.json',
                     self.run_dir / 'qc/Toy.Seed1_1.fastq.gz',
                     self.run_dir / 'qc/Toy.Seed1_2.fastq.gz']
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in protected}
        self.gene.write_bytes(b'>toy\r\nACGTAC\r\nGTACGT\r\n')
        self.stage3('--force', '--skip-qc', '--skip-assemble')
        with (self.run_dir / 'map/target_call_summary.tsv').open() as f:
            windows = next(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(windows['TargetBreadthPct'], linux['TargetBreadthPct'])
        self.assertEqual(windows['TargetCoveredBp'], linux['TargetCoveredBp'])
        self.assertEqual(len(self.calls()), 1)
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in protected})


if __name__ == '__main__':
    unittest.main()
