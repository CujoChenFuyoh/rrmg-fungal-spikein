"""Offline input-preparation and active-environment tests; no network/PBS."""
import csv
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]


@unittest.skipUnless(sys.platform.startswith("linux"), "Linux shell utilities required")
class PortabilityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="rrmg portable ")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.env = dict(os.environ)
        for key in ("JAVA_HOME", "JAVA_BIN", "BBMAP_DIR", "RANDOMREADSMG_SH", "RANDOMREADS_SH", "RRMG_CONFIG"):
            self.env.pop(key, None)

    def run_tool(self, name, *args, ok=True):
        command = [sys.executable if name.endswith('.py') else 'bash', str(REPO / name), *map(str, args)]
        r = subprocess.run(command, cwd=self.root, env=self.env, capture_output=True, text=True)
        self.assertEqual(r.returncode == 0, ok, r.stdout + r.stderr)
        return r

    def table(self, path):
        with path.open(newline='') as fh:
            return list(csv.DictReader(fh, delimiter='\t'))

    def prepare(self):
        out = self.root / 'input/rrmg'
        target = self.root / 'ref/target.fna'
        self.run_tool('rrmg_prepare_manifest.py', '--input', REPO / 'examples/taxa_demo.tsv',
                      '--outdir', out, '--topN', '2', '--target-name', 'Demo target',
                      '--target-taxid', '330879', '--target-genome-id', 'Afum', '--target-fasta', target)
        path = out / 'community_manifest.tsv'
        rows = self.table(path)
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[-1]['role'], 'target')
        self.assertEqual(float(rows[-1]['rel_abund']), 0)
        self.assertAlmostEqual(sum(float(r['rel_abund']) for r in rows), 1)
        return path, rows

    def fake_datasets(self, fail=False):
        bindir = self.root / 'bin'
        bindir.mkdir(exist_ok=True)
        script = bindir / 'datasets'
        script.write_text('#!' + sys.executable + '\n' + (
            'import sys\nsys.exit(3)\n' if fail else
            'import json, sys, zipfile\n'
            'if sys.argv[1] == "rehydrate": sys.exit(0)\n'
            'assert sys.argv[1:4] == ["download", "genome", "taxon"]\n'
            'out = sys.argv[sys.argv.index("--filename") + 1]\n'
            'with zipfile.ZipFile(out, "w") as z:\n'
            '    z.writestr("ncbi_dataset/data/assembly_data_report.jsonl", json.dumps({"accession":"GCF_TEST.1", "assemblyInfo":{"assemblyLevel":"Complete Genome"}}) + "\\n")\n'
            '    z.writestr("ncbi_dataset/data/GCF_TEST.1/test_genomic.fna", ">test\\n" + "ACGT" * 2500 + "\\n")\n'))
        script.chmod(0o755)
        self.env['PATH'] = str(bindir) + os.pathsep + self.env['PATH']

    def test_prepare_download_and_plans_offline(self):
        path, _ = self.prepare()
        self.fake_datasets()
        downloaded = path.with_name('community_manifest.downloaded.tsv')
        self.run_tool('rrmg_download_genomes.sh', '--manifest', path, '--out-manifest', downloaded)
        rows = self.table(downloaded)
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[-1]['role'], 'target')  # no trailing CR from shell parsing
        for r in rows:
            self.assertTrue(Path(r['fasta']).read_text().startswith('>test'))
        # Rerun with the client failing: existing FASTAs must be reused.
        self.fake_datasets(fail=True)
        self.run_tool('rrmg_download_genomes.sh', '--manifest', path, '--out-manifest', downloaded)
        plan = self.root / 'plan.tsv'
        self.run_tool('rrmg_make_plan.py', '--manifest', downloaded, '--out', plan,
                      '--sample-prefix', 'Pilot', '--target-name', 'Demo target',
                      '--target-taxid', '330879', '--target-genome-id', 'Afum',
                      '--target-fasta', rows[-1]['fasta'], '--target-abundances', '0.01',
                      '--total-pairs', '10000', '--seeds', '1-2', '--readlen', '150')
        self.assertEqual([r['run_index'] for r in self.table(plan)], ['1', '2'])
        bg = self.root / 'temp/rrmg_background'
        self.run_tool('rrmg_make_background_plan.py', '--plan', plan, '-w', self.root, '--outdir', bg)
        self.assertEqual(len(self.table(bg / 'background_plan.tsv')), 1)

    def test_target_download_failure_is_not_success(self):
        path, rows = self.prepare()
        for r in rows[:-1]:
            fa = Path(r['fasta'])
            fa.parent.mkdir(parents=True, exist_ok=True)
            fa.write_text('>existing\nACGT\n')
        self.fake_datasets(fail=True)
        out = self.root / 'downloaded.tsv'
        result = self.run_tool('rrmg_download_genomes.sh', '--manifest', path, '--out-manifest', out, ok=False)
        self.assertIn('Target genome download failed', result.stderr)
        self.assertFalse(out.exists())

    def test_java_active_path_symlink_and_explicit_override(self):
        active = self.root / 'active'
        real = active / 'lib/jvm/bin/java'
        real.parent.mkdir(parents=True)
        real.write_text('#!/bin/bash\nexit 0\n')
        real.chmod(0o755)
        (active / 'bin').mkdir()
        (active / 'bin/java').symlink_to(real)
        self.env['PATH'] = str(active / 'bin') + os.pathsep + self.env['PATH']
        helper = shlex.quote(str(REPO / 'rrmg_environment.sh'))
        command = 'set -eu; source ' + helper + '; printf "%s\\n%s\\n" "$JAVA_HOME" "$(command -v java)"'
        r = subprocess.run(['bash', '-c', command], env=self.env, text=True, capture_output=True, check=True)
        self.assertEqual(r.stdout.splitlines(), [str(active / 'lib/jvm'), str(active / 'bin/java')])
        other = self.root / 'explicit/bin/java'
        other.parent.mkdir(parents=True)
        other.write_text('#!/bin/bash\nexit 0\n')
        other.chmod(0o755)
        self.env['JAVA_BIN'] = str(other)
        r = subprocess.run(['bash', '-c', command], env=self.env, text=True, capture_output=True, check=True)
        self.assertEqual(r.stdout.splitlines(), [str(other.parent.parent), str(other)])


if __name__ == '__main__':
    unittest.main()
