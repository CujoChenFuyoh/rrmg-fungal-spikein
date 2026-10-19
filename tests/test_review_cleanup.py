"""Review regressions: fail-closed plans and portable launchers; no real mapping."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]


@unittest.skipUnless(sys.platform.startswith('linux'), 'Linux shell utilities required')
class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='rrmg review ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.env = dict(os.environ)
        self.env.pop('RRMG_CONFIG', None)
        self.fa = self.root / 'toy.fasta'
        self.fa.write_text('>toy\nACGTACGTACGT\n')
        self.manifest = self.root / 'manifest.tsv'
        self.manifest.write_text('genome_ID\tname\ttaxid\trel_abund\tfasta\trole\n'
                                 f'B\tbackground\t1\t1\t{self.fa}\tbackground\n'
                                 f'T\ttarget\t2\t0\t{self.fa}\ttarget\n')

    def plan_command(self, *extra):
        return [sys.executable, str(REPO / 'rrmg_make_plan.py'), '--manifest', str(self.manifest),
                '--out', str(self.root / 'plan.tsv'), '--sample-prefix', 'Toy',
                '--target-name', 'target', '--target-taxid', '2', '--target-genome-id', 'T',
                '--target-fasta', str(self.fa), '--target-abundances', '0.01',
                '--total-pairs', '10000', '--seeds', '1-2', *extra]

    def invoke(self, command, ok=True, env=None):
        r = subprocess.run(command, cwd=self.root, env=env or self.env,
                           capture_output=True, text=True)
        self.assertEqual(r.returncode == 0, ok, r.stdout + r.stderr)
        return r

    def test_duplicate_seeds_and_designs_rejected_before_plan_write(self):
        for extra in (['--seeds', '1-2,2'], ['--target-abundances', '0.01,0.01'],
                      ['--total-pairs', '10000,10000']):
            r = self.invoke(self.plan_command(*extra), ok=False)
            self.assertIn('Duplicate', r.stdout + r.stderr)
            self.assertFalse((self.root / 'plan.tsv').exists())

    def test_read_fraction_expectations_remain_unchanged(self):
        self.invoke(self.plan_command())
        with (self.root / 'plan.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual([r['run_index'] for r in rows], ['1', '2'])
        self.assertEqual(len({r['run_id'] for r in rows}), 2)
        self.assertEqual([float(r['target_expected_pairs']) for r in rows], [100, 100])

    def test_copy_number_refused_in_planning_background_and_step2(self):
        r = self.invoke(self.plan_command('--abundance-mode', 'copy_number'), ok=False)
        self.assertIn('copy_number', r.stderr)
        self.assertFalse((self.root / 'plan.tsv').exists())
        self.invoke(self.plan_command())
        plan = self.root / 'plan.tsv'
        plan.write_text(plan.read_text().replace('read_fraction', 'copy_number'))
        r = self.invoke([sys.executable, str(REPO / 'rrmg_make_background_plan.py'),
                         '--plan', str(plan), '-w', str(self.root)], ok=False)
        self.assertIn('Only read_fraction', r.stdout + r.stderr)
        index, tool = self.root / 'background_index.tsv', self.root / 'randomreads.sh'
        index.write_text('header\n')
        tool.write_text('# fixture\n')
        env = dict(self.env, RANDOMREADS_SH=str(tool))
        r = self.invoke(['bash', str(REPO / 'rrmg_step2_make_target_and_merge.sh'), str(plan), '1',
                         '-w', str(self.root), '--background-index', str(index)], ok=False, env=env)
        self.assertIn('Only read_fraction', r.stdout + r.stderr)

    def test_duplicate_run_id_stops_dry_run_before_work(self):
        plan = self.root / 'plan.tsv'
        plan.write_text('run_index\trun_id\n1\tSame\n2\tSame\n')
        config = self.root / 'rrmg_config.sh'
        config.write_text(f'PIPELINE_DIR={shlex.quote(str(REPO))}\nWD={shlex.quote(str(self.root))}\n'
                          f'PLAN={shlex.quote(str(plan))}\n')
        env = dict(self.env, RRMG_CONFIG=str(config))
        r = self.invoke(['bash', str(REPO / 'rrmg_submit.sh'), 'step2', '--dry-run'], ok=False, env=env)
        self.assertIn('Duplicate run_id', r.stdout + r.stderr)
        self.assertFalse((self.root / 'logs').exists())

    def test_vamb_tool_check_does_not_require_exported_conda_name(self):
        text = (REPO / 'rrmg_step3b_vamb_binning.sh').read_text()
        block = text[text.index('ACTIVE_ENV='):text.index('VAMB_VERSION=')]
        bins = self.root / 'bin'
        bins.mkdir()
        for name in ('minimap2', 'samtools', 'vamb'):
            tool = bins / name
            tool.write_text('#!/bin/bash\nexit 0\n')
            tool.chmod(0o755)
        env = dict(self.env, PATH=str(bins) + os.pathsep + self.env['PATH'])
        env.pop('RRMG_CONDA_ENV', None)
        env.pop('CONDA_DEFAULT_ENV', None)
        r = self.invoke(['bash', '-uc', 'VAMB_COMMAND=vamb\n' + block + '\nprintf "%s" "$ACTIVE_ENV"'], env=env)
        self.assertEqual(r.stdout, 'active-PATH')

    def test_windows_r_path_detection_and_optional_library(self):
        bins = self.root / 'bin'
        bins.mkdir()
        result = self.root / 'r_call.json'
        rscript = bins / 'Rscript.exe'
        rscript.write_text('#!' + sys.executable + '\nimport json, os, sys\n'
                           'with open(os.environ["R_CALL"],"w") as f:\n'
                           ' json.dump({"args":sys.argv[1:],"lib":os.getenv("R_LIBS_USER"),"wslenv":os.getenv("WSLENV")},f)\n')
        rscript.chmod(0o755)
        wslpath = bins / 'wslpath'
        wslpath.write_text('#!/bin/bash\nprintf "Z:%s\\n" "$2"\n')
        wslpath.chmod(0o755)
        env = dict(self.env, PATH=str(bins) + os.pathsep + self.env['PATH'],
                   R_CALL=str(result), R_LIBS_USER='/linux/library')
        env.pop('WIN_RSCRIPT', None)
        env.pop('WIN_R_LIB', None)
        command = ['bash', str(REPO / 'rrmg_windows_r.sh'), '--check']
        self.invoke(command, env=env)
        self.assertIsNone(json.loads(result.read_text())['lib'])
        env.update(WIN_RSCRIPT=str(rscript), WIN_R_LIB='C:/R libraries/user')
        self.invoke(command, env=env)
        call = json.loads(result.read_text())
        self.assertEqual(call['lib'], 'C:/R libraries/user')
        self.assertIn('R_LIBS_USER/w', call['wslenv'])
        fake_r = self.root / 'plot script.R'
        fake_r.write_text('# fixture\n')
        self.invoke(['bash', str(REPO / 'rrmg_windows_r.sh'), str(fake_r), str(self.root / 'data dir')], env=env)
        self.assertEqual(json.loads(result.read_text())['args'][1], 'Z:' + str(fake_r))

    def test_paf_defaults_and_public_runbook_are_preserved(self):
        self.assertIn('MAKE_PAF=1', (REPO / 'rrmg_step4_reads_to_gene_depth.sh').read_text())
        runbook = (REPO / 'runbook.sh').read_text()
        active = '\n'.join(line for line in runbook.splitlines() if not line.lstrip().startswith('#'))
        self.assertNotIn('--no-paf', active)
        for private in ('q6321', '/mnt/e/', '/mnt/f/', 'E:/github_repo', 'H:/genedata'):
            self.assertNotIn(private, runbook)
            self.assertNotIn(private, (REPO / 'rrmg_windows_r.sh').read_text())
        self.assertIn('step3 1-2 -t 4 -m 16 --jobs 1\n', runbook)
        self.assertFalse((REPO / 'legacy/rrmg_add_target_only_cyp51a_to_figure_data.py').exists())
        names = (REPO / 'src/rrmg_cli/core_files.txt').read_text().splitlines()
        self.assertNotIn('rrmg_add_target_only_cyp51a_to_figure_data.py', names)

    def test_removed_wrappers_are_absent_and_generic_templates_remain(self):
        names = (REPO / 'src/rrmg_cli/core_files.txt').read_text().splitlines()
        for name in ('rrmg_step1_background_array.pbs', 'rrmg_step2_target_merge_array.pbs',
                     'rrmg_step3_downstream_array.pbs', 'rrmg_step3b_vamb_array.pbs',
                     'rrmg_step4_5_supplementary_downstream_array.pbs'):
            self.assertFalse((REPO / name).exists())
            self.assertNotIn(name, names)
        for name in ('pbs_example/pbs.pbs', 'pbs_example/pbs_array.pbs', 'rrmg_prepare_alignment_inputs.py'):
            self.assertTrue((REPO / name).is_file())

    def test_build_cannot_reintroduce_removed_resources_from_old_cache(self):
        build = self.root / 'build'
        assets = build / 'rrmg_cli/_workflow'
        assets.mkdir(parents=True)
        obsolete = assets / 'rrmg_step1_background_array.pbs'
        obsolete.write_text('# obsolete generated resource\n')
        old_patch = assets / 'rrmg_add_target_only_cyp51a_to_figure_data.py'
        old_patch.write_text('# obsolete generated resource\n')
        result = subprocess.run([sys.executable, 'setup.py', 'build_py', '--build-lib', str(build)],
                                cwd=REPO, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse(obsolete.exists())
        self.assertFalse(old_patch.exists())
        for name in ('pbs_example/pbs.pbs', 'pbs_example/pbs_array.pbs', 'rrmg_prepare_alignment_inputs.py'):
            self.assertTrue((assets / name).is_file())


if __name__ == '__main__':
    unittest.main()
