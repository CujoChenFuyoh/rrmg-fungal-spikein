"""Unit tests for installation/dispatch. No external analyses or real PBS jobs."""
import contextlib
import io
import json
import os
from pathlib import Path
import shlex
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrmg_cli import cli


class CLITests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="rrmg-cli-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.out = io.StringIO()
        self.addCleanup(self.out.close)

    def call(self, *args):
        with contextlib.redirect_stdout(self.out), contextlib.redirect_stderr(self.out):
            return cli.main(list(map(str, args)))

    def test_core_allowlist_is_complete_and_excludes_private_data(self):
        source = cli.core_root()
        names = cli.core_names()
        self.assertEqual(len(names), len(set(names)))
        self.assertNotIn("rrmg_config.sh", names)
        self.assertIn("rrmg_run_local.sh", names)
        for p in source.glob("rrmg_*"):
            if p.suffix in (".sh", ".pbs", ".py") and p.name != "rrmg_config.sh":
                self.assertIn(p.name, names)
        self.assertEqual((source / "VERSION").read_text().strip(), cli.__version__)

    def test_export_preserves_original_core_bytes(self):
        target = self.root / "export"
        self.assertEqual(self.call("export", target), 0)
        expected = cli.inventory(cli.core_root())
        self.assertEqual(cli.inventory(target), expected)
        self.assertEqual(json.loads((target / "rrmg_core_manifest.json").read_text())["sha256"], expected)
        self.assertFalse((target / "rrmg_config.sh").exists())

    def test_export_never_merges_existing_directory(self):
        target = self.root / "export"
        target.mkdir()
        marker = target / "keep.txt"
        marker.write_text("keep")
        self.assertEqual(self.call("export", target), 1)
        self.assertEqual(marker.read_text(), "keep")
        self.assertEqual(len(list(target.iterdir())), 1)

    def test_init_keeps_input_plan_results_and_config(self):
        wd = self.root / "existing project"
        p = wd / "input/rrmg/rrmg_permutation_plan.tsv"
        p.parent.mkdir(parents=True)
        p.write_text("original plan\n")
        with patch.object(cli, "detected_settings", return_value={"RRMG_CONDA_ENV": "fungamr_meta"}):
            self.assertEqual(self.call("init", "--workdir", wd), 0)
            config = wd / "rrmg_config.sh"
            before = config.read_bytes()
            self.assertEqual(self.call("init", "--workdir", wd), 1)
        self.assertEqual(config.read_bytes(), before)
        self.assertEqual(p.read_text(), "original plan\n")
        self.assertFalse((wd / "temp").exists())
        self.assertFalse((wd / "result").exists())
        self.assertTrue((wd / ".rrmg" / ("core-" + cli.__version__) / "rrmg_run_local.sh").is_file())

    def test_init_can_bind_existing_original_checkout(self):
        wd = self.root / "work"
        with patch.object(cli, "detected_settings", return_value={}):
            self.assertEqual(self.call("init", "--workdir", wd, "--pipeline-dir", cli.core_root()), 0)
        self.assertFalse((wd / ".rrmg").exists())
        text = (wd / "rrmg_config.sh").read_text()
        self.assertIn("PIPELINE_DIR=" + shlex.quote(str(cli.core_root())), text)

    def test_init_quotes_spaces_and_shell_metacharacters(self):
        wd = self.root / "a b's;$literal"
        with patch.object(cli, "detected_settings", return_value={}):
            self.assertEqual(self.call("init", "--workdir", wd, "--pipeline-dir", cli.core_root()), 0)
        line = next(s for s in (wd / "rrmg_config.sh").read_text().splitlines() if s.startswith("WD="))
        self.assertEqual(shlex.split(line), ["WD=" + str(wd.resolve())])

    def test_unresolved_paths_are_visible_not_invented(self):
        args = cli.parser().parse_args(["init", "--workdir", str(self.root)])
        with patch.dict(os.environ, {}, clear=True), patch.object(cli.shutil, "which", return_value=None):
            settings = cli.detected_settings(args)
        self.assertIn("/EDIT_ME/", settings["CONDA_SH"])
        self.assertIn("/EDIT_ME/", settings["JAVA_BIN"])

    def test_explicit_tool_paths_are_retained(self):
        args = cli.parser().parse_args(["init", "--workdir", str(self.root), "--conda-sh", str(self.root / "conda.sh"),
                                        "--randomreadsmg", str(self.root / "bbmap/randomreadsmg.sh"),
                                        "--java-bin", str(self.root / "java/bin/java")])
        settings = cli.detected_settings(args)
        self.assertEqual(settings["BBMAP_DIR"], str((self.root / "bbmap").resolve()))
        self.assertEqual(settings["JAVA_HOME"], str((self.root / "java").resolve()))

    def test_run_dispatches_original_runner_without_qsub(self):
        config = self.root / "site config.sh"
        config.write_text("# test\n")
        with patch.object(cli.sys, "platform", "linux"), patch.object(cli.shutil, "which", return_value="/bin/bash"), \
             patch.object(cli.os, "execvpe") as run:
            self.assertEqual(self.call("run", "step2", "451-480", "--jobs", "1", "--dry-run", "--config", config), 0)
        name, argv, env = run.call_args.args
        self.assertEqual(name, "bash")
        self.assertEqual(argv[4:], [str(config.resolve()), "rrmg_run_local.sh", "step2", "451-480", "--jobs", "1", "--dry-run"])
        self.assertEqual(env["RRMG_CONFIG"], str(config.resolve()))
        self.assertNotIn("qsub", argv)
        self.assertNotIn("--force", argv)

    def test_force_is_never_implicit(self):
        args = cli.parser().parse_args(["run", "step3", "1"])
        self.assertFalse(args.force)
        self.assertTrue(cli.parser().parse_args(["run", "step3", "1", "--force"]).force)

    def test_config_priority_is_explicit_then_env_then_workdir(self):
        first, second = self.root / "first.sh", self.root / "second.sh"
        first.write_text("first")
        second.write_text("second")
        args = cli.parser().parse_args(["run", "step3", "1", "--config", str(first)])
        with patch.dict(os.environ, {"RRMG_CONFIG": str(second)}):
            self.assertEqual(cli.config_path(args), first.resolve())
            args.config = None
            self.assertEqual(cli.config_path(args), second.resolve())
        args.workdir = self.root
        (self.root / "rrmg_config.sh").write_text("local")
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(cli.config_path(args), self.root / "rrmg_config.sh")

    def test_invalid_selection_fails_before_launch(self):
        config = self.root / "config.sh"
        config.write_text("# test")
        with patch.object(cli.sys, "platform", "linux"), patch.object(cli.shutil, "which", return_value="bash"), \
             patch.object(cli.os, "execvpe") as execute:
            self.assertEqual(self.call("run", "step3", "1,2", "--config", config), 1)
            self.assertEqual(self.call("run", "step3", "1", "--jobs", "0", "--config", config), 1)
        execute.assert_not_called()

    def test_script_arguments_and_failure_status_are_preserved(self):
        with patch.object(cli.subprocess, "call", return_value=7) as call:
            self.assertEqual(self.call("script", "rrmg_make_plan.py", "--", "--manifest", "a b.tsv"), 7)
        self.assertEqual(call.call_args.args[0][-2:], ["--manifest", "a b.tsv"])

    def test_script_cannot_escape_allowlist(self):
        with patch.object(cli.subprocess, "call") as call:
            self.assertEqual(self.call("script", "../../foreign.py"), 1)
        call.assert_not_called()

    def test_core_lookup_reports_broken_installation(self):
        with patch.object(cli, "core_names", return_value=["NOT_A_REAL_CORE_FILE"]):
            with self.assertRaises(ValueError):
                cli.core_root()

    def test_doctor_dispatches_original_preflight(self):
        config = self.root / "site.sh"
        config.write_text("# test")
        with patch.object(cli.sys, "platform", "linux"), patch.object(cli.shutil, "which", return_value="bash"), \
             patch.object(cli.os, "execvpe") as execute:
            self.assertEqual(self.call("doctor", "--config", config), 0)
        self.assertEqual(execute.call_args.args[1][-1], "rrmg_preflight.sh")


if __name__ == "__main__":
    unittest.main()
