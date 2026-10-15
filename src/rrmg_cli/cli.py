"""Portable installation and thin launchers, not a replacement workflow engine."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

from . import __version__

STEPS = ("step1", "step2", "step3", "step3b", "step4", "step5", "step45")


def core_names():
    path = Path(__file__).with_name("core_files.txt")
    return [s for s in path.read_text(encoding="utf-8").splitlines() if s and not s.startswith("#")]


def core_root():
    package = Path(__file__).resolve().parent
    candidates = (package / "_workflow", package.parents[1])
    for candidate in candidates:
        if all((candidate / name).is_file() for name in core_names()):
            return candidate
    raise ValueError("Core scripts missing from installation; reinstall a complete wheel or checkout")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda: f.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def inventory(root):
    return {name: digest(root / name) for name in core_names()}


def export_core(destination):
    """Copy an allowlisted, self-contained set; never merge with existing code."""
    destination = destination.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError("Destination already exists; choose a NEW code directory: " + str(destination))
    source = core_root()
    expected = inventory(source)
    destination.mkdir(parents=True)
    for name in core_names():
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / name, target)
        if target.suffix in (".sh", ".pbs"):
            target.chmod(target.stat().st_mode | 0o111)
    if inventory(destination) != expected:
        raise ValueError("Export checksum mismatch; incomplete directory retained for inspection")
    record = dict(version=__version__, created_utc=datetime.now(timezone.utc).isoformat(), sha256=expected)
    with (destination / "rrmg_core_manifest.json").open("x", encoding="utf-8") as f:
        json.dump(record, f, indent=2)
    return destination.resolve()


def resolved(value):
    return str(Path(value).expanduser().resolve())


def detected_settings(args):
    """Detect only real local paths; unresolved placeholders stay conspicuous."""
    conda_sh = args.conda_sh
    if not conda_sh and os.environ.get("CONDA_EXE"):
        candidate = Path(os.environ["CONDA_EXE"]).resolve().parent.parent / "etc/profile.d/conda.sh"
        if candidate.is_file():
            conda_sh = str(candidate)
    if not conda_sh and shutil.which("conda"):
        try:
            result = subprocess.run(["conda", "info", "--base"], capture_output=True, text=True, timeout=15, check=True)
            candidate = Path(result.stdout.strip()) / "etc/profile.d/conda.sh"
            if candidate.is_file():
                conda_sh = str(candidate)
        except (OSError, subprocess.SubprocessError):
            pass
    # Defaults stay portable: resolve active PATH at runtime, not at installation.
    # Explicit overrides are local deployment choices, never release defaults.
    randomreadsmg = args.randomreadsmg
    randomreads = args.randomreads
    if args.bbmap_dir:
        randomreadsmg = randomreadsmg or str(Path(args.bbmap_dir).expanduser() / "randomreadsmg.sh")
        randomreads = randomreads or str(Path(args.bbmap_dir).expanduser() / "randomreads.sh")
    if not randomreads and randomreadsmg:
        randomreads = str(Path(randomreadsmg).expanduser().parent / "randomreads.sh")
    if not randomreadsmg and randomreads:
        randomreadsmg = str(Path(randomreads).expanduser().parent / "randomreadsmg.sh")
    bbmap = args.bbmap_dir or (str(Path(randomreadsmg).resolve().parent) if randomreadsmg else None)
    java = args.java_bin
    settings = {
        "CONDA_SH": resolved(conda_sh) if conda_sh else "",
        "RRMG_CONDA_ENV": args.env,
        "BBMAP_DIR": resolved(bbmap) if bbmap else "",
        "RANDOMREADSMG_SH": resolved(randomreadsmg) if randomreadsmg else "",
        "RANDOMREADS_SH": resolved(randomreads) if randomreads else "",
        "JAVA_HOME": str(Path(java).resolve().parent.parent) if java else "",
        "JAVA_BIN": resolved(java) if java else "",
    }
    return settings


def init_project(args):
    wd = args.workdir.expanduser().resolve()
    config = args.config.expanduser().absolute() if args.config else wd / "rrmg_config.sh"
    if config.exists() or config.is_symlink():
        raise ValueError("Configuration exists; keep it and use rrmg run --config " + str(config))
    if args.pipeline_dir:
        pipeline = args.pipeline_dir.expanduser().resolve()
        for name in ("rrmg_run_local.sh", "rrmg_config.example.sh", "rrmg_preflight.sh"):
            if not (pipeline / name).is_file():
                raise ValueError("Not an existing RRMG code directory: " + str(pipeline))
    else:
        pipeline = wd / ".rrmg" / ("core-" + __version__)
        if pipeline.exists() or pipeline.is_symlink():
            raise ValueError("Core snapshot exists; do not overwrite it: " + str(pipeline))
    settings = detected_settings(args)
    settings.update(PIPELINE_DIR=str(pipeline), WD=str(wd))
    template = (core_root() / "rrmg_config.example.sh").read_text(encoding="utf-8")
    for key, value in settings.items():
        if key in ("JAVA_BIN", "JAVA_HOME", "BBMAP_DIR", "RANDOMREADSMG_SH", "RANDOMREADS_SH") and not value:
            # Preserve template inheritance so PBS-job exports remain effective.
            continue
        if "\n" in value or "\r" in value:
            raise ValueError("Newlines not allowed in configuration paths/names")
        template, n = re.subn(r"^" + key + r"=.*$", lambda _: key + "=" + shlex.quote(value), template, flags=re.M)
        if n != 1:
            raise ValueError("Unexpected configuration template key: " + key)
    wd.mkdir(parents=True, exist_ok=True)
    if not args.pipeline_dir:
        export_core(pipeline)
    config.parent.mkdir(parents=True, exist_ok=True)
    with config.open("x", encoding="utf-8", newline="\n") as f:
        f.write("# Initialized by rrmg " + __version__ + "; input data and plans were not modified.\n")
        f.write(template)
    print("[CONFIG] " + str(config))
    print("[CORE]   " + str(pipeline))
    for key, value in settings.items():
        if "/EDIT_ME/" in value:
            print("[EDIT]   " + key + "=" + value)
    print("[NEXT] Review paths/resources, then: rrmg doctor --config " + shlex.quote(str(config)))
    print("[NOTE] No plan, reads, reference, environment or analysis result was created/replaced.")


def config_path(args):
    value = args.config or os.environ.get("RRMG_CONFIG")
    path = Path(value).expanduser().resolve() if value else (args.workdir or Path.cwd()).expanduser().resolve() / "rrmg_config.sh"
    if not path.is_file():
        raise ValueError("Configuration not found: {}. Use --config or rrmg init.".format(path))
    return path


def bash_command(config, script, arguments):
    # Config is trusted, executable shell code, just as in the original scripts.
    # All user arguments are separate argv entries, never interpolated in code.
    program = ('set -e\nsource "$1"\n'
               ': "${PIPELINE_DIR:?Set PIPELINE_DIR in the configuration}"\n'
               'exec bash "$PIPELINE_DIR/$2" "${@:3}"')
    return ["bash", "-c", program, "rrmg", str(config), script, *arguments]


def launch(args):
    if not sys.platform.startswith("linux"):
        raise ValueError("Analysis/doctor requires Linux (or WSL); export/init/help work on other platforms")
    config = config_path(args)
    if not shutil.which("bash"):
        raise ValueError("GNU Bash is required")
    script = "rrmg_preflight.sh" if args.command == "doctor" else "rrmg_run_local.sh"
    forwarded = []
    if args.command == "run":
        if args.rows and not re.fullmatch(r"[1-9][0-9]*(?:-[1-9][0-9]*)?", args.rows):
            raise ValueError("Use an original row index or a START-END range")
        if args.jobs is not None and args.jobs < 1:
            raise ValueError("--jobs must be positive")
        if args.threads is not None and args.threads < 1:
            raise ValueError("--threads must be positive")
        if args.memory is not None and args.memory < 1:
            raise ValueError("--memory must be a positive integer in GiB")
        if args.memory is not None and args.step not in ("step1", "step2", "step3"):
            raise ValueError("--memory is supported for step1/2/3 only; set scheduler memory for this step")
        forwarded = [args.step]
        if args.rows:
            forwarded.append(args.rows)
        if args.threads is not None:
            forwarded += ["--threads", str(args.threads)]
        if args.memory is not None:
            forwarded += ["--memory", str(args.memory)]
        if args.jobs is not None:
            forwarded += ["--jobs", str(args.jobs)]
        if args.dry_run:
            forwarded.append("--dry-run")
        if args.force:
            forwarded.append("--force")
    env = dict(os.environ, RRMG_CONFIG=str(config))
    # exec keeps PBS signals/exit codes attached to the original runner.
    command = bash_command(config, script, forwarded)
    os.execvpe(command[0], command, env)


def run_script(args):
    """Planning/collection utilities retain their original arguments and cwd."""
    if args.name not in core_names() or Path(args.name).suffix not in (".py", ".sh", ".pbs"):
        raise ValueError("Unknown core script; run rrmg path to inspect available scripts")
    path = core_root() / args.name
    forwarded = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
    if path.suffix == ".py":
        return subprocess.call([sys.executable, str(path), *forwarded])
    if not sys.platform.startswith("linux"):
        raise ValueError("Shell/PBS scripts require Linux (or WSL)")
    # Running a .pbs file with bash does not ask PBS to submit a job.
    return subprocess.call(["bash", str(path), *forwarded])


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--version", action="version", version="rrmg " + __version__)
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("init", help="Write site config; preserve all existing data/configuration")
    q.add_argument("--workdir", type=Path, required=True)
    q.add_argument("--config", type=Path, help="New config filename (default WORKDIR/rrmg_config.sh)")
    q.add_argument("--pipeline-dir", type=Path, help="Use existing checkout; otherwise export a versioned core snapshot")
    q.add_argument("--env", default="fungamr_meta")
    for flag in ("conda-sh", "bbmap-dir", "randomreadsmg", "randomreads", "java-bin"):
        q.add_argument("--" + flag, help="Override auto-detected local path")
    q = sub.add_parser("export", help="Copy the original core scripts to a NEW directory")
    q.add_argument("destination", type=Path)
    sub.add_parser("path", help="Print installed core script directory")
    for name in ("run", "doctor"):
        q = sub.add_parser(name, help="Call original " + ("local runner (no qsub)" if name == "run" else "preflight"))
        q.add_argument("--config", type=Path)
        q.add_argument("--workdir", type=Path, help="Directory containing rrmg_config.sh")
        if name == "run":
            q.add_argument("step", choices=STEPS)
            q.add_argument("rows", nargs="?", help="Original row or START-END; omitted means all rows")
            q.add_argument("-t", "--threads", type=int, help="Threads per row; overrides the step config")
            q.add_argument("-m", "--memory", type=int, metavar="GB", help="Per-row GiB for step1/2 Java or step3 SPAdes")
            q.add_argument("--jobs", type=int)
            q.add_argument("--dry-run", action="store_true")
            q.add_argument("--force", action="store_true", help="Original destructive overwrite behavior; use deliberately")
    q = sub.add_parser("script", help="Run one original utility with its original arguments")
    q.add_argument("name")
    q.add_argument("arguments", nargs=argparse.REMAINDER)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == "init":
            init_project(args)
        elif args.command == "export":
            print(export_core(args.destination))
        elif args.command == "path":
            print(core_root())
        elif args.command in ("run", "doctor"):
            launch(args)
        elif args.command == "script":
            return run_script(args)
        return 0
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print("[ERROR] " + str(error), file=sys.stderr)
        return 1
