"""
Thin wrapper around the ORIGINAL clean_sample.py and classify.py (kept unchanged).

The GUI calls run_pipeline(); this module runs the two original scripts in a private
temporary folder, using the bundled references.fasta, targets.txt and BLAST+ programs,
and returns the results as plain Python objects.
"""
import contextlib
import csv
import io
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path


class PipelineError(Exception):
    """A problem the user can understand (bad file, missing BLAST, ...)."""


def resource_dir() -> Path:
    """Folder that holds data/ and blast_bin/ (the repo root, or the PyInstaller bundle)."""
    if getattr(sys, "frozen", False):
        return Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
    return Path(__file__).resolve().parent.parent


# --- Make the original scripts' BLAST calls quiet on Windows --------------------------
# classify.py calls subprocess.run(...) directly. In a windowed app that would flash a
# console window for every BLAST call, so we add CREATE_NO_WINDOW and a closed stdin.
_original_run = subprocess.run


def _quiet_run(*args, **kwargs):
    if sys.platform == "win32":
        kwargs["creationflags"] = kwargs.get("creationflags", 0) | subprocess.CREATE_NO_WINDOW
    if "input" not in kwargs:
        kwargs.setdefault("stdin", subprocess.DEVNULL)
    return _original_run(*args, **kwargs)


subprocess.run = _quiet_run

import clean_sample  # noqa: E402  (original script, unchanged)
import classify  # noqa: E402  (original script, unchanged)


@dataclass
class PipelineResult:
    sample_name: str
    header: list
    rows: list
    counts: list  # [(label, n), ...] most common first
    total: int
    clean_log: str = ""
    classify_log: str = ""
    csv_text: str = field(default="", repr=False)


def _prepare_blast_path():
    blast_dir = resource_dir() / "blast_bin"
    if blast_dir.is_dir():
        os.environ["PATH"] = str(blast_dir) + os.pathsep + os.environ.get("PATH", "")
    missing = [t for t in ("blastn", "makeblastdb") if shutil.which(t) is None]
    if missing:
        raise PipelineError(
            "BLAST+ programs not found: " + ", ".join(missing)
            + ". The app installation looks incomplete."
        )


def _run_script(main_fn, argv, what):
    """Run an original script's main() as if started from the command line."""
    buf = io.StringIO()
    old_argv = sys.argv
    sys.argv = argv
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            try:
                main_fn()
            except SystemExit as e:
                if e.code not in (0, None):
                    msg = e.code if isinstance(e.code, str) else buf.getvalue().strip()
                    raise PipelineError(f"{what}: {msg}") from None
            except subprocess.CalledProcessError as e:
                detail = e.stderr.strip() if isinstance(e.stderr, str) else ""
                raise PipelineError(f"{what}: BLAST failed. {detail}") from None
    finally:
        sys.argv = old_argv
    return buf.getvalue()


def _input_without_bom(sample_path: Path, work: Path) -> Path:
    """Files saved by Windows Notepad may start with a hidden BOM that breaks the cleaner."""
    with open(sample_path, "rb") as f:
        head = f.read(3)
    if head != b"\xef\xbb\xbf":
        return sample_path
    fixed = work / "input_no_bom.txt"
    with open(sample_path, "rb") as src, open(fixed, "wb") as dst:
        src.read(3)
        shutil.copyfileobj(src, dst)
    return fixed


def run_pipeline(sample_path, status=lambda text: None) -> PipelineResult:
    sample_path = Path(sample_path).resolve()  # absolute, because we change folder below
    if not sample_path.is_file():
        raise PipelineError(f"File not found: {sample_path}")

    res = resource_dir()
    refs = res / "data" / "references.fasta"
    targets = res / "data" / "targets.txt"
    for p in (refs, targets):
        if not p.is_file():
            raise PipelineError(f"Missing app data file: {p.name}. The app installation looks incomplete.")
    _prepare_blast_path()

    work = Path(tempfile.mkdtemp(prefix="stool16s_"))
    old_cwd = os.getcwd()
    try:
        # Work inside the temp folder with short relative names. BLAST mis-reads database
        # paths that contain spaces (e.g. C:\Users\First Last\...), so we avoid long paths.
        shutil.copy(refs, work / "references.fasta")
        shutil.copy(targets, work / "targets.txt")
        os.chdir(work)

        status("Cleaning the sample...")
        source = _input_without_bom(sample_path, work)
        clean_log = _run_script(
            clean_sample.main,
            ["clean_sample.py", "--input", str(source), "--output", "clean_sample.fasta"],
            "Cleaning",
        )
        with open("clean_sample.fasta") as f:
            kept = sum(1 for line in f if line.startswith(">"))
        if kept == 0:
            raise PipelineError(
                "No sequences were left after cleaning, so there is nothing to classify.\n\n"
                + clean_log.strip()
            )

        status(f"Classifying {kept} sequences with BLAST...")
        classify_log = _run_script(
            classify.main,
            ["classify.py", "--targets", "targets.txt", "--sample", "clean_sample.fasta",
             "--out", "results.csv", "--skip-fetch"],
            "Classification",
        )

        with open("results.csv", newline="", encoding="utf-8") as f:
            csv_text = f.read()
        parsed = list(csv.reader(io.StringIO(csv_text, newline="")))
        header, rows = parsed[0], parsed[1:]
        label_col = header.index("label")
        counts = Counter(r[label_col] for r in rows).most_common()
        return PipelineResult(
            sample_name=sample_path.name,
            header=header,
            rows=rows,
            counts=counts,
            total=len(rows),
            clean_log=clean_log,
            classify_log=classify_log,
            csv_text=csv_text,
        )
    finally:
        os.chdir(old_cwd)
        shutil.rmtree(work, ignore_errors=True)
