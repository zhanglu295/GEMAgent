"""Runtime configuration for the standalone gemagent package.

The package intentionally does not vendor gapseq or its databases.  All paths
that vary by installation are resolved from environment variables first.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

DEFAULT_DB_DIR = Path(os.environ.get("GEMAGENT_DB_DIR", "")) if os.environ.get("GEMAGENT_DB_DIR") else None


def env(name: str, default: str = "") -> str:
    return os.environ.get(name, default).strip()


def python() -> str:
    return env("GEMAGENT_PYTHON", sys.executable)


def rscript() -> str:
    candidates = [env("GEMAGENT_RSCRIPT"), env("RSCRIPT")]
    py_bin = Path(sys.executable).resolve().parent
    candidates.append(str(py_bin / "Rscript"))
    candidates.append(shutil.which("Rscript") or "")
    for cand in candidates:
        if cand and os.access(cand, os.X_OK):
            return cand
    return ""


def gapseq() -> str:
    candidates = [env("GAPSEQ_BIN"), shutil.which("gapseq") or ""]
    env_root = Path(sys.executable).resolve().parent.parent
    envs_dir = env_root.parent
    if envs_dir.is_dir():
        for entry in sorted(envs_dir.iterdir()):
            cand = entry / "bin" / "gapseq"
            if cand.is_file() and os.access(cand, os.X_OK):
                candidates.append(str(cand))
    for cand in candidates:
        if cand and os.access(cand, os.X_OK):
            return cand
    return ""


def gapseq_dat_dir() -> str:
    explicit = env("GEMAGENT_GAPSEQ_DAT_DIR")
    if explicit:
        return explicit
    gs = gapseq()
    if gs:
        # .../envs/gapseq/share/gapseq/gapseq -> .../share/gapseq/dat
        cand = Path(gs).resolve().parent / "dat"
        if cand.is_dir():
            return str(cand)
    return ""


def db_dir() -> str:
    explicit = env("GEMAGENT_DB_DIR")
    if explicit:
        return explicit
    # optional repository-local database for advanced users
    local = Path(__file__).resolve().parents[2] / "database"
    return str(local) if local.is_dir() else ""


def configure(db_dir_override: str = "", gapseq_bin: str = "", rscript_bin: str = "",
              eval_root: str = "") -> dict[str, str]:
    """Set the environment consumed by the copied review implementations."""
    if db_dir_override:
        os.environ["GEMAGENT_DB_DIR"] = db_dir_override
    if gapseq_bin:
        os.environ["GAPSEQ_BIN"] = gapseq_bin
    if rscript_bin:
        os.environ["GEMAGENT_RSCRIPT"] = rscript_bin
    if eval_root:
        os.environ["GEMAGENT_EVAL_ROOT"] = eval_root
    gs = gapseq()
    if gs and not env("GAPSEQ_BIN"):
        os.environ["GAPSEQ_BIN"] = gs
    rs = rscript()
    if rs and not env("GEMAGENT_RSCRIPT"):
        os.environ["GEMAGENT_RSCRIPT"] = rs
    dat = gapseq_dat_dir()
    if dat and not env("GEMAGENT_GAPSEQ_DAT_DIR"):
        os.environ["GEMAGENT_GAPSEQ_DAT_DIR"] = dat
    return {
        "python": python(),
        "rscript": rs,
        "gapseq": gs,
        "gapseq_dat_dir": dat,
        "db_dir": db_dir(),
        "eval_root": env("GEMAGENT_EVAL_ROOT"),
    }


def require_tools() -> None:
    got = configure()
    missing = [k for k in ("rscript", "gapseq") if not got[k]]
    if missing:
        raise RuntimeError("Missing required executable(s): " + ", ".join(missing))