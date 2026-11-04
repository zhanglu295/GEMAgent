"""Small, dependency-free run manifest and hashing helpers."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256(path: str | os.PathLike[str]) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _json_default(value: Any):
    if isinstance(value, Path):
        return str(value)
    return str(value)


def write_manifest(workspace: str | os.PathLike[str], **payload: Any) -> Path:
    ws = Path(workspace)
    ws.mkdir(parents=True, exist_ok=True)
    data = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        **payload,
    }
    out = ws / "run_manifest.json"
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(data, indent=2, default=_json_default), encoding="utf-8")
    tmp.replace(out)
    return out


def read_manifest(workspace: str | os.PathLike[str]) -> dict:
    return json.loads((Path(workspace) / "run_manifest.json").read_text(encoding="utf-8"))