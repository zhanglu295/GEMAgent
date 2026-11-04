"""Thin, non-invasive adapters around the installed official gapseq CLI.

No gapseq R source is copied into this package.  The only gapfill-specific
operation is filtering the candidate RDS tables before invoking official
``gapseq fill``.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Iterable

from .. import config


class GapseqError(RuntimeError):
    pass


def _log_tail(path: Path, n: int = 20) -> str:
    if not path.exists():
        return ""
    lines = path.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n:])


def run_gapseq(args: list[str], log_path: str | os.PathLike[str] | None = None,
               dry_run: bool = False, cwd: str | os.PathLike[str] | None = None) -> subprocess.CompletedProcess:
    bin_path = config.gapseq()
    if not bin_path:
        raise GapseqError("gapseq executable not found; set GAPSEQ_BIN or add gapseq to PATH")
    cmd = [bin_path, *map(str, args)]
    if dry_run:
        print("[dry-run] " + " ".join(cmd))
        return subprocess.CompletedProcess(cmd, 0, "", "")
    env = os.environ.copy()
    rs = config.rscript()
    if rs:
        env["PATH"] = str(Path(rs).parent) + os.pathsep + env.get("PATH", "")
        env.setdefault("RSCRIPT", rs)
    log = Path(log_path) if log_path else None
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        fh = open(log, "w", encoding="utf-8")
    else:
        fh = None
    try:
        proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True,
                              env=env,
                              stdout=fh if fh else subprocess.PIPE,
                              stderr=subprocess.STDOUT if fh else subprocess.PIPE,
                              check=False)
    finally:
        if fh:
            fh.close()
    if proc.returncode != 0:
        tail = _log_tail(log) if log else ((proc.stdout or "") + "\n" + (proc.stderr or ""))[-2000:]
        raise GapseqError(f"command failed ({proc.returncode}): {' '.join(cmd)}\n{tail}")
    return proc


def version(dry_run: bool = False) -> str:
    if dry_run:
        return "[dry-run]"
    proc = run_gapseq(["-v"])
    return (proc.stdout or proc.stderr or "").strip().splitlines()[0]


def find(genome: str, output_dir: str, genome_id: str, threads: int = 1,
         dry_run: bool = False) -> bool:
    log = Path(output_dir) / "logs" / "find.log"
    run_gapseq(["find", "-p", "all", "-b", "200", "-t", "auto",
                "-K", str(threads), "-v", "1", str(genome)],
               log_path=log, dry_run=dry_run)
    return dry_run or (Path(output_dir) / f"{genome_id}-all-Pathways.tbl").exists()


def find_transport(genome: str, output_dir: str, genome_id: str,
                   threads: int = 1, dry_run: bool = False) -> bool:
    log = Path(output_dir) / "logs" / "transport.log"
    run_gapseq(["find-transport", "-b", "200", "-K", str(threads),
                "-v", "1", str(genome)], log_path=log, dry_run=dry_run)
    return dry_run or (Path(output_dir) / f"{genome_id}-Transporter.tbl").exists()


def draft(genome: str, output_dir: str, genome_id: str,
          high_evi_bitscore: int = 200, min_bitscore_core: int = 100,
          dry_run: bool = False) -> bool:
    out = Path(output_dir)
    reactions = out / f"{genome_id}-all-Reactions.tbl"
    transporters = out / f"{genome_id}-Transporter.tbl"
    pathways = out / f"{genome_id}-all-Pathways.tbl"
    log = out / "logs" / "draft.log"
    run_gapseq(["draft", "-r", str(reactions), "-t", str(transporters),
                "-c", str(genome), "-p", str(pathways),
                "-u", str(high_evi_bitscore), "-l", str(min_bitscore_core),
                "-b", "auto", "-n", genome_id, "-f", str(out) + "/"],
               log_path=log, dry_run=dry_run)
    return dry_run or (out / f"{genome_id}-draft.RDS").exists()


def medium(draft_rds: str, pathways: str, output_csv: str,
           dry_run: bool = False) -> bool:
    log = Path(output_csv).parent / "logs" / "medium.log"
    out = Path(output_csv)
    run_gapseq(["medium", "-m", str(draft_rds), "-p", str(pathways),
                "-f", str(out.parent) + "/", "-o", out.name],
               log_path=log, dry_run=dry_run)
    return dry_run or Path(output_csv).exists()


def _filter_candidates(weights_in: Path, genes_in: Path, exclude_file: Path,
                       weights_out: Path, genes_out: Path, log_path: Path | None = None,
                       dry_run: bool = False) -> None:
    ex = Path(exclude_file)
    if not ex.exists() or ex.stat().st_size <= len("seed\n"):
        shutil.copy2(weights_in, weights_out)
        shutil.copy2(genes_in, genes_out)
        return
    helper = Path(__file__).with_name("gapfill_filter.R")
    rscript = config.rscript()
    if not rscript:
        raise GapseqError("Rscript not found for gapfill candidate filtering")
    cmd = [rscript, str(helper), str(weights_in), str(genes_in), str(ex),
           str(weights_out), str(genes_out)]
    if dry_run:
        print("[dry-run] " + " ".join(cmd))
        return
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise GapseqError("gapfill candidate filter failed:\n" + ((proc.stdout or "") + (proc.stderr or ""))[-2000:])


def _write_minimal_evidence(draft_xml: Path, gapfilled_xml: Path, out_json: Path) -> None:
    import json
    import re
    def reactions(path: Path) -> set[str]:
        text = path.read_text(errors="replace")
        ids = set(re.findall(r'<reaction\b[^>]*\bid="([^"]+)"', text))
        return {re.sub(r'^R_', '', re.sub(r'_c0$', '', i)) for i in ids if not i.startswith('EX_')}
    added = sorted(reactions(gapfilled_xml) - reactions(draft_xml))
    payload = [{"seed": s, "geneless": True, "spontaneous": False} for s in added]
    out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_gapfill_evidence(workspace: Path, genome_id: str, draft_rds: Path,
                            gapfilled_rds: Path, weights: Path, genes: Path,
                            dry_run: bool = False) -> None:
    out_json = workspace / f"{genome_id}-gapfill-added.json"
    draft_xml = workspace / f"{genome_id}-draft.xml"
    gapfilled_xml = workspace / f"{genome_id}-gapfilled.xml"
    helper = Path(__file__).with_name("gapfill_evidence.R")
    rscript = config.rscript()
    seed_db = Path(config.db_dir()) / "reaction" / "seed_reactions_corrected.tsv" if config.db_dir() else None
    if rscript and helper.exists() and draft_rds.exists() and gapfilled_rds.exists() and seed_db is not None and seed_db.exists():
        cmd = [rscript, str(helper), str(draft_rds), str(gapfilled_rds),
               str(weights), str(genes), str(seed_db), str(out_json)]
        if dry_run:
            print("[dry-run] " + " ".join(cmd))
            return
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode == 0 and out_json.exists():
            return
    if draft_xml.exists() and gapfilled_xml.exists() and not dry_run:
        _write_minimal_evidence(draft_xml, gapfilled_xml, out_json)
    elif not dry_run:
        out_json.write_text("[]", encoding="utf-8")


def fill_with_excludes(workspace: str, genome_id: str, draft_rds: str,
                       medium_file: str, exclude_file: str | None = None,
                       bcore: int = 100, dry_run: bool = False) -> bool:
    """Invoke official gapseq fill with an optional candidate exclusion list."""
    ws = Path(workspace)
    wt = ws / f"{genome_id}-rxnWeights.RDS"
    gt = ws / f"{genome_id}-rxnXgenes.RDS"
    if not dry_run and (not wt.exists() or not gt.exists()):
        raise GapseqError(f"gapfill weights missing: {wt} / {gt}")
    tmp = ws / ".tmp" / f"gapfill-{uuid.uuid4().hex[:8]}"
    if not dry_run:
        tmp.mkdir(parents=True, exist_ok=True)
        input_model = tmp / "gapfill_input.RDS"
        shutil.copy2(draft_rds, input_model)
        weights = wt
        genes = gt
        if exclude_file:
            weights = tmp / "filtered-rxnWeights.RDS"
            genes = tmp / "filtered-rxnXgenes.RDS"
            _filter_candidates(wt, gt, Path(exclude_file), weights, genes,
                               log_path=ws / "logs" / "gapfill-filter.log",
                               dry_run=dry_run)
        log = ws / "logs" / "gapfill.log"
        run_gapseq(["fill", "-m", str(input_model), "-n", str(medium_file),
                    "-c", str(weights), "-g", str(genes), "-b", str(bcore),
                    "-f", str(tmp) + "/", "-w", "TRUE"],
                   log_path=log, dry_run=dry_run)
        out_rds = tmp / "gapfill_input.RDS"
        out_xml = tmp / "gapfill_input.xml"
        if not out_rds.exists() or not out_xml.exists():
            # official gapseq may retain a different stem; choose the newest model
            rds = sorted(tmp.glob("*.RDS"), key=lambda p: p.stat().st_mtime)
            xml = sorted(tmp.glob("*.xml"), key=lambda p: p.stat().st_mtime)
            if not rds or not xml:
                raise GapseqError("official gapseq fill produced no RDS/XML output")
            out_rds, out_xml = rds[-1], xml[-1]
        shutil.copy2(out_rds, ws / f"{genome_id}-gapfilled.RDS")
        shutil.copy2(out_xml, ws / f"{genome_id}-gapfilled.xml")
        _write_gapfill_evidence(ws, genome_id, Path(draft_rds),
                                ws / f"{genome_id}-gapfilled.RDS",
                                weights, genes, dry_run=dry_run)
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        print("[dry-run] official gapseq fill with candidate filtering for", genome_id)
    return True