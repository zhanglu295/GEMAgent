"""Public orchestration API for gemagent-core."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from . import config, manifest
from .adapters import gapseq
from .skills import SKILLS

__version__ = "0.2.0"


def _as_path(value: str | os.PathLike[str]) -> Path:
    return Path(value).expanduser().resolve()


def _genome_id(genome: Path) -> str:
    name = genome.name
    for suffix in (".gz", ".fna", ".fa", ".fasta"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def _stage_input(genome: Path, workspace: Path) -> Path:
    input_dir = workspace / "_input"
    input_dir.mkdir(parents=True, exist_ok=True)
    target = input_dir / genome.name
    if not target.exists():
        try:
            target.symlink_to(genome)
        except OSError:
            shutil.copy2(genome, target)
    return input_dir


def _pipeline():
    # Import after configure() has populated GAPSEQ_BIN and data paths.
    from . import pipeline_legacy
    return pipeline_legacy


def _infer_genome_id(workspace: Path, explicit: str = "") -> str:
    if explicit:
        return explicit
    drafts = sorted(workspace.glob("*-draft.RDS"))
    if drafts:
        return drafts[0].name[:-len("-draft.RDS")]
    finals = sorted(workspace.glob("*-gapfilled-gprchecked.xml"))
    if finals:
        return finals[0].name[:-len("-gapfilled-gprchecked.xml")]
    return ""


def _common_configure(db_dir: str, gapseq_bin: str, rscript: str,
                      eval_root: str) -> dict[str, str]:
    return config.configure(db_dir_override=db_dir, gapseq_bin=gapseq_bin,
                            rscript_bin=rscript, eval_root=eval_root)


def _configure_exchange_map(value: str) -> None:
    if not value:
        return
    path = Path(value).expanduser().resolve()
    data_dir = path.parent if path.is_file() else path
    os.environ["GEMAGENT_US_REVIEW_DIR"] = str(data_dir)

def _require_review_data() -> None:
    if not config.db_dir() or not Path(config.db_dir()).is_dir():
        raise RuntimeError(
            "Review skills require GEMAGENT_DB_DIR (the GEMAgent-curated data "
            "directory containing biomass/, media/, mimedb/, metabolite/, reaction/)."
        )



def _run_prune(ws: Path, gid: str, claims_dir: str = "", microbe_id: str = "",
               org_name: str = "", gtdb_summary: str = "", mimedb_tsv: str = "",
               features_only: bool = False) -> dict:
    model = ws / f"{gid}-gapfilled-gprchecked-usreviewed.xml"
    if not model.exists():
        model = ws / f"{gid}-gapfilled-gprchecked.xml"
    if not model.exists():
        model = ws / f"{gid}-gapfilled.xml"
    if not model.exists():
        raise RuntimeError(f"prune input model missing for {gid}")
    claims = Path(claims_dir) / f"{microbe_id or gid}.tsv" if claims_dir else Path("")
    out_dir = ws / "us_review_prune"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_xml = ws / f"{gid}-gapfilled-gprchecked-usreviewed-pruned.xml"
    out_rds = ws / f"{gid}-gapfilled-gprchecked-usreviewed-pruned.RDS"
    script = Path(__file__).resolve().parent / "impl" / "agent" / "agent_us_prune.py"
    cmd = [sys.executable, str(script), "--workspace", str(ws), "--genome-id", gid,
           "--model", str(model), "--draft", str(ws / f"{gid}-draft.xml"),
           "--gapfilled", str(ws / f"{gid}-gapfilled.xml"),
           "--gprchecked", str(ws / f"{gid}-gapfilled-gprchecked.xml"),
           "--gapfill-added", str(ws / f"{gid}-gapfill-added.json"),
           "--out-dir", str(out_dir), "--out-xml", str(out_xml),
           "--out-rds", str(out_rds), "--rscript", config.rscript(),
           "--org-name", org_name or ""]
    optional = {
        "--claims": claims,
        "--gtdb-summary": Path(gtdb_summary) if gtdb_summary else None,
        "--mimedb-tsv": Path(mimedb_tsv) if mimedb_tsv else None,
        "--pathways": ws / f"{gid}-all-Pathways.tbl",
        "--reactions": ws / f"{gid}-all-Reactions.tbl",
        "--transporter": ws / f"{gid}-Transporter.tbl",
        "--gpr-decisions": ws / "gpr_check" / "decisions.tsv",
        "--kofam": ws / "kofam" / f"{gid}_kofam.tsv",
    }
    for flag, path in optional.items():
        if path is not None and Path(path).exists():
            cmd += [flag, str(path)]
    if features_only:
        cmd.append("--features-only")
    subprocess.run(cmd, check=True)
    return {"status": "features-only" if features_only else "completed",
            "model": str(out_xml), "rds": str(out_rds),
            "features": str(out_dir / "features.tsv"),
            "decisions": str(out_dir / "decisions.tsv")}

def build(genome: str, out: str, recipe: str = "full", *, org_name: str = "",
          taxonomy: str = "", medium: str = "", eval_medium: str = "",
          claims_dir: str = "", us_microbe_id: str = "", kofam_tsv: str = "",
          db_dir: str = "", gapseq_bin: str = "", rscript: str = "",
          eval_root: str = "", threads: int = 1, force: bool = False,
          from_draft: bool = False, biomass_review: bool = True,
          gapfill_review: bool = True, gpr_check: bool = True,
          us_review: bool = True, us_prune: bool = True, dry_run: bool = False,
          gpr_min_or: int = 2, gtdb_summary: str = "", mimedb_tsv: str = "",
          exchange_map: str = "") -> dict[str, Any]:
    if recipe not in {"full", "base-only"}:
        raise ValueError(f"unknown recipe: {recipe}")
    genome_path = _as_path(genome)
    ws = _as_path(out)
    if not genome_path.is_file():
        raise FileNotFoundError(f"genome not found: {genome_path}")
    cfg = _common_configure(db_dir, gapseq_bin, rscript, eval_root)
    _configure_exchange_map(exchange_map)
    if dry_run:
        if recipe == "base-only":
            gid = _genome_id(genome_path)
            if not from_draft:
                gapseq.find(str(genome_path), str(ws), gid, threads, dry_run=True)
                gapseq.find_transport(str(genome_path), str(ws), gid, threads, dry_run=True)
                gapseq.draft(str(genome_path), str(ws), gid, dry_run=True)
            gapseq.medium(str(ws / f"{gid}-draft.RDS"),
                          str(ws / f"{gid}-all-Pathways.tbl"),
                          str(ws / f"{gid}-medium.csv"), dry_run=True)
            gapseq.fill_with_excludes(str(ws), gid,
                                      str(ws / f"{gid}-draft.RDS"),
                                      str(ws / f"{gid}-medium.csv"),
                                      None, bcore=100, dry_run=True)
        else:
            print("[dry-run] full recipe: official gapseq find/find-transport/draft/medium/fill + biomass/gapfill/GPR/expand/prune")
            print("[dry-run] requires GAPSEQ_BIN, Rscript, GEMAGENT_DB_DIR for review evidence")
        return {"status": "dry-run", "recipe": recipe, "config": cfg}
    ws.mkdir(parents=True, exist_ok=True)
    input_dir = _stage_input(genome_path, ws)
    if recipe == "full":
        _require_review_data()
    P = _pipeline()
    if threads > 0:
        P.n_threads = threads
    gid = ""
    status = "failed"
    error = ""
    try:
        if recipe == "base-only":
            gid = _genome_id(genome_path)
            if not from_draft:
                if not P.step_find_pathways(str(genome_path), str(ws), gid):
                    raise RuntimeError("gapseq find failed")
                if not P.step_predict_transporters(str(genome_path), str(ws), gid):
                    raise RuntimeError("gapseq find-transport failed")
                if not P.step_build_draft_gapseq(str(ws), gid, str(genome_path), force=force):
                    raise RuntimeError("gapseq draft failed")
            medium_file = medium or str(ws / f"{gid}-medium.csv")
            if not medium and not P.step_predict_medium(str(ws), gid):
                raise RuntimeError("official gapseq medium failed")
            if not gapseq.fill_with_excludes(str(ws), gid,
                                             str(ws / f"{gid}-draft.RDS"),
                                             medium_file, None, bcore=100):
                raise RuntimeError("official gapseq fill failed")
        else:
            ok = P.process_genome(
                input_dir=str(input_dir), output_dir=str(ws),
                org_name=org_name or None, taxonomy_tree=taxonomy or None,
                reuse_evidence=False, force_rebuild=force,
                from_draft=from_draft,
                gapfill_medium=medium or None, eval_medium=eval_medium or None,
                kofam_tsv=kofam_tsv or None, gpr_min_or=gpr_min_or,
                gene_pos_table=None, gtdb_fallback=True, allow_network=True,
                biomass_review=biomass_review, gapfill_review=gapfill_review,
                gpr_check=gpr_check, us_review=(us_review and bool(claims_dir)),
                us_claims_dir=claims_dir or None,
                us_microbe_id=us_microbe_id or None,
            )
            if not ok:
                raise RuntimeError("full recipe failed")
            gid = _infer_genome_id(ws)
            if us_prune:
                _run_prune(ws, gid, claims_dir, us_microbe_id, org_name,
                           gtdb_summary, mimedb_tsv)
        status = "completed"
    except Exception as exc:
        error = str(exc)
        raise
    finally:
        outputs = sorted(str(p) for p in ws.glob("*") if p.is_file())
        manifest.write_manifest(
            ws, status=status, recipe=recipe, genome=str(genome_path),
            genome_id=gid, error=error,
            config={"gapseq": config.gapseq(), "rscript": config.rscript(),
                    "db_dir": config.db_dir(), "gapseq_dat_dir": config.gapseq_dat_dir()},
            outputs=outputs,
        )
    return {"status": status, "recipe": recipe, "out": str(ws), "genome_id": gid}


def run_skill(name: str, workspace: str, *, org_name: str = "",
              taxonomy: str = "", genome_id: str = "", medium: str = "",
              kofam_tsv: str = "", claims_dir: str = "", microbe_id: str = "",
              db_dir: str = "", eval_root: str = "", force: bool = False,
              gpr_min_or: int = 2, gtdb_summary: str = "", mimedb_tsv: str = "",
              features_only: bool = False,
              exchange_map: str = "", base: str = "") -> dict[str, Any]:
    if name not in SKILLS:
        raise ValueError(f"unknown skill: {name}; choose: {', '.join(SKILLS)}")
    _common_configure(db_dir, "", "", eval_root)
    _configure_exchange_map(exchange_map)
    _require_review_data()
    ws = _as_path(workspace)
    P = _pipeline()
    gid = _infer_genome_id(ws, genome_id)
    if not gid:
        raise RuntimeError("cannot infer genome_id from workspace")
    status = "completed"
    try:
        if name == "biomass":
            medium_file = medium or str(ws / f"{gid}-medium.csv")
            gram = P.step_read_biomass_template(str(ws), gid)
            if not P.step_biomass_check(str(ws), gid, gram, org_name=org_name,
                                        taxonomy_tree=taxonomy, force=force,
                                        medium_csv=medium_file,
                                        kofam_tsv=kofam_tsv or None):
                raise RuntimeError("biomass review failed")
        elif name == "gapfill-review":
            medium_file = medium or str(ws / f"{gid}-medium.csv")
            P.step_review_gapfill_added(str(ws), gid, org_name=org_name,
                                        taxonomy_tree=taxonomy,
                                        gapfill_medium=medium_file,
                                        eval_medium=None)
            if not (ws / f"{gid}-gapfilled.xml").exists():
                raise RuntimeError("gapfill review produced no gapfilled model")
        elif name == "gpr":
            if not P.step_gpr_check(str(ws), gid, org_name=org_name,
                                    taxonomy_tree=taxonomy, kofam_tsv=kofam_tsv or None,
                                    min_or=gpr_min_or):
                raise RuntimeError("GPR review failed")
        elif name in ("us-review", "us-review-expand"):
            if not claims_dir:
                status = "skipped"
            else:
                P.step_us_review(str(ws), gid, microbe_id or gid, claims_dir,
                                 positive_only=True, base=base or "gapseq")
        elif name == "us-review-prune":
            _run_prune(ws, gid, claims_dir, microbe_id, org_name,
                       gtdb_summary, mimedb_tsv, features_only=features_only)
        else:
            raise AssertionError(name)
    except Exception as exc:
        status = "failed"
        manifest.write_manifest(ws, status=status, skill=name, error=str(exc))
        raise
    manifest.write_manifest(ws, status=status, skill=name, genome_id=gid,
                            config={"db_dir": config.db_dir(),
                                    "eval_root": config.env("GEMAGENT_EVAL_ROOT")})
    return {"status": status, "skill": name, "workspace": str(ws), "genome_id": gid}


def run_recipe(name: str, genome: str, out: str, **kwargs: Any) -> dict[str, Any]:
    return build(genome, out, recipe=name, **kwargs)
