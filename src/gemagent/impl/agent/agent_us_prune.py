#!/usr/bin/env python3
"""Model-consistency pruning for GEMAgent uptake/secretion calls.

The skill audits capabilities already present in the expanded model and emits
KEEP, BLOCK, or UNCERTAIN. It never consumes BacDive/NJC19 truth, and the
feature table never contains an observed label.

Evidence comes only from the GEMAgent model lineage (draft, gapfilled,
GPR-checked, expanded), gapfill provenance, GPR/KOfam, gapseq pathway tables,
network structure, explicit literature claims, GTDB-Tk lineage, and MIMEdb.
GTDB/MIMEdb are context and priors, not negative truth.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import cobra
from cobra import Configuration
from cobra.flux_analysis import flux_variability_analysis

HERE = Path(__file__).resolve()
AGENT_DIR = HERE.parent
IMPL_DIR = AGENT_DIR.parent
sys.path.insert(0, str(IMPL_DIR))
sys.path.insert(0, str(AGENT_DIR))
import paths
import taxonomy

Configuration().solver = "glpk"
TOL = 1e-7
SKILL_VERSION = "0.2.0"
norm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())
canon_gene = lambda s: re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")
EXPLICIT_NEGATIVE = re.compile(
    r"(?:\bnot\b|\bdid not\b|\bdoes not\b|\bcannot\b|\bunable to\b|\bno\b)"
    r"\s+(?:utiliz|utilis|assimilat|grow|produc|secret|excret|ferment|consum|catabol|detect)"
    r"|\bno growth\b|\bnegative for\b|\bnot detected\b|\bno production\b|"
    r"\babsent\b|\bincapable\b", re.I)
KEEP_DECISION = {"KEEP", "BLOCK", "UNCERTAIN"}
# Compounds whose panel name and model abbreviation differ; kept in sync with the
# expand skill's synonym table so a positive claim protects the exchange it created.
COMPOUND_SYNONYMS = {
    "cellulose": "cellul",
    "lichenan": "lichn",
    "glucomannan": "glcmannan",
    "xylan": "xylan",
    "inulin": "inulin",
    "pectin": "pect",
    "citruspectin": "pect",
    "tryptophan": "trp_L",
}


def sha256(path):
    if not path or not Path(path).is_file():
        return ""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(out_dir, payload):
    data = {"generated_at": datetime.now(timezone.utc).isoformat(), **payload}
    target = Path(out_dir) / "run_manifest.json"
    target.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
    return target


def first_existing(*values):
    for value in values:
        if value and Path(value).exists():
            return Path(value)
    return Path("")


def first_glob(workspace, pattern):
    hits = sorted(Path(workspace).glob(pattern))
    return hits[0] if hits else Path("")


def load_model(path):
    model = cobra.io.read_sbml_model(str(path))
    for met in model.metabolites:
        if met.formula and ";" in met.formula:
            met.formula = met.formula.split(";")[0].strip()
    return model


def norm_rxn(x):
    return re.sub(r"^R_", "", re.sub(r"_c0$", "", str(x)))


def xml_reactions(path):
    if not path or not Path(path).is_file():
        return set()
    text = Path(path).read_text(errors="replace")
    return {norm_rxn(x) for x in re.findall(r'<reaction\b[^>]*\bid="([^"]+)"', text)}


def exchange_reactions(model):
    """Real exchange reactions only.

    Expand may add single-metabolite helper reactions (GEMPATCH_*) that also look
    like boundary reactions.  They are part of a transport route, not a
    metabolite exchange, and must never be blocked as if they were secretion.
    """
    out = {}
    for rxn in model.reactions:
        rid = str(rxn.id)
        if not (rid.startswith("EX_") or rid.startswith("R_EX_")):
            continue
        if len(rxn.metabolites) != 1:
            continue
        met = next(iter(rxn.metabolites))
        comp = str(met.compartment).lower()
        # Foreign models use native compartment prefixes (CarveMe: C_e) while
        # AGORA/gapseq use e/e0.  Treat any explicit extracellular suffix alike.
        if not (comp in {"e", "e0"} or comp.endswith("_e") or comp.endswith("_e0")):
            continue
        out[rxn.id] = met
    return out


def support_reactions(model, exchange_rxn, met):
    return [r for r in model.reactions
            if r.id != exchange_rxn.id and not r.boundary and met in r.metabolites]


def read_added(path):
    if not path or not Path(path).is_file():
        return {}
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}
    return {norm_rxn(x.get("seed", "")): x for x in data if isinstance(x, dict)}


def read_gtdb(path, genome_id, workspace):
    if path and Path(path).is_file():
        p = Path(path)
        lines = [x.strip() for x in p.read_text(encoding="utf-8", errors="replace").splitlines()
                 if x.strip() and not x.lstrip().startswith("#")]
        for line in lines:
            compact = line.replace("|", ";")
            if "d__" in compact and "\t" not in line:
                return compact[compact.index("d__"):]
        try:
            df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
            cols = {str(c).lower().strip(): c for c in df.columns}
            taxcol = (cols.get("gtdb_taxonomy") or cols.get("gtdb_lineage")
                   or cols.get("taxonomy"))
            idcol = cols.get("user_genome") or cols.get("genome") or cols.get("genome_id")
            if taxcol:
                row = df.iloc[0]
                if idcol:
                    hit = df[df[idcol].astype(str).str.contains(genome_id, regex=False, na=False)]
                    if len(hit):
                        row = hit.iloc[0]
                return str(row[taxcol]).replace("|", ";")
            if len(df) == 1:
                values = [str(x) for x in df.iloc[0].tolist() if "d__" in str(x)]
                if values:
                    return values[0].replace("|", ";")
        except Exception:  # noqa: BLE001
            pass
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 2 and "d__" in parts[-1]:
                if genome_id in parts[0] or len(parts) == 2:
                    return parts[-1].replace("|", ";")
    return taxonomy.lineage_of(str(workspace)) if workspace else ""


def read_mimedb(path, org_name):
    if not path or not Path(path).is_file():
        return {}
    p = Path(path)
    sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    try:
        df = pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)
    except Exception:  # noqa: BLE001
        return {}
    if not len(df):
        return {}
    cols = {norm(c): c for c in df.columns}
    row = None
    q = norm(org_name)
    if q:
        for wanted in ("name", "species", "organism", "organismname"):
            col = cols.get(wanted)
            if col:
                hit = df[df[col].astype(str).map(norm).str.contains(q, regex=False, na=False)]
                if len(hit):
                    row = hit.iloc[0]
                    break
    if row is None and len(df) == 1:
        row = df.iloc[0]
    if row is None:
        return {}
    aliases = {
        "gram": ("gram", "gramstain", "gramstaining", "gramreaction"),
        "oxygen_requirement": ("oxygenrequirement", "oxygen", "oxygenrequirements"),
        "metabolism": ("metabolism", "metabolictype"),
        "energy_source": ("energysource", "energysources"),
        "habitat": ("habitat", "environment", "isolation"),
        "optimal_temperature": ("optimaltemperature", "temperature", "temperaturerange", "temp"),
        "description": ("description", "phenotype", "notes"),
    }
    out = {}
    for key, names in aliases.items():
        value = "unknown"
        for name in names:
            col = cols.get(name)
            if col and str(row.get(col, "")).strip():
                candidate = str(row.get(col)).strip()
                if candidate.lower() not in {"null", "nan", "none", "na"}:
                    value = candidate
                    break
        out[key] = value
    return out


def find_mapping_file(workspace):
    root = os.environ.get("GEMAGENT_EVAL_ROOT", "")
    candidates = []
    if root:
        candidates += [Path(root) / "us_review" / "data" / "exchange_map.tsv",
                       Path(root) / "data" / "exchange_map.tsv"]
    candidates += [Path(workspace) / "exchange_map.tsv"]
    return first_existing(*candidates)


def read_name_map(workspace, model):
    names = {}
    model_ids = {reaction.id for reaction in model.reactions}
    mapping = find_mapping_file(workspace)
    if mapping and mapping.is_file():
        df = pd.read_csv(mapping, sep="\t", dtype=str, keep_default_na=False)
        for _, r in df.iterrows():
            g = str(r.get("gapseq_exchange", "") or "").strip()
            abbr = str(r.get("abbr", "") or "").strip()
            local = "EX_%s_e0" % abbr if abbr else ""
            targets = set()
            if g:
                targets.add(g)
            if local and local in model_ids:
                targets.add(local)
            if not targets:
                continue
            for key in ("abbr", "compound_name", "exchange"):
                value = r.get(key)
                if isinstance(value, str) and value.strip():
                    names.setdefault(norm(value), set()).update(targets)
    for name, abbr in COMPOUND_SYNONYMS.items():
        local = "EX_%s_e0" % abbr
        if local in model_ids:
            names.setdefault(name, set()).add(local)
    for rxn in model.reactions:
        names.setdefault(norm(rxn.name), set()).add(rxn.id)
        short = re.sub(r"\s*[-_ ]?e0\s*$", "", str(rxn.name), flags=re.I)
        names.setdefault(norm(short), set()).add(rxn.id)
    return names


def read_claims(path, names, exchange_ids):
    evidence = {}
    if not path or not Path(path).is_file():
        return evidence
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for _, c in df.iterrows():
        compound = str(c.get("compound", "") or "")
        polarity = str(c.get("polarity", "") or "").lower()
        action = str(c.get("action", "") or "").lower()
        quote = str(c.get("quote", "") or "")
        if polarity not in ("positive", "negative"):
            continue
        targets = set()
        for key in ("gapseq_exchange", "exchange_id", "exchange"):
            value = str(c.get(key, "") or "").strip()
            if value in exchange_ids:
                targets.add(value)
        for hit in names.get(norm(compound), ()):
            if hit in exchange_ids:
                targets.add(hit)
        if not targets:
            key = norm(compound)
            if len(key) >= 5:
                for name, hits in names.items():
                    if len(name) >= 5 and (key.startswith(name) or name.startswith(key)):
                        targets.update(hit for hit in hits if hit in exchange_ids)
        directions = ["uptake", "secretion"] if action == "ambiguous" else [action]
        for direction in directions:
            if direction not in ("uptake", "secretion"):
                continue
            for exid in targets:
                evidence.setdefault((exid, direction), []).append(
                    {"polarity": polarity, "quote": quote,
                     "doi": str(c.get("doi", "") or ""), "compound": compound})
    return evidence


def read_table(path):
    if not path or not Path(path).is_file():
        return pd.DataFrame()
    lines = [x for x in Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
             if x.strip() and not x.lstrip().startswith("#")]
    if not lines:
        return pd.DataFrame()
    try:
        return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str,
                           keep_default_na=False)
    except Exception:  # noqa: BLE001
        return pd.DataFrame()


def read_pathways(path):
    df = read_table(path)
    info = {}
    reaction_to_pathways = {}
    if df.empty:
        return info, reaction_to_pathways
    cols = {norm(c): c for c in df.columns}
    idcol = cols.get("id") or df.columns[0]
    namecol = cols.get("name")
    ccol = cols.get("completeness")
    rfcol = cols.get("reactionsfound")
    kfcol = cols.get("keyreactionsfound")
    for _, r in df.iterrows():
        pid = str(r.get(idcol, "") or "").strip().strip("|")
        if not pid:
            continue
        try:
            completeness = float(str(r.get(ccol, "0") or "0"))
        except ValueError:
            completeness = 0.0
        info[pid] = {"pathway_id": pid,
                     "pathway_name": str(r.get(namecol, "") or "") if namecol else "",
                     "completeness": completeness}
        for col in (rfcol, kfcol):
            if not col:
                continue
            for token in re.split(r"[\s,]+", str(r.get(col, "") or "")):
                token = token.strip("|")
                if token:
                    reaction_to_pathways.setdefault(norm_rxn(token), set()).add(pid)
    return info, reaction_to_pathways


def read_reaction_annotations(path):
    df = read_table(path)
    out = {}
    if df.empty:
        return out
    for _, r in df.iterrows():
        rid = norm_rxn(r.get("rxn", ""))
        if not rid:
            continue
        pws = [x.strip("|") for x in re.split(r"[\s|,]+", str(r.get("pathway", "") or ""))
               if x.strip("|")]
        out[rid] = {"pathways": pws,
                    "status": str(r.get("status", "") or ""),
                    "pathway_status": str(r.get("pathway.status", "") or ""),
                    "pident": float(r.get("pident", 0) or 0),
                    "bitscore": float(r.get("bitscore", 0) or 0),
                    "qcovs": float(r.get("qcovs", 0) or 0)}
    return out


def read_transporter_annotations(path):
    df = read_table(path)
    out = {}
    if df.empty:
        return out
    for _, r in df.iterrows():
        try:
            pident = float(r.get("pident", 0) or 0)
            bitscore = float(r.get("bitscore", 0) or 0)
            qcovs = float(r.get("qcovs", 0) or 0)
        except ValueError:
            continue
        for rid in str(r.get("rea", "") or "").split(","):
            key = norm_rxn(rid.strip())
            if not key:
                continue
            item = out.setdefault(key, {"pident": 0.0, "bitscore": 0.0, "qcovs": 0.0})
            item["pident"] = max(item["pident"], pident)
            item["bitscore"] = max(item["bitscore"], bitscore)
            item["qcovs"] = max(item["qcovs"], qcovs)
    return out


def read_kofam(path):
    if not path or not Path(path).is_file():
        return {}
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    for i, line in enumerate(lines):
        if "gene name" in line and "\tKO" in line:
            start = i
            break
    if start is None:
        return {}
    header = lines[start].lstrip("#").lstrip("\t").split("\t")
    body = "\n".join(["\t".join(header)] + lines[start + 1:])
    try:
        df = pd.read_csv(io.StringIO(body), sep="\t", dtype=str, keep_default_na=False)
    except Exception:  # noqa: BLE001
        return {}
    out = {}
    for _, r in df.iterrows():
        gene = canon_gene(r.get("gene name", ""))
        if not gene:
            continue
        try:
            score = float(r.get("score", 0) or 0)
        except ValueError:
            score = 0.0
        out[gene] = {"ko": str(r.get("KO", "") or ""), "score": score}
    return out


def read_gpr(path):
    by_rxn = {}
    if not path or not Path(path).is_file():
        return by_rxn
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for _, r in df.iterrows():
        rid = norm_rxn(r.get("reaction_id", ""))
        if not rid:
            continue
        verdict = str(r.get("verdict", "") or "").upper()
        item = by_rxn.setdefault(rid, {"KEEP": 0, "REMOVE": 0, "UNCERTAIN": 0,
                                        "genes": []})
        if verdict in item:
            item[verdict] += 1
        if r.get("gene"):
            item["genes"].append(str(r.get("gene")))
    return by_rxn


def fva_exchanges(model):
    ex = exchange_reactions(model)
    for rid in ex:
        model.reactions.get_by_id(rid).bounds = (-1000.0, 1000.0)
    fva = flux_variability_analysis(
        model, reaction_list=list(ex), fraction_of_optimum=0.0, processes=1
    )
    return ex, fva


def signal_categories(signals):
    categories = set()
    for signal in signals:
        if signal in {"gapfill_new_exchange", "gapfill_added_support", "expand_new_exchange"}:
            categories.add("lineage")
        elif signal in {"no_gpr_support", "gpr_all_removed", "weak_homology", "no_homology"}:
            categories.add("gpr")
        elif signal in {"pathway_missing", "incomplete_pathway"}:
            categories.add("pathway")
        elif signal in {"empty_exchange", "no_direction_support", "no_production_path"}:
            categories.add("network")
        elif signal == "explicit_negative_claim":
            categories.add("literature_negative")
    return categories


def build_features(args):
    model = load_model(args.model)
    pathway_annotated = bool(args.pathways and Path(args.pathways).is_file())
    homology_annotated = bool(
        (args.reactions and Path(args.reactions).is_file()) or
        (args.transporter and Path(args.transporter).is_file()))
    gpr_table = bool(args.gpr_decisions and Path(args.gpr_decisions).is_file())
    model_annotated = pathway_annotated or homology_annotated
    ex, ex_fva = fva_exchanges(model)
    generational = {
        "draft": xml_reactions(args.draft),
        "gapfill": xml_reactions(args.gapfilled),
        "gpr": xml_reactions(args.gprchecked),
        "expand": xml_reactions(args.model),
    }
    added = read_added(args.gapfill_added)
    pathways, _ = read_pathways(args.pathways)
    reaction_ann = read_reaction_annotations(args.reactions)
    transporter_ann = read_transporter_annotations(args.transporter)
    kofam = read_kofam(args.kofam)
    gpr = read_gpr(args.gpr_decisions)
    names = read_name_map(args.workspace, model)
    claims = read_claims(args.claims, names, set(ex))
    gtdb = read_gtdb(args.gtdb_summary, args.genome_id, args.workspace)
    organism = args.org_name or taxonomy.org_name_of(args.workspace)
    mimedb_path = args.mimedb_tsv or str(paths.db("mimedb", "mimedb_microbes_v2.csv"))
    mimedb = read_mimedb(mimedb_path, organism)

    prelim = []
    analysis_ids = set()
    for exid, met in ex.items():
        vals = ex_fva.loc[exid]
        directions = []
        if float(vals["minimum"]) < -TOL:
            directions.append("uptake")
        if float(vals["maximum"]) > TOL:
            directions.append("secretion")
        supports = support_reactions(model, model.reactions.get_by_id(exid), met)
        cpd_id = met.id[:-3] + "_c0" if met.id.endswith("_e0") else ""
        internal = []
        if cpd_id and cpd_id in model.metabolites:
            internal = [r for r in model.reactions
                        if r.id != exid and not r.boundary and cpd_id in r.metabolites]
        for direction in directions:
            prelim.append({"exchange_id": exid, "direction": direction, "met": met,
                           "cpd_id": cpd_id, "flux_min": float(vals["minimum"]),
                           "flux_max": float(vals["maximum"]), "supports": supports,
                           "internal": internal})
            analysis_ids.update(r.id for r in supports)
            analysis_ids.update(r.id for r in internal)

    fva_lookup = {}
    if analysis_ids:
        fva2 = flux_variability_analysis(
            model, reaction_list=sorted(analysis_ids),
            fraction_of_optimum=0.0, processes=1,
        )
        for rid, row in fva2.iterrows():
            fva_lookup[norm_rxn(rid)] = (float(row["minimum"]), float(row["maximum"]))

    rows = []
    for item in prelim:
        exid = item["exchange_id"]
        direction = item["direction"]
        supports = item["supports"]
        internal = item["internal"]
        support_ids = [r.id for r in supports]
        internal_ids = [r.id for r in internal]
        support_norm = {norm_rxn(r.id) for r in supports}
        genes = set()
        for r in supports:
            genes.update(g.id for g in r.genes)
        gpr_keep = sum(gpr.get(x, {}).get("KEEP", 0) for x in support_norm)
        gpr_remove = sum(gpr.get(x, {}).get("REMOVE", 0) for x in support_norm)
        gpr_uncertain = sum(gpr.get(x, {}).get("UNCERTAIN", 0) for x in support_norm)
        gpr_all_removed = bool(gpr_remove and not gpr_keep and not gpr_uncertain)
        pidents, bitscores = [], []
        for rid in support_ids + internal_ids:
            key = norm_rxn(rid)
            if key in transporter_ann:
                pidents.append(transporter_ann[key]["pident"])
                bitscores.append(transporter_ann[key]["bitscore"])
            if key in reaction_ann:
                if reaction_ann[key]["pident"] > 0:
                    pidents.append(reaction_ann[key]["pident"])
                if reaction_ann[key]["bitscore"] > 0:
                    bitscores.append(reaction_ann[key]["bitscore"])
        max_pident = max(pidents) if pidents else 0.0
        max_bitscore = max(bitscores) if bitscores else 0.0
        weak_homology = bool((max_pident and max_pident < 30.0) or
                             (max_bitscore and max_bitscore < 50.0))
        no_homology = not pidents and not bitscores
        ko_hits = []
        for gene in genes:
            hit = kofam.get(canon_gene(gene))
            if hit:
                ko_hits.append(hit)
        ko_scores = [x["score"] for x in ko_hits if x["score"] > 0]

        pathway_ids = set()
        for rid in support_ids + internal_ids:
            pathway_ids.update(reaction_ann.get(norm_rxn(rid), {}).get("pathways", []))
        for rid in support_norm:
            meta = added.get(rid, {})
            value = meta.get("pathway", "") if isinstance(meta, dict) else ""
            if value:
                pathway_ids.update(x.strip("|") for x in str(value).split("|") if x.strip("|"))
        pathway_ids = {x for x in pathway_ids if x}
        pathway_names = [pathways[x]["pathway_name"] for x in pathway_ids if x in pathways]
        completeness = [pathways[x]["completeness"] for x in pathway_ids if x in pathways]
        max_completeness = max(completeness) if completeness else 0.0
        if not pathway_ids:
            pathway_status = "missing"
        elif max_completeness >= 80:
            pathway_status = "complete"
        elif max_completeness >= 50:
            pathway_status = "partial"
        elif max_completeness > 0:
            pathway_status = "incomplete"
        else:
            pathway_status = "incomplete_unknown"

        support_flux = any(
            abs(fva_lookup.get(norm_rxn(rid), (0.0, 0.0))[0]) > TOL or
            abs(fva_lookup.get(norm_rxn(rid), (0.0, 0.0))[1]) > TOL
            for rid in support_ids)
        production = []
        internal_use = []
        cpd = item["cpd_id"]
        cpd_met = model.metabolites.get_by_id(cpd) if cpd and cpd in model.metabolites else None
        for r in internal:
            if cpd_met is None:
                continue
            key = norm_rxn(r.id)
            coefficient = float(r.metabolites.get(cpd_met) or 0.0)
            lo, hi = fva_lookup.get(key, (0.0, 0.0))
            if direction == "secretion" and coefficient > 0 and hi > TOL:
                production.append(r.id)
            if direction == "uptake" and coefficient < 0 and lo < -TOL:
                internal_use.append(r.id)

        present = {stage: norm_rxn(exid) in ids for stage, ids in generational.items()}
        if present["draft"]:
            first_seen = "draft"
        elif present["gapfill"]:
            first_seen = "gapfill"
        elif present["gpr"]:
            first_seen = "gpr"
        elif present["expand"]:
            first_seen = "expand"
        else:
            first_seen = "unknown"
        support_generations = set()
        for rid in support_ids:
            key = norm_rxn(rid)
            if key in generational["draft"]:
                support_generations.add("draft")
            elif key in generational["gapfill"]:
                support_generations.add("gapfill")
            elif key in generational["gpr"]:
                support_generations.add("gpr")
            elif key in generational["expand"]:
                support_generations.add("expand")
        gapfill_new_exchange = (not present["draft"]) and present["gapfill"]
        expand_new_exchange = (not present["gpr"]) and present["expand"]
        gapfill_added_support = bool(support_norm & set(added)) or "gapfill" in support_generations

        claim_ev = claims.get((exid, direction), [])
        positive = [x for x in claim_ev if x["polarity"] == "positive"]
        negative = [x for x in claim_ev if x["polarity"] == "negative" and
                    EXPLICIT_NEGATIVE.search(x["quote"])]
        signals = []
        if gapfill_new_exchange:
            signals.append("gapfill_new_exchange")
        if gapfill_added_support:
            signals.append("gapfill_added_support")
        if expand_new_exchange:
            signals.append("expand_new_exchange")
        if not genes:
            signals.append("no_gpr_support")
        if gpr_all_removed:
            signals.append("gpr_all_removed")
        if weak_homology:
            signals.append("weak_homology")
        if no_homology and homology_annotated:
            signals.append("no_homology")
        if not pathway_ids and pathway_annotated:
            signals.append("pathway_missing")
        elif max_completeness and max_completeness < 50:
            signals.append("incomplete_pathway")
        if not supports:
            signals.append("empty_exchange")
        if not support_flux:
            signals.append("no_direction_support")
        if direction == "secretion" and not production:
            signals.append("no_production_path")
        if negative:
            signals.append("explicit_negative_claim")
        categories = signal_categories(signals)
        strong_homology = max_pident >= 30.0 or max_bitscore >= 50.0
        strong_gpr = bool(genes) and not gpr_all_removed and (
            gpr_keep > 0 or strong_homology or not gpr_table)
        pathway_good = bool(pathway_ids) and (max_completeness >= 50 or
                                               pathway_status == "incomplete_unknown")
        protected = bool(positive) or (
            strong_gpr and (pathway_good or strong_homology or not model_annotated))
        eligible = bool(negative) or len(categories - {"literature_negative"}) >= 2
        triage = "PROTECTED" if protected else ("BLOCK_ELIGIBLE" if eligible else "UNCERTAIN")

        rows.append({
            "exchange_id": exid,
            "direction": direction,
            "compound_name": getattr(model.reactions.get_by_id(exid), "name", ""),
            "flux_min": item["flux_min"],
            "flux_max": item["flux_max"],
            "present_in_draft": present["draft"],
            "present_after_gapfill": present["gapfill"],
            "present_after_gpr": present["gpr"],
            "present_in_expanded": present["expand"],
            "first_seen_stage": first_seen,
            "gapfill_new_exchange": gapfill_new_exchange,
            "expand_new_exchange": expand_new_exchange,
            "gapfill_added_support": gapfill_added_support,
            "support_reactions": "|".join(support_ids),
            "support_generations": "|".join(sorted(support_generations)),
            "gene_count": len(genes),
            "genes": "|".join(sorted(genes)),
            "gpr_support": bool(genes),
            "gpr_keep_count": gpr_keep,
            "gpr_remove_count": gpr_remove,
            "gpr_uncertain_count": gpr_uncertain,
            "gpr_all_removed": gpr_all_removed,
            "min_pident": min(pidents) if pidents else 0.0,
            "max_pident": max_pident,
            "max_bitscore": max_bitscore,
            "weak_homology": weak_homology,
            "no_homology": no_homology,
            "kofam_ko_count": len(ko_hits),
            "kofam_kos": "|".join(sorted({x["ko"] for x in ko_hits if x["ko"]})),
            "kofam_score_min": min(ko_scores) if ko_scores else 0.0,
            "kofam_score_max": max(ko_scores) if ko_scores else 0.0,
            "pathway_ids": "|".join(sorted(pathway_ids)),
            "pathway_names": "|".join(pathway_names),
            "pathway_completeness_max": max_completeness,
            "pathway_status": pathway_status,
            "pathway_present": bool(pathway_ids),
            "internal_reactions": "|".join(internal_ids),
            "production_reactions": "|".join(production),
            "precursor_use_reactions": "|".join(internal_use),
            "support_has_flux": support_flux,
            "network_consistency": ("ok" if support_flux and
                                    (direction == "uptake" or bool(production))
                                    else "inconsistent"),
            "empty_exchange": not supports,
            "positive_claim": bool(positive),
            "positive_claim_count": len(positive),
            "positive_claim_quotes": " || ".join(x["quote"][:220] for x in positive[:3]),
            "negative_claim": bool(negative),
            "negative_claim_quotes": " || ".join(x["quote"][:220] for x in negative[:3]),
            "gtdb_lineage": gtdb or "unknown",
            "mimedb_gram": mimedb.get("gram", "unknown"),
            "mimedb_oxygen": mimedb.get("oxygen_requirement", "unknown"),
            "mimedb_metabolism": mimedb.get("metabolism", "unknown"),
            "mimedb_energy_source": mimedb.get("energy_source", "unknown"),
            "mimedb_habitat": mimedb.get("habitat", "unknown"),
            "mimedb_temperature": mimedb.get("optimal_temperature", "unknown"),
            "mimedb_description": str(mimedb.get("description", "unknown"))[:600],
            "risk_signals": "|".join(signals),
            "risk_categories": "|".join(sorted(categories)),
            "triage": triage,
        })
    return model, pd.DataFrame(rows), gtdb, mimedb


def prompt(features, organism, gtdb, mimedb):
    eligible = features[features.triage == "BLOCK_ELIGIBLE"].copy()
    cols = ["exchange_id", "direction", "compound_name", "flux_min", "flux_max",
            "first_seen_stage", "gapfill_new_exchange", "gapfill_added_support",
            "support_reactions", "gene_count", "gpr_keep_count", "gpr_remove_count",
            "min_pident", "max_pident", "max_bitscore", "weak_homology",
            "pathway_ids", "pathway_completeness_max", "pathway_status",
            "internal_reactions", "production_reactions", "network_consistency",
            "positive_claim", "negative_claim_quotes", "risk_signals", "risk_categories"]
    table = eligible[[x for x in cols if x in eligible.columns]].to_csv(sep="\t", index=False)
    lines = [
        "You are auditing existing positive metabolite-exchange calls in a GEMAgent model.",
        "Use only the supplied GEMAgent model evidence, taxonomy context, and MIMEdb context.",
        "A positive literature claim is a protection set. Absence of literature is neutral.",
        "GTDB-Tk and MIMEdb are context/priors, not negative truth.",
        "Only an explicit negative statement (not utilized, no growth, not produced) is negative evidence.",
        "Return KEEP unless the evidence supports a broken or unsupported capability.",
        f"Organism: {organism}",
        f"GTDB-Tk lineage: {gtdb or 'unknown'}",
        "MIMEdb context: " + json.dumps(mimedb, ensure_ascii=False) if mimedb else "MIMEdb context: unknown",
        "",
        "Candidate evidence table (evaluation labels excluded):",
        table,
        "",
        'Return JSON {"decisions":[{"exchange_id":"...","direction":"uptake|secretion",',
        '"decision":"KEEP|BLOCK|UNCERTAIN","confidence":"high|medium|low","reason":"...",',
        '"evidence_refs":["field or value used"]}]}.',
        "Every BLOCK must include non-empty evidence_refs and must not contradict a positive claim.",
    ]
    return "\n".join(lines)


def ask_llm(features, organism, gtdb, mimedb):
    if features.empty or not (features.triage == "BLOCK_ELIGIBLE").any():
        return []
    try:
        from llm_config import llm
        resp = llm.bind(response_format={"type": "json_object"}).invoke([
            {"role": "system", "content": "You are a conservative GEM consistency auditor."},
            {"role": "user", "content": prompt(features, organism, gtdb, mimedb)},
        ])
        data = json.loads(resp.content)
        return [d for d in data.get("decisions", []) if isinstance(d, dict)]
    except Exception as exc:  # noqa: BLE001
        print(f"  [prune] LLM unavailable or invalid response: {exc}", file=sys.stderr)
        return []


def deterministic_decisions(features):
    rows = []
    for _, r in features.iterrows():
        if r["triage"] == "PROTECTED":
            decision = "KEEP"
            reason = "positive claim or strong model support"
        else:
            decision = "UNCERTAIN"
            reason = ("model risks are insufficient for pruning"
                      if r["triage"] == "UNCERTAIN"
                      else "eligible for audit; no LLM verdict")
        rows.append({"exchange_id": r["exchange_id"], "direction": r["direction"],
                     "decision": decision,
                     "confidence": "low" if decision == "UNCERTAIN" else "medium",
                     "reason": reason, "evidence_refs": [],
                     "triage": r["triage"], "source": "deterministic",
                     "applied": False})
    return rows


def merge_llm_decisions(base, llm_rows, features):
    eligible = {(r.exchange_id, r.direction) for _, r in
                features[features.triage == "BLOCK_ELIGIBLE"].iterrows()}
    protected = {(r.exchange_id, r.direction) for _, r in
                 features[features.triage == "PROTECTED"].iterrows()}
    feature_by_key = {(r.exchange_id, r.direction): r for _, r in features.iterrows()}
    by_key = {(x["exchange_id"], x["direction"]): x for x in base}
    n_weak_rejected = 0
    for row in llm_rows:
        exid = str(row.get("exchange_id", "") or "")
        direction = str(row.get("direction", "") or "").lower()
        key = (exid, direction)
        if key not in eligible:
            continue
        decision = str(row.get("decision", "") or "").upper()
        if decision not in KEEP_DECISION:
            continue
        refs = row.get("evidence_refs", [])
        if isinstance(refs, str):
            refs = [refs] if refs.strip() else []
        if decision == "BLOCK":
            if key in protected or not refs:
                continue
            feat = feature_by_key.get(key)
            if feat is None:
                continue
            def truthy(field):
                return str(feat.get(field, "")).strip().lower() in {"true", "1", "yes"}
            explicit_negative = truthy("negative_claim")
            hard_model_artifact = any(truthy(x) for x in (
                "empty_exchange", "gapfill_new_exchange", "expand_new_exchange",
                "gpr_all_removed"))
            # Missing GPR/pathway/production evidence is neutral, not negative.
            if not (explicit_negative or hard_model_artifact):
                n_weak_rejected += 1
                continue
        by_key[key].update({"decision": decision,
                            "confidence": str(row.get("confidence", "low") or "low"),
                            "reason": str(row.get("reason", "") or ""),
                            "evidence_refs": refs,
                            "source": "llm"})
    print(f"  [prune] rejected weak-evidence LLM blocks: {n_weak_rejected}")
    return list(by_key.values())


def capabilities(model, pairs):
    unique = sorted({exid for exid, _ in pairs})
    if not unique:
        return {}
    fva = flux_variability_analysis(
        model, reaction_list=unique, fraction_of_optimum=0.0, processes=1
    )
    out = {}
    for exid, direction in pairs:
        row = fva.loc[exid]
        out[(exid, direction)] = (
            float(row['minimum']) < -TOL if direction == 'uptake'
            else float(row['maximum']) > TOL
        )
    return out

def write_rds(xml_path, rds_path, model_rds, rscript):
    if not rds_path:
        return None
    rds_path = Path(rds_path)
    if not rscript:
        raise RuntimeError("RDS output requested but no Rscript executable was found")
    script = AGENT_DIR / "agent_gpr_check" / "write_rds.R"
    cmd = [rscript, str(script), "-x", str(xml_path), "-o", str(rds_path)]
    if model_rds and Path(model_rds).exists():
        cmd += ["-m", str(model_rds)]
    envs = dict(os.environ)
    genv = str(Path(sys.executable).resolve().parent.parent / "lib")
    if Path(genv).is_dir():
        envs["LD_LIBRARY_PATH"] = genv + (
            (":" + envs["LD_LIBRARY_PATH"]) if envs.get("LD_LIBRARY_PATH") else "")
    result = subprocess.run(cmd, env=envs)
    if result.returncode != 0 or not rds_path.exists() or rds_path.stat().st_size == 0:
        # Some community models carry non-stoichiometric formulas that sybil's
        # validator rejects.  The pruned SBML and the decisions are still valid,
        # so a failed convenience export must not fail the whole skill.
        print(f"WARNING: RDS export failed ({rds_path}); SBML and decisions are unaffected")
        return None
    return rds_path
def apply_decisions(args, features, decisions):
    model = load_model(args.model)
    base_growth = model.slim_optimize()
    positive = features[features.positive_claim == True]
    protected = [(r.exchange_id, r.direction) for _, r in positive.iterrows()]
    initially_capable = capabilities(model, protected)
    applied_rows = []
    for row in decisions:
        if str(row.get("decision", "")).upper() != "BLOCK":
            row["growth_before"] = base_growth
            row["growth_after"] = base_growth
            row["guard"] = ""
            applied_rows.append(row)
            continue
        exid = str(row.get("exchange_id", "") or "")
        direction = str(row.get("direction", "") or "").lower()
        if direction not in ("uptake", "secretion") or exid not in model.reactions:
            row.update({"applied": False, "guard": "invalid exchange/direction"})
            applied_rows.append(row)
            continue
        rxn = model.reactions.get_by_id(exid)
        old = tuple(rxn.bounds)
        if direction == "uptake":
            rxn.bounds = (max(0.0, old[0]), old[1])
        else:
            rxn.bounds = (old[0], min(0.0, old[1]))
        growth = model.slim_optimize()
        growth_ok = (math.isfinite(float(base_growth)) and base_growth > TOL and
                     math.isfinite(float(growth)) and growth >= 0.9 * base_growth)
        current_capable = capabilities(model, protected)
        positive_ok = all(current_capable[pair] == initial
                          for pair, initial in initially_capable.items())
        applied = bool(growth_ok and positive_ok)
        if not applied:
            rxn.bounds = old
            growth = model.slim_optimize()
        row.update({"applied": applied,
                    "bounds_before": f"{old[0]},{old[1]}",
                    "bounds_after": f"{rxn.bounds[0]},{rxn.bounds[1]}",
                    "growth_before": base_growth,
                    "growth_after": growth,
                    "guard": "" if applied else
                    ("growth_guard" if not growth_ok else "positive_claim_guard")})
        applied_rows.append(row)
    out_xml = Path(args.out_xml)
    out_xml.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(out_xml))
    return model, applied_rows, base_growth


def input_manifest(args):
    values = {
        "model": args.model, "draft": args.draft, "gapfilled": args.gapfilled,
        "gprchecked": args.gprchecked, "gapfill_added": args.gapfill_added,
        "pathways": args.pathways, "reactions": args.reactions,
        "transporter": args.transporter, "gpr_decisions": args.gpr_decisions,
        "kofam": args.kofam, "claims": args.claims,
        "gtdb_summary": args.gtdb_summary, "mimedb_tsv": args.mimedb_tsv,
    }
    return {key: {"path": str(value), "sha256": sha256(value)}
            for key, value in values.items() if value and Path(value).is_file()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workspace", required=True)
    ap.add_argument("--genome-id", required=True)
    ap.add_argument("--model", default="")
    ap.add_argument("--draft", default="")
    ap.add_argument("--gapfilled", default="")
    ap.add_argument("--gprchecked", default="")
    ap.add_argument("--gapfill-added", default="")
    ap.add_argument("--pathways", default="")
    ap.add_argument("--reactions", default="")
    ap.add_argument("--transporter", default="")
    ap.add_argument("--gpr-decisions", default="")
    ap.add_argument("--kofam", default="")
    ap.add_argument("--claims", default="")
    ap.add_argument("--gtdb-summary", default="")
    ap.add_argument("--mimedb-tsv", default="")
    ap.add_argument("--org-name", default="")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-xml", default="")
    ap.add_argument("--out-rds", default="")
    ap.add_argument("--rscript", default="")
    ap.add_argument("--features-only", action="store_true")
    args = ap.parse_args()

    ws = Path(args.workspace)
    gid = args.genome_id
    args.model = str(first_existing(args.model,
                                    ws / f"{gid}-gapfilled-gprchecked-usreviewed.xml",
                                    ws / f"{gid}-gapfilled-gprchecked.xml"))
    args.draft = str(first_existing(args.draft, ws / f"{gid}-draft.xml"))
    args.gapfilled = str(first_existing(args.gapfilled, ws / f"{gid}-gapfilled.xml"))
    args.gprchecked = str(first_existing(args.gprchecked,
                                         ws / f"{gid}-gapfilled-gprchecked.xml"))
    args.gapfill_added = str(first_existing(args.gapfill_added,
                                            ws / f"{gid}-gapfill-added.json"))
    args.pathways = str(first_existing(args.pathways,
                                       first_glob(ws, f"{gid}-all-Pathways.tbl")))
    args.reactions = str(first_existing(args.reactions,
                                        first_glob(ws, f"{gid}-all-Reactions.tbl")))
    args.transporter = str(first_existing(args.transporter,
                                          first_glob(ws, f"{gid}-Transporter.tbl")))
    args.gpr_decisions = str(first_existing(args.gpr_decisions,
                                            ws / "gpr_check" / "decisions.tsv"))
    args.kofam = str(first_existing(args.kofam,
                                    first_glob(ws, f"kofam/{gid}_kofam.tsv")))
    if not args.model or not Path(args.model).exists():
        raise RuntimeError("prune input model missing")
    if not args.out_xml:
        args.out_xml = str(ws /
                           f"{gid}-gapfilled-gprchecked-usreviewed-pruned.xml")
    if not args.out_rds:
        args.out_rds = str(Path(args.out_xml).with_suffix(".RDS"))
    args.rscript = args.rscript or os.environ.get("GEMAGENT_RSCRIPT", "") or "Rscript"

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    model, features, gtdb, mimedb = build_features(args)
    fpath = out_dir / "features.tsv"
    features.to_csv(fpath, sep="\t", index=False)
    print(f"FEATURES_WRITTEN {fpath} rows={len(features)}")
    counts = {k: int((features.triage == k).sum())
              for k in ("PROTECTED", "BLOCK_ELIGIBLE", "UNCERTAIN")}
    common = {"skill": "us-review-prune", "skill_version": SKILL_VERSION,
              "workspace": str(ws), "genome_id": gid, "features": str(fpath),
              "triage_counts": counts, "gtdb_lineage": gtdb or "unknown",
              "mimedb": mimedb, "inputs": input_manifest(args)}
    if args.features_only:
        write_manifest(out_dir, {"status": "features-only", **common})
        return

    decisions = deterministic_decisions(features)
    llm_rows = ask_llm(features, args.org_name, gtdb, mimedb)
    decisions = merge_llm_decisions(decisions, llm_rows, features)
    model, decisions, base_growth = apply_decisions(args, features, decisions)
    djson = out_dir / "decisions.json"
    dtsv = out_dir / "decisions.tsv"
    djson.write_text(json.dumps({"decisions": decisions,
                                 "counts": {k: sum(1 for x in decisions
                                                   if x.get("decision") == k)
                                            for k in sorted(KEEP_DECISION)}},
                                indent=2), encoding="utf-8")
    pd.DataFrame(decisions).to_csv(dtsv, sep="\t", index=False)
    rds_source = Path(args.gprchecked).with_suffix(".RDS")
    rds_written = write_rds(args.out_xml, args.out_rds, str(rds_source), args.rscript)
    write_manifest(out_dir, {
        "status": "completed", **common, "decisions_json": str(djson),
        "decisions_tsv": str(dtsv), "model_out": str(args.out_xml),
        "rds_out": str(rds_written or ""), "base_growth": base_growth,
        "counts": {k: sum(1 for x in decisions if x.get("decision") == k)
                   for k in sorted(KEEP_DECISION)},
        "applied_blocks": sum(1 for x in decisions
                              if x.get("decision") == "BLOCK" and x.get("applied")),
    })
    print(f"PRUNE_DONE decisions={dtsv} applied="
          f"{sum(bool(x.get('applied')) for x in decisions)}")


if __name__ == "__main__":
    main()
