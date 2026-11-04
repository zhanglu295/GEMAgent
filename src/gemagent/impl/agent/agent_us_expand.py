#!/usr/bin/env python3
"""us-review-expand: truth-free, namespace-agnostic uptake/secretion expansion.

Canonical implementation of the ``us-review-expand`` skill.  Positive
literature claims are resolved to native exchange reactions through the shared
metabolite crosswalk and applied as ADD (native transporter) or CREATE
(exchange + transport patch) decisions.

This deliberately does not read the BacDive evaluation panel.  It is intended
for AGORA2, gapseq, CarveMe and ModelSEED drafts; the caller supplies the base
tag and the exchange crosswalk.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
from scipy.optimize import linprog

try:  # installed package (pip install -e gemagent-core)
    from gemagent import namespace as nsm
except ImportError:  # running straight from a checkout
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from gemagent import namespace as nsm

DECODE = [("__40__", "("), ("__41__", ")"), ("__91__", "["), ("__93__", "]"),
          ("__45__", "-"), ("__95__", "_"), ("__46__", "."), ("__58__", ":")]

norm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def dec(s: str) -> str:
    for a, b in DECODE:
        s = s.replace(a, b)
    return s


def attrs(elem: str, name: str) -> str:
    m = re.search(r'\b%s="([^"]*)"' % re.escape(name), elem)
    return m.group(1) if m else ""


def parse_sbml(path: str):
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    params = {}
    for p in re.findall(r"<parameter\b[^>]*>", text):
        pid, value = attrs(p, "id"), attrs(p, "value")
        try:
            if pid and value != "":
                params[pid] = float(value)
        except ValueError:
            pass
    species = {}
    for block in re.findall(r"<species\b[^>]*>", text):
        sid = attrs(block, "id")
        if sid:
            species[sid] = {"name": attrs(block, "name"),
                            "compartment": attrs(block, "compartment")}
    reactions = []
    for block in re.findall(r"<reaction\b[^>]*>.*?</reaction>", text, re.S):
        rid = attrs(block, "id")
        lb = params.get(attrs(block, "fbc:lowerFluxBound"), -1000.0)
        ub = params.get(attrs(block, "fbc:upperFluxBound"), 1000.0)
        stoich = {}
        for section, sign in (("listOfReactants", -1.0), ("listOfProducts", 1.0)):
            found = re.search(rf"<{section}>(.*?)</{section}>", block, re.S)
            if not found:
                continue
            for sr in re.findall(r"<speciesReference\b[^>]*>", found.group(1)):
                sid = attrs(sr, "species")
                try:
                    coeff = float(attrs(sr, "stoichiometry") or "1")
                except ValueError:
                    continue
                stoich[sid] = stoich.get(sid, 0.0) + sign * coeff
        reactions.append({"id": rid, "name": attrs(block, "name"),
                          "lb": lb, "ub": ub, "stoich": stoich})
    objective = ""
    m = re.search(r'fbc:activeObjective="([^"]+)"', text)
    if m:
        obj = re.search(r'<fbc:objective\b[^>]*id="%s".*?</fbc:objective>' %
                        re.escape(m.group(1)), text, re.S)
        if obj:
            flux = []
            for fo in re.findall(r'<fbc:fluxObjective\b[^>]*/>', obj.group(0)):
                rid = attrs(fo, "fbc:reaction")
                try:
                    coefficient = float(attrs(fo, "fbc:coefficient") or "0")
                except ValueError:
                    coefficient = 0.0
                if rid:
                    flux.append((coefficient, rid))
            positive = [rid for coefficient, rid in flux if coefficient > 0]
            objective = positive[0] if positive else (flux[0][1] if flux else "")
    return species, reactions, objective


def exchange_records(species, reactions):
    out = {}
    for rxn in reactions:
        if len(rxn["stoich"]) != 1:
            continue
        sid = next(iter(rxn["stoich"]))
        decoded = dec(sid)
        comp = decoded.rsplit("_", 1)[-1] if "_" in decoded else ""
        if not (decoded.endswith("_e") or decoded.endswith("_e0") or
                decoded.endswith("(e)") or decoded.endswith("[e]")):
            continue
        rid = rxn["id"]
        canonical = rid[2:] if rid.startswith("R_") else rid
        canonical = dec(canonical)
        out[canonical] = {
            "canonical": canonical,
            "reaction_id": rid,
            "name": rxn["name"],
            "species": sid,
            "species_decoded": decoded,
            "compartment": comp,
            "reaction": rxn,
        }
    return out


def build_linear(species, reactions):
    ids = [r["id"] for r in reactions]
    species_ids = sorted({s for r in reactions for s in r["stoich"]})
    idx = {s: i for i, s in enumerate(species_ids)}
    S = np.zeros((len(species_ids), len(reactions)), dtype=float)
    for j, r in enumerate(reactions):
        for s, c in r["stoich"].items():
            S[idx[s], j] += c
    lb = np.array([r["lb"] for r in reactions], dtype=float)
    ub = np.array([r["ub"] for r in reactions], dtype=float)
    return ids, idx, S, lb, ub


def fva_one(S, lb, ub, j, growth_row=None, growth_min=None):
    n = S.shape[1]
    # Uptake/secretion calls are feasibility calls.  Do not condition them on
    # an artificial biomass optimum obtained after opening every exchange.
    # A growth constraint can still be supplied explicitly by another caller.
    if growth_row is None or growth_min is None:
        A_ub = None
        b_ub = None
    else:
        A_ub = -growth_row.reshape(1, -1)
        b_ub = np.array([-growth_min])
    c = np.zeros(n)
    c[j] = 1.0
    hi = linprog(-c, A_eq=S, b_eq=np.zeros(S.shape[0]), A_ub=A_ub, b_ub=b_ub,
                 bounds=list(zip(lb, ub)), method="highs")
    lo = linprog(c, A_eq=S, b_eq=np.zeros(S.shape[0]), A_ub=A_ub, b_ub=b_ub,
                 bounds=list(zip(lb, ub)), method="highs")
    return (float(lo.fun) if lo.success else float("nan"),
            float(-hi.fun) if hi.success else float("nan"))


def scan(path: str, selected=None, open_all=True):
    species, reactions, objective = parse_sbml(path)
    records = exchange_records(species, reactions)
    if objective not in {r["id"] for r in reactions}:
        raise RuntimeError(f"objective reaction not found: {objective!r}")
    # Preserve biomass under the model's own medium for QC/reporting.  This is
    # distinct from the intentionally permissive exchange-feasibility scan.
    ids0, _, S0, lb0, ub0 = build_linear(species, reactions)
    oi0 = ids0.index(objective)
    growth0 = np.zeros(len(reactions))
    growth0[oi0] = 1.0
    baseline = linprog(-growth0, A_eq=S0, b_eq=np.zeros(S0.shape[0]),
                       bounds=list(zip(lb0, ub0)), method="highs")
    baseline_growth = float(-baseline.fun) if baseline.success else float("nan")
    if open_all:
        for rec in records.values():
            rec["reaction"]["lb"] = -1000.0
    ids, _, S, lb, ub = build_linear(species, reactions)
    oi = ids.index(objective)
    growth = np.zeros(len(reactions))
    growth[oi] = 1.0
    targets = set(records) if selected is None else (set(records) & set(selected))
    calls = {}
    for canonical in sorted(targets):
        rec = records[canonical]
        j = ids.index(rec["reaction_id"])
        lo, hi = fva_one(S, lb, ub, j)
        calls[canonical] = {"minimum": lo, "maximum": hi,
                            "uptake": bool(lo < -1e-6),
                            "secretion": bool(hi > 1e-6)}
    return {"calls": calls, "growth": baseline_growth, "species": species,
            "reactions": reactions, "objective": objective, "records": records}


def load_crosswalk(path: str):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return df.to_dict("records")


def clean_panel_id(ex: str) -> str:
    ex = dec(ex)
    ex = re.sub(r"^EX_", "", ex)
    ex = re.sub(r"(?:\(e\)|\[e\]|_e0|_e)$", "", ex)
    return ex


def candidate_ids(row):
    values = []
    for key in ("gapseq_exchange", "exchange"):
        value = str(row.get(key, "") or "")
        if value:
            values.append(value)
            values.append(re.sub(r"^R_", "", value))
    for key in ("abbr", "cpd"):
        value = str(row.get(key, "") or "")
        if value:
            values.extend([f"EX_{value}_e", f"EX_{value}_e0",
                           f"R_EX_{value}_e", f"R_EX_{value}_e0",
                           f"EX_{value}(e)", f"R_EX_{value}(e)"])
    panel = clean_panel_id(str(row.get("exchange", "") or ""))
    if panel:
        values.extend([f"EX_{panel}_e", f"EX_{panel}_e0", f"EX_{panel}(e)"])
    answer = []
    for value in values:
        value = re.sub(r"^R_", "", value)
        if value not in answer:
            answer.append(value)
    return answer


def row_aliases(row):
    out = []
    for key in ("abbr", "compound_name"):
        value = str(row.get(key, "") or "")
        if value:
            out.append(value)
    out.append(clean_panel_id(str(row.get("exchange", "") or "")))
    return [norm(x) for x in out if norm(x)]


def aliases_dict(crosswalk):
    out = {}
    for row in crosswalk:
        for alias in row_aliases(row):
            out.setdefault(alias, []).append(row)
    return out


def resolve_rows(compound, crosswalk, alias_map):
    key = norm(compound)
    if not key:
        return []
    rows = list(alias_map.get(key, []))
    if rows:
        return rows
    hits = []
    if len(key) >= 4:
        for alias, vals in alias_map.items():
            if len(alias) >= 4 and (key.startswith(alias) or alias.startswith(key)):
                hits.extend(vals)
    unique = []
    seen = set()
    for row in hits:
        token = tuple((k, str(row.get(k, ""))) for k in
                      ("exchange", "abbr", "cpd", "gapseq_exchange"))
        if token not in seen:
            seen.add(token)
            unique.append(row)
    return unique


def resolve_exchange_ids(compound, records, crosswalk, alias_map):
    ids = []
    by_norm = {}
    for key in records:
        by_norm.setdefault(norm(key), key)
    rows = resolve_rows(compound, crosswalk, alias_map)
    for row in rows:
        for cand in candidate_ids(row):
            hit = cand if cand in records else by_norm.get(norm(cand), "")
            if hit and hit not in ids:
                ids.append(hit)
    return ids, rows


PHENO = re.compile(r"utiliz|utilis|assimilat|grow|growth|catabol|consum|ferment|produc|secret|excret|synthesi|detect|degrad", re.I)
RECIPE = re.compile(r"cultured in|grown in|growth medium|medium|broth|agar|supplemented with|yeast extract|peptone|buffer|w/v|v/v", re.I)
GENET = re.compile(r"mutant|knockout|knock-out|deletion|overexpress|expression|operon|plasmid|recombinant|heterologous|transform|gene cluster", re.I)


def evidence_gate(claim):
    form = str(claim.get("evidence_form", "") or "").lower()
    subject = str(claim.get("subject_level", "") or "").lower()
    attribution = str(claim.get("attribution", "") or "").lower()
    quote = str(claim.get("quote", "") or "")
    if not str(claim.get("doi", "") or "").strip() or not quote.strip():
        return False, "missing DOI or source quote"
    if form != "phenotype_test":
        return False, "evidence_form is not phenotype_test"
    ok_pair = ((subject in {"strain", "species"} and attribution in {"strain", "species"})
               or (subject == "genus" and attribution == "genus"))
    if not ok_pair:
        return False, "subject/attribution not strain- or species-consistent"
    if not PHENO.search(quote):
        return False, "quote lacks a phenotype verb"
    if RECIPE.search(quote) and not PHENO.search(quote):
        return False, "recipe-only evidence"
    if GENET.search(quote) and not PHENO.search(quote):
        return False, "genetic-context-only evidence"
    if claim.get("quote_compound_consistent") is False:
        return False, "quote-compound inconsistency"
    return True, "passed evidence gate"


def load_claims(path, aliases):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    text = (df.get("quote", pd.Series("", index=df.index)).astype(str) + " " +
            df.get("paper", pd.Series("", index=df.index)).astype(str))
    contradictory = {
        "mannitol": {"sorbitol"}, "sorbitol": {"mannitol"},
        "glucose": {"galactose"}, "galactose": {"glucose"},
        "fructose": {"mannose"}, "mannose": {"fructose"},
    }
    out = []
    for idx, row in df.iterrows():
        item = row.to_dict()
        ckey = norm(item.get("compound", ""))
        chem_key = ckey[1:] if len(ckey) > 2 and ckey[0] in {"d", "l"} else ckey
        qkey = norm(text.iloc[idx])
        conflict = [x for x in contradictory.get(chem_key, set()) if norm(x) in qkey]
        consistent = (not conflict) and (ckey in qkey or ckey in aliases)
        item["quote_compound_consistent"] = bool(consistent)
        item["consistency_note"] = ("contradictory compound in quote: " + ",".join(conflict)
                                     if conflict else "")
        out.append(item)
    return pd.DataFrame(out) if out else pd.DataFrame()


def model_metabolite_id(model, sid, ns):
    """Map a raw SBML species id onto the in-memory cobra id."""
    return nsm.model_metabolite_id(model, sid, ns)


def find_cytosolic(model, rec, ns):
    """Cytosolic partner of an exchange species, in the model's namespace."""
    return nsm.cytosolic_partner(model, rec["species"], ns)


def role_of_species(species, sid):
    decoded = dec(sid)
    if any(decoded.endswith(x) for x in ("_e", "_e0", "(e)", "[e]")):
        return "external"
    return "internal"


def add_native_transport(model, rec, label, ns):
    ext = model_metabolite_id(model, rec["species"], ns)
    if not ext:
        return None, "no external metabolite in model namespace"
    cyto = find_cytosolic(model, rec, ns)
    if not cyto:
        return None, "no cytosolic partner metabolite"
    rid = "ADDTRAN_" + hashlib.sha1((label + rec["canonical"]).encode()).hexdigest()[:12]
    if rid in model.reactions:
        return model.reactions.get_by_id(rid), "already present"
    rxn = cobra.Reaction(rid, name="GEMAgent literature-supported transport " + rec["name"])
    rxn.bounds = (-1000.0, 1000.0)
    rxn.add_metabolites({model.metabolites.get_by_id(ext): -1.0,
                         model.metabolites.get_by_id(cyto): 1.0})
    model.add_reactions([rxn])
    return rxn, "added"


def add_create(model, row, direction, label, ns):
    token = str(row.get("cpd", "") or clean_panel_id(str(row.get("exchange", ""))))
    token = re.sub(r"[^A-Za-z0-9_]", "_", token) or "claim"
    # the id suffix comes from the family, the compartment attribute from the
    # model's own vocabulary (e.g. CarveMe keeps C_e, SEED keeps e0)
    ext_comp = "e0" if ns.family == "seed" else "e"
    ext_attr = ns.compartment_attr("e")
    ext_id = f"{token}_{ext_comp}"
    if ext_id not in model.metabolites:
        model.add_metabolites([cobra.Metabolite(ext_id, name=token,
                                               compartment=ext_attr)])
    exid = "EX_" + token + "_" + ext_comp
    if exid not in model.reactions:
        ex = cobra.Reaction(exid, name="GEMAgent literature-created exchange " + token)
        ex.bounds = (-1000.0, 1000.0)
        ex.add_metabolites({model.metabolites.get_by_id(ext_id): -1.0})
        model.add_reactions([ex])
    patch_id = "GEMPATCH_" + hashlib.sha1((label + token + direction).encode()).hexdigest()[:12]
    if patch_id not in model.reactions:
        patch = cobra.Reaction(patch_id, name="GEMAgent literature transport patch")
        patch.bounds = (0.0, 1000.0)
        sign = -1.0 if direction == "uptake" else 1.0
        patch.add_metabolites({model.metabolites.get_by_id(ext_id): sign})
        model.add_reactions([patch])
    return exid, patch_id


def resolve_exchange_map(explicit: str) -> str:
    """Locate the panel -> native exchange crosswalk.

    Precedence: explicit ``--exchange-map``, then ``GEMAGENT_US_REVIEW_DIR``
    (the directory holding ``exchange_map.tsv``), then the evaluation tree
    under ``GEMAGENT_EVAL_ROOT``.
    """
    candidates: list[Path] = []
    if explicit:
        given = Path(explicit)
        candidates.append(given / "exchange_map.tsv" if given.is_dir() else given)
    env_dir = os.environ.get("GEMAGENT_US_REVIEW_DIR", "").strip()
    if env_dir:
        base = Path(env_dir)
        candidates += [base / "exchange_map.tsv", base / "data" / "exchange_map.tsv"]
    eval_root = os.environ.get("GEMAGENT_EVAL_ROOT", "").strip()
    if eval_root:
        candidates.append(Path(eval_root) / "metabolite_uptake_secretion"
                          / "us_review" / "data" / "exchange_map.tsv")
    for cand in candidates:
        if cand.is_file():
            return str(cand)
    raise SystemExit(
        "us-review-expand: exchange map not found; pass --exchange-map or set "
        "GEMAGENT_US_REVIEW_DIR")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--base", default="",
                    help="base name: gapseq/carveme/modelseed/agora2, or auto")
    fam = ap.add_mutually_exclusive_group()
    fam.add_argument("--gapseq", dest="base_flag", action="store_const", const="gapseq")
    fam.add_argument("--carveme", dest="base_flag", action="store_const", const="carveme")
    fam.add_argument("--modelseed", dest="base_flag", action="store_const", const="modelseed")
    fam.add_argument("--agora2", dest="base_flag", action="store_const", const="agora2")
    ap.add_argument("--claims", required=True)
    ap.add_argument("--exchange-map", default="",
                    help="panel crosswalk TSV; default from GEMAGENT_US_REVIEW_DIR")
    ap.add_argument("--out-xml", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--allow-create", action="store_true")
    args = ap.parse_args()

    requested = args.base_flag or args.base or "auto"
    if requested == "auto":
        ns = nsm.detect(args.model)
    else:
        try:
            ns = nsm.for_base(requested)
        except nsm.UnknownNamespace:
            # a custom label (e.g. "gapseq-ft") still gets the right adapter
            ns = nsm.detect(args.model)
    args.base = requested if requested != "auto" else ns.name

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exchange_map = resolve_exchange_map(args.exchange_map)
    crosswalk = load_crosswalk(exchange_map)
    alias_map = aliases_dict(crosswalk)
    species0, reactions0, _ = parse_sbml(args.model)
    records = exchange_records(species0, reactions0)
    claims = load_claims(args.claims, set(alias_map))
    claim_rows = []
    decisions = []
    selected = set()

    for _, claim in claims.iterrows():
        ids, rows = resolve_exchange_ids(claim.get("compound", ""), records,
                                         crosswalk, alias_map)
        selected.update(ids)
        direction = str(claim.get("action", "") or "").lower()
        polarity = str(claim.get("polarity", "") or "").lower()
        consistent = bool(claim.get("quote_compound_consistent", False))
        base_row = claim.to_dict()
        if ids:
            for exid in ids:
                item = dict(base_row)
                item["exchange_id"] = records[exid]["reaction_id"]
                claim_rows.append(item)
        else:
            item = dict(base_row)
            item["exchange_id"] = ""
            claim_rows.append(item)
        gate_ok, gate_note = evidence_gate(claim.to_dict())
        if not consistent or polarity not in ("positive", "negative"):
            continue
        if polarity != "positive" or direction not in ("uptake", "secretion") or not gate_ok:
            continue
        if ids:
            for exid in ids:
                decisions.append({"exchange_id": exid, "action": "ADD",
                                  "direction": direction,
                                  "compound": claim.get("compound", ""),
                                  "row": rows[0] if rows else {}, "claim": claim.to_dict(),
                                  "gate": gate_note})
        elif args.allow_create and rows:
            decisions.append({"exchange_id": "", "action": "CREATE",
                              "direction": direction,
                              "compound": claim.get("compound", ""),
                              "row": rows[0], "claim": claim.to_dict()})

    # Collapse duplicate literature sentences onto one edit per reaction/action.
    unique_decisions = {}
    for item in decisions:
        row = item.get("row", {})
        target = item.get("exchange_id", "") or f"{row.get('cpd','')}|{row.get('exchange','')}|{item.get('compound','')}"
        key = (target, item.get("action", ""), item.get("direction", ""))
        if key in unique_decisions:
            unique_decisions[key]["n_claim_sentences"] += 1
        else:
            item["n_claim_sentences"] = 1
            unique_decisions[key] = item
    decisions = list(unique_decisions.values())

    # Densify claims so prune can see the actual native exchange IDs.
    pd.DataFrame(claim_rows).to_csv(out_dir / "claims_resolved.tsv", sep="\t", index=False)

    # Re-scan the selected targets after opening the exchange bounds.
    checked = scan(args.model, selected=selected, open_all=True)
    pending_add = []
    for dec_item in decisions:
        if dec_item["action"] != "ADD":
            continue
        call = checked["calls"].get(dec_item["exchange_id"], {})
        if not call.get(dec_item["direction"], False):
            pending_add.append(dec_item)
        else:
            dec_item["action"] = "NO_ACTION"

    species, reactions, objective = parse_sbml(args.model)
    model = cobra.io.read_sbml_model(args.model)
    for met in model.metabolites:
        if met.formula and ";" in met.formula:
            met.formula = met.formula.split(";")[0].strip()

    added_native = 0
    created = 0
    for dec_item in pending_add:
        rec = checked["records"][dec_item["exchange_id"]]
        rxn, note = add_native_transport(model, rec, args.base + "_" + args.organism, ns)
        if rxn is not None:
            dec_item["action"] = "ADD"
            dec_item["status"] = note
            added_native += 1
        else:
            dec_item["action"] = "REFUSED"
            dec_item["status"] = note
    for dec_item in decisions:
        if dec_item["action"] != "CREATE":
            continue
        try:
            exid, patch = add_create(model, dec_item["row"], dec_item["direction"],
                                     args.base + "_" + args.organism, ns)
            dec_item["exchange_id"] = exid
            dec_item["status"] = "created exchange+patch " + patch
            created += 1
        except Exception as exc:  # noqa: BLE001
            dec_item["action"] = "REFUSED"
            dec_item["status"] = "create failed: " + str(exc)

    growth0 = checked["growth"]
    growth1 = model.slim_optimize()
    guard = ""
    if created and (not np.isfinite(growth1) or growth1 < 0.5 * growth0):
        model = cobra.io.read_sbml_model(args.model)
        for dec_item in decisions:
            if dec_item["action"] == "CREATE":
                dec_item["action"] = "ROLLED_BACK"
        growth1 = model.slim_optimize()
        guard = "create growth guard: <50% baseline"

    Path(args.out_xml).parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, args.out_xml)
    rows_out = []
    for d in decisions:
        c = d.get("claim", {})
        rows_out.append({"exchange_model": d.get("exchange_id", ""),
                         "exchange_panel": d.get("row", {}).get("exchange", ""),
                         "action": d.get("action", ""), "direction": d.get("direction", ""),
                         "compound": d.get("compound", ""), "doi": c.get("doi", ""),
                         "quote": str(c.get("quote", ""))[:240],
                         "evidence_form": c.get("evidence_form", ""),
                         "subject_level": c.get("subject_level", ""),
                         "attribution": c.get("attribution", ""),
                         "quote_compound_consistent": c.get("quote_compound_consistent", ""),
                         "n_claim_sentences": d.get("n_claim_sentences", 1),
                         "source": "rule", "status": d.get("status", "")})
    pd.DataFrame(rows_out).to_csv(out_dir / "decisions.tsv", sep="\t", index=False)
    summary = {"base": args.base, "namespace": ns.family, "organism": args.organism,
               "n_claims": int(len(claims)), "n_resolved_claims": int(sum(bool(x) for x in pd.DataFrame(claim_rows).get("exchange_id", []))),
               "n_decisions": len(decisions), "n_added_transport": added_native,
               "n_created_exchanges": created, "growth_before": growth0,
               "growth_after": growth1, "growth_guard": guard,
               "truth_used": False,
               "inputs": {"model": args.model, "claims": args.claims,
                          "exchange_map": exchange_map}}
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("EXPAND_DONE " + json.dumps(summary))


if __name__ == "__main__":
    main()
