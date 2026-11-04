#!/usr/bin/env python3
"""gapseq_scan.py - drive a gapseq SBML model (same protocol as the AGORA2 track).

Differences from agora_scan.py: gapseq species are M_cpdXXXXX_e0 / M_cpdXXXXX_c0, exchange
reactions are R_EX_cpdXXXXX_e0, and the FT medium is already baked into the fbc bounds.

usage: gapseq_scan.py <model.xml> [MicrobeID] [compare_column]
"""
import os
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog

US_DATA_DIR = os.environ.get("GEMAGENT_US_DATA_DIR", "").strip()
if not US_DATA_DIR:
    eval_root = os.environ.get("GEMAGENT_EVAL_ROOT", "").strip()
    if eval_root:
        US_DATA_DIR = str(Path(eval_root).expanduser() / "metabolite_uptake_secretion")
DECODE = [("__40__", "("), ("__41__", ")"), ("__91__", "["), ("__93__", "]"),
          ("__45__", "-"), ("__95__", "_")]

def dec(s):
    for a, b in DECODE:
        s = s.replace(a, b)
    return s

def attrs(elem, name):
    m = re.search(r'\b%s="([^"]*)"' % re.escape(name), elem)
    return m.group(1) if m else ""

def parse(path):
    t = open(path, encoding="utf-8", errors="ignore").read()
    params = {}
    for p in re.findall(r"<parameter\b[^>]*>", t):
        i, v = attrs(p, "id"), attrs(p, "value")
        if i and v:
            try:
                params[i] = float(v)
            except ValueError:
                pass
    spec = {}
    for s in re.findall(r"<species\b[^>]*>", t):
        i = attrs(s, "id")
        if i:
            spec[i] = attrs(s, "name")
    rxs = []
    for r in re.findall(r"<reaction\b[^>]*>.*?</reaction>", t, re.S):
        rid = attrs(r, "id")
        lb = params.get(attrs(r, "fbc:lowerFluxBound"), -1000.0)
        ub = params.get(attrs(r, "fbc:upperFluxBound"), 1000.0)
        stoich = {}
        for block, sign in ((re.search(r"<listOfReactants>(.*?)</listOfReactants>", r, re.S), -1.0),
                            (re.search(r"<listOfProducts>(.*?)</listOfProducts>", r, re.S), 1.0)):
            if not block:
                continue
            for sr in re.findall(r"<speciesReference\b[^>]*>", block.group(1)):
                sid = attrs(sr, "species")
                try:
                    stoich[sid] = stoich.get(sid, 0.0) + sign * float(attrs(sr, "stoichiometry") or "1")
                except ValueError:
                    pass
        rxs.append(dict(id=rid, lb=lb, ub=ub, stoich=stoich))
    m = re.search(r'fbc:activeObjective="([^"]+)"', t)
    obj = ""
    if m:
        blk = re.search(r'<fbc:objective\b[^>]*id="%s".*?</fbc:objective>' % re.escape(m.group(1)),
                        t, re.S)
        if blk:
            fo = re.search(r'fbc:reaction="([^"]+)"', blk.group(0))
            obj = fo.group(1) if fo else ""
    return spec, rxs, obj

def build(spec, rxs):
    sp = sorted({s for r in rxs for s in r["stoich"]})
    idx = {s: i for i, s in enumerate(sp)}
    S = np.zeros((len(sp), len(rxs)))
    for j, r in enumerate(rxs):
        for s, c in r["stoich"].items():
            S[idx[s], j] += c
    return sp, idx, S, np.array([r["lb"] for r in rxs], float), np.array([r["ub"] for r in rxs], float)

def fva(S, lb, ub, j, growth_row=None, growth_min=None):
    n = S.shape[1]
    A_ub = b_ub = None
    if growth_row is not None and growth_min is not None:
        A_ub = -growth_row.reshape(1, -1)
        b_ub = np.array([-growth_min])
    c = np.zeros(n); c[j] = 1.0
    hi = linprog(-c, A_eq=S, b_eq=np.zeros(S.shape[0]), A_ub=A_ub, b_ub=b_ub,
                 bounds=list(zip(lb, ub)), method="highs")
    lo = linprog(c, A_eq=S, b_eq=np.zeros(S.shape[0]), A_ub=A_ub, b_ub=b_ub,
                 bounds=list(zip(lb, ub)), method="highs")
    return (lo.fun if lo.success else np.nan, -hi.fun if hi.success else np.nan)

def is_exchange(r):
    """gapseq exchange = reaction whose only species is extracellular (_e0)."""
    refs = list(r["stoich"])
    if len(refs) != 1:
        return None
    sid = dec(refs[0])
    if sid.endswith("_e0") or sid.endswith("(e)") or sid.endswith("[e]"):
        rid = r["id"]
        return rid[2:] if rid.startswith("R_") else rid
    return None

def scan(model_path, growth_frac=0.01, open_all=True):
    spec, rxs, obj = parse(model_path)
    if open_all:
        # The AGORA2 track asks a pure network question: can this model move flux through
        # the exchange at all?  The gapseq export bakes the gapfilling medium into the
        # bounds, which would otherwise fold "is it in the medium" into the answer and
        # make the two tracks incomparable.  Opening every exchange restores the same
        # semantics on both sides (the ADD arm then means "add a route", not "add a
        # nutrient").
        n_open = 0
        for r in rxs:
            if is_exchange(r) and r["lb"] > -1000.0:
                r["lb"] = -1000.0
                n_open += 1
        print("opened %d exchange bounds for the uniform protocol" % n_open)
    ids = [r["id"] for r in rxs]
    sp, idx, S, lb, ub = build(spec, rxs)
    oi = ids.index(obj) if obj in ids else None
    if oi is None:
        raise SystemExit("objective reaction %r not found" % obj)
    grow = np.zeros(len(rxs)); grow[oi] = 1.0
    g = linprog(-grow, A_eq=S, b_eq=np.zeros(S.shape[0]),
                bounds=list(zip(lb, ub)), method="highs")
    gmax = -g.fun if g.success else float("nan")
    gmin = 0.0 if not np.isfinite(gmax) else max(0.0, growth_frac * gmax)
    out = {}
    for j, r in enumerate(rxs):
        name = is_exchange(r)
        if not name:
            continue
        lo, hi = fva(S, lb, ub, j, growth_row=grow, growth_min=gmin)
        out[name] = dict(j=j, lo=lo, hi=hi, lb=lb[j], ub=ub[j],
                         can_uptake=bool(lo < -1e-6), can_secrete=bool(hi > 1e-6))
    return out, gmax, spec, rxs, obj

def main():
    model = sys.argv[1]
    org = sys.argv[2] if len(sys.argv) > 2 else ""
    col = sys.argv[3] if len(sys.argv) > 3 else "gapseq"
    ex, gmax, spec, rxs, obj = scan(model)
    print("model: %s | exchanges %d | objective %s | max biomass %.4f"
          % (model.split("/")[-1], len(ex), obj, gmax))
    npnl = sum(1 for v in ex.values() if v["lb"] < -1e-9)
    print("exchanges open for uptake (lb<0): %d of %d" % (npnl, len(ex)))
    if org:
        if not US_DATA_DIR:
            raise SystemExit(
                "Set GEMAGENT_US_DATA_DIR or GEMAGENT_EVAL_ROOT before comparing "
                "a model with the published panel."
            )
        data_dir = Path(US_DATA_DIR)
        panel = pd.read_csv(data_dir / "agora2_ref" / "agora2_S4c_panels.tsv", sep="\t")
        mp = pd.read_csv(data_dir / "us_review" / "data" / "exchange_map.tsv", sep="\t")
        p2g = {}
        for _, r in mp.iterrows():
            g = r.gapseq_exchange if isinstance(r.gapseq_exchange, str) else ""
            if g.strip():
                p2g[r.exchange] = g.strip()
        p = panel[(panel.source == "bacdive_metabolites") & (panel.Organism == org)]
        p = p[p[col].notna() & (p[col].astype(str).str.strip() != "")]
        mapped = sum(1 for e in p.Exchange.unique() if e in p2g)
        print("nomenclature map: %d of %d panel exchanges resolved to a gapseq compound"
              % (mapped, p.Exchange.nunique()))
        tp = fp = tn = fn = 0
        agree = n = 0
        for _, o in p.iterrows():
            e = ex.get(p2g.get(o.Exchange, ""))
            up = str(o.Direction).lower().startswith("up")
            call = (e["can_uptake"] if up else e["can_secrete"]) if e else False
            obs = int(o.observed) == 1
            mine = "TP" if (call and obs) else "FP" if (call and not obs) else \
                   "FN" if ((not call) and obs) else "TN"
            tp += mine == "TP"; fp += mine == "FP"; tn += mine == "TN"; fn += mine == "FN"
            agree += mine == str(o[col]).strip(); n += 1
        if n:
            print("ops=%d | ours TP %d FP %d TN %d FN %d acc %.4f | agreement with published %s: %d/%d"
                  % (n, tp, fp, tn, fn, (tp + tn) / n, col, agree, n))
    out_path = Path(tempfile.gettempdir()) / ("gapseq_scan_%s.tsv" % (org or "model"))
    pd.DataFrame([dict(exchange=k, **v) for k, v in ex.items()]).to_csv(
        out_path, sep="\t", index=False)
    print(f"wrote {out_path}")

if __name__ == "__main__":
    main()
