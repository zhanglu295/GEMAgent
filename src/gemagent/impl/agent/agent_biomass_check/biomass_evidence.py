#!/usr/bin/env python3
"""
biomass_evidence.py — 为 biomass 逐前体审查收集证据（v3：gap 前置、无参考依赖）。

设计要点(回答"真实场景没有 ground truth / 参考模型怎么办")：

  1. 必须跑在 **draft(未 gapfill)** 模型上，因为 gapfill 会读取 biomass reaction
     来决定补哪些反应。若 biomass 过度要求某个前体，gapfill 就会去补一条本来
     不该存在的通路(错误 bypass 掩盖必需基因 = FN)，并人为制造该前体合成基因
     的"必需性"(= FP)。所以审计必须发生在 gapfill 之前。

  2. 施加 **预测培养基**({genome_id}-medium.csv，gapfill 正是用它)：把对应
     EX_cpdXXXX_e0 的 uptake 下界设为 -maxFlux。

  3. 主判据 = **可生产性(gap 检测)**，纯 FBA，无需 ground truth：
       producible = 给该前体加 demand 探针，看 draft+medium 能否合成。
       若不可生产且培养基未提供 → 这是一个 gap，gapfill 会去补。
       交给 LLM 判断：这个 gap 是「真缺失」(保留，让 gapfill 补) 还是
       「过度要求」(删掉，别让 gapfill 乱补)。

  4. 基因注释只用 **KOfam**(模型自身基因 ID → KO definition)，不用任何
     参考模型的 gene symbol。

输出: <out_dir>/biomass_evidence.json
用法:
  python3 biomass_evidence.py \
      --model_xml   <genome_id>-draft.xml \
      --biomass_json <database>/biomass/biomass_Gram_neg.json \
      --medium_csv  <genome_id>-medium.csv \
      --out_dir     <output_dir>/biomass_check \
      [--kofam <genome_id>_kofam.tsv] [--template_dir <biomass dir>]
"""
import argparse, csv, json, os, re, sys
from collections import defaultdict

# ── 角色分类 ─────────────────────────────────────────────────────────────
MACROMOLECULE_GROUPS = {"DNA", "RNA", "Protein"}

NOT_PRECURSOR_NAMES = {
    "h2o", "water", "h+", "proton", "hydron", "ppi", "diphosphate",
    "co2", "carbon dioxide", "o2", "oxygen", "nh4", "ammonium", "nh3",
    "ammonia", "biomass", "acetate", "formate",
}
ENERGY_CURRENCY_KEYWORDS = (
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "ctp", "dttp", "dgtp",
    "nadh", "nadph", "nad", "nadp", "fadh2", "fad", "fmn", "coenzyme a",
    "coa", "acetyl-coa", "tetrahydrofolate", "thf", "s-adenosyl",
)
REDOX_POOL_KEYWORDS = (
    "menaquinon", "ubiquinon", "demethylmenaquinon", "plastoquinon",
    "glutathion", "spermidin", "putrescin", "spermin", "cadaverin",
    "thioredoxin", "glutaredoxin", "mycothiol", "ergothionein",
)
REDOX_POOL_EXACT = {"gsh", "gssg", "glutathione", "glutathione disulfide"}
STRUCTURAL_KEYWORDS = (
    "peptidoglycan", "murein", "lipid a", "lipid-a", "kdo", "lipo",
    "teichoic", "wall", "phosphatidyl", "cardiolipin", "lysylphosphatidyl",
    "glycan", "colanic", "capsul", "mannan", "glucan", "chitin", "s-layer",
    "polysaccharide", "oligosaccharide",
)
COFACTOR_KEYWORDS = (
    "heme", "siroheme", "protoheme", "pyridoxal", "thiamin", "riboflavin",
    "cobalamin", "cobinamid", "cobamid", "biotin", "pantothen", "folate",
    "dihydrofolate", "tetrahydrobiopterin", "molybdopterin", "lipoate",
    "lipoic", "ascorbate", "quinone", "bactoprenyl", "undecaprenyl",
    "dimethylbenzimidazole", "acyl-carrier",
)
COFACTOR_EXACT = {"tpp", "thiamine pyrophosphate", "acp", "apo-acp", "calomide"}


# 下游 apply 阶段的硬闸：这几类 role 的 DROP 永不落地（只记红旗）。
# 依据 prompt 里的 [CRITICAL ASYMMETRY]：判错 DROP 只会变成假阴性基因，而这几类
# 是「判错代价最大、合法该删概率最低」的一组（氨基酸/核苷酸/离子/ATP-NAD-FAD-CoA
# 类货币）。注意 role 不再参与裁决，硬闸是独立于 prompt 的第二道防线。
HARD_GUARD_ROLES = {"MACROMOLECULE_MONOMER", "INORGANIC", "ENERGY_REDOX_CURRENCY"}


def classify_matches(group_name, met_name, cpd_id=""):
    """返回该前体命中的**全部** role 候选（按历史优先级排序）。

    保留这个顺序只为让 classify() 的返回值与 v3 保持一致（向后兼容、便于对比）；
    v4 的 prompt 不再用 role 选规则，因此这个顺序不再影响任何判定结果。
    多个命中通常来自本体嵌套（例：Menaquinone-8 同时命中 menaquinon 和 quinone）。
    """
    g = (group_name or "").strip()
    n = (met_name or "").lower().strip()
    n = re.sub(r"-c\d+$", "", n)
    hits = []

    def add(r):
        if r not in hits:
            hits.append(r)

    if cpd_id == "cpd11416" or n in ("biomass",):
        add("NOT_A_PRECURSOR")
    if g in MACROMOLECULE_GROUPS:
        add("MACROMOLECULE_MONOMER")
    if g == "Inorganics":
        add("INORGANIC")
    if n in NOT_PRECURSOR_NAMES:
        add("NOT_A_PRECURSOR")
    if any(k in n for k in ENERGY_CURRENCY_KEYWORDS):
        add("ENERGY_REDOX_CURRENCY")
    if n in REDOX_POOL_EXACT or any(k in n for k in REDOX_POOL_KEYWORDS):
        add("REDOX_CATALYTIC_POOL")
    if any(k in n for k in STRUCTURAL_KEYWORDS):
        add("STRUCTURAL_ENVELOPE")
    if n in COFACTOR_EXACT or any(k in n for k in COFACTOR_KEYWORDS):
        add("COFACTOR_SMALL")
    if not hits:
        add("OTHER")
    return hits


def classify(group_name, met_name, cpd_id=""):
    """单个 role（= classify_matches 的第一个）。保留供 v3 路径与对比使用。"""
    return classify_matches(group_name, met_name, cpd_id)[0]


def audit_tier(role):
    if role in ("MACROMOLECULE_MONOMER", "INORGANIC", "ENERGY_REDOX_CURRENCY",
                "NOT_A_PRECURSOR"):
        return "AUTO_KEEP"
    if role == "REDOX_CATALYTIC_POOL":
        return "AUDIT_POOL"
    if role == "COFACTOR_SMALL":
        return "AUDIT_COFACTOR"
    if role == "STRUCTURAL_ENVELOPE":
        return "AUDIT_STRUCTURAL"
    return "AUDIT_OTHER"


def norm(s):
    return re.sub(r"[\s\-_]+", "", str(s or "").lower())


# ── 模板 / 培养基 / 其它模板 / KOfam ──────────────────────────────────────
def load_template(json_path):
    d = json.load(open(json_path, encoding="utf-8"))
    by_cpd, by_name = {}, {}
    for g in d.get("met_groups", []):
        for c in g.get("components", []):
            rec = {"group": g.get("group_name"), "group_mass": g.get("mass"),
                   "unit": g.get("unit_components"), "cpd": c.get("id"),
                   "name": c.get("name"), "coef": c.get("coef"),
                   "comp": c.get("comp"), "link": c.get("link")}
            if c.get("id"):
                by_cpd[norm(c["id"])] = rec
            if c.get("name"):
                by_name[norm(c["name"])] = rec
    return d, by_cpd, by_name


def scan_other_templates(biomass_dir, exclude_json):
    out = defaultdict(list)
    if not biomass_dir or not os.path.isdir(biomass_dir):
        return out
    for fn in sorted(os.listdir(biomass_dir)):
        if not fn.endswith(".json") or fn == os.path.basename(exclude_json):
            continue
        try:
            d = json.load(open(os.path.join(biomass_dir, fn), encoding="utf-8"))
        except Exception:
            continue
        tid = d.get("id") or fn
        for g in d.get("met_groups", []):
            for c in g.get("components", []):
                if c.get("id"):
                    out[norm(c["id"])].append(tid)
    return out


def load_medium(medium_csv):
    """返回 {norm(cpd): maxFlux}。"""
    provided = {}
    if not medium_csv or not os.path.exists(medium_csv):
        return provided
    with open(medium_csv, encoding="utf-8") as fh:
        rows = list(csv.reader((l for l in fh if not l.startswith("#"))))
    if not rows:
        return provided
    hdr = [h.strip().lower() for h in rows[0]]

    def col(*names):
        for nm in names:
            if nm in hdr:
                return hdr.index(nm)
        return None
    c_cpd, c_name = col("compounds", "compound", "cpd", "id", "modelseed"), col("name")
    c_flux = col("maxflux", "maxFlux", "flux", "uptake")
    for r in rows[1:]:
        if len(r) <= 1:
            continue
        flux = 0.0
        if c_flux is not None and c_flux < len(r):
            try:
                flux = float(r[c_flux] or 0)
            except ValueError:
                flux = 0.0
        for idx in (c_cpd, c_name):
            if idx is not None and idx < len(r) and r[idx].strip():
                provided[norm(r[idx])] = flux
    return provided


def load_kofam(kofam):
    """返回 {model_gene_id: [KO def, ...]}。"""
    info = defaultdict(list)
    if not kofam or not os.path.exists(kofam):
        return info
    for line in open(kofam, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if p and p[0] == "*" and len(p) >= 7:
            gid = p[1].replace(".", "_")
            # KO 定义字段本身含空格并被引号包裹，必须取 p[6:] 后再剥引号。
            # 只取 p[6] 会把定义截成一个词（例：K02552 只剩 '"menaquinone-specific'），
            # 而 v4 把 related_kos 当作承重证据，截断会直接削弱判定依据。
            definition = " ".join(p[6:]).strip().strip('"').strip()
            info[gid].append(f"{p[2]} {definition}")
    return info


# ── 培养基施加 + 可生产性 ────────────────────────────────────────────────
def apply_medium(model, medium):
    """把预测培养基的 maxFlux 设为对应 EX_cpdXXXX_e0 的 uptake 下界。"""
    n = 0
    for cpd_norm, flux in medium.items():
        if flux <= 0:
            continue
        for suffix in ("_e0", "_c0", "_p0"):
            exid = f"EX_{cpd_norm}{suffix}"
            if exid in model.reactions:
                model.reactions.get_by_id(exid).lower_bound = -flux
                n += 1
                break
    return n


# ── 同类类似物（链长 / 酰基链变体）───────────────────────────────────────
# 用于回答「模板要的是这一类的错误形态」：例如模板要求 menaquinone-8，而该生物
# 实际合成 menaquinone-7。核心词归一后**完全相等**才算同类，故意保守：
# "Menaquinone 8" -> "menaquinone"、"Menaquinone 7" -> "menaquinone"（同类）；
# "2-Demethylmenaquinone 8" -> "demethylmenaquinone"（不同类，且它是 MK 的合成
# 前体，绝不能被当成替代品）。
_CHAIN_STOP = {"n", "c", "b", "p", "link", "linked", "type", "form",
               "subtilis", "coli", "aureus", "aerobic", "anaerobic"}


def class_core(name):
    s = str(name or "").lower()
    s = re.sub(r"-[cpe]\d+$", "", s)                    # 区室后缀
    s = re.sub(r"\(.*?\)", " ", s)                      # 括号限定
    s = re.sub(r"\bn\s*=\s*\d+\b", " ", s)            # n=24
    s = re.sub(r"^\s*\d+(?:\s*,\s*\d+)*\s*-", " ", s)  # 前置定位数字 2-
    s = re.sub(r"[\d\s,]+", " ", s)                     # 全部数字
    toks = [t for t in re.split(r"[^a-z]+", s) if t and t not in _CHAIN_STOP]
    return "".join(toks)


def build_class_index(model):
    """{核心词: [(cpd, name)]}，取自模型全部代谢物。"""
    idx = defaultdict(list)
    for met in model.metabolites:
        core = class_core(met.name or met.id)
        if not core or len(core) < 6:
            continue
        idx[core].append((re.sub(r"_[a-z]\d+$", "", met.id), met.name or met.id))
    return idx


def class_analogues(model, class_idx, met):
    """同类、非自身、模型里确实存在的类似物，附可生产性（demand 探针）。"""
    core = class_core(met.name or met.id)
    base = re.sub(r"_[a-z]\d+$", "", met.id)
    out, seen = [], {norm(base)}
    for cpd, name in class_idx.get(core, []):
        if norm(cpd) in seen:
            continue
        seen.add(norm(cpd))
        pm = None
        for suffix in ("_c0", "_e0", "_p0"):
            try:
                pm = model.metabolites.get_by_id(cpd + suffix)
                break
            except Exception:
                continue
        prod = None
        if pm is not None:
            try:
                with model:
                    dm = model.add_boundary(pm, type="demand", reaction_id="DM_analog")
                    model.objective = dm
                    v = model.slim_optimize()
                prod = bool(v is not None and v > 1e-6)
            except Exception:
                prod = None
        out.append({"cpd": cpd, "name": name, "producible": prod})
    return out[:6]


def precursor_gap(model, met):
    """给前体加 demand 探针，看 draft+medium 能否合成。返回 (producible, val)。"""
    try:
        with model:
            dm = model.add_boundary(met, type="demand", reaction_id="DM_probe")
            model.objective = dm
            val = model.slim_optimize()
        return bool(val is not None and val > 1e-6), val
    except Exception:
        return None, None


# ── 相关 KO（弱参考信号，非通路归属）───────────────────────────────────────
_STOP = {"the", "with", "from", "into", "acid", "group", "dna", "rna",
         "protein", "polymer", "subunit", "core", "lipid", "phosphatidylethanolamine",
         "phosphatidylglycerol", "oligosaccharide", "dioctadecanoyl",
         "diisoheptadecanoyl", "dianteisoheptadecanoyl", "isoheptadecanoyl",
         "anteisoheptadecanoyl", "stearoyl", "bactoprenyl"}

# 通用化学/酶类词：会命中大量与该化合物无关的酶注释。实测 Pyridoxal-phosphate
# 用 "phosphate" 检索，8 条命中里 0 条是 PLP 合成酶（全是 G6PDH、磷酸转运蛋白、
# GlcNAc-6-P 脱乙酰酶之类）。这类词一律不作为检索关键词。
_GENERIC_CHEM = {
    "phosphate", "phosphates", "pyrophosphate", "diphosphate", "triphosphate",
    "monophosphate", "methyl", "methylene", "hydroxy", "hydroxyl", "amino",
    "oxo", "carbon", "oxide", "oxygen", "hydrogen", "sulfur", "sulphur",
    "sulfate", "sulphate", "nitrate", "nitrite", "ammonia", "iron", "copper",
    "zinc", "calcium", "magnesium", "potassium", "sodium", "chloride",
    "complex", "dependent", "containing", "putative", "family", "chain",
    "domain", "type", "related", "specific", "like", "unknown", "probable",
    "synthesis", "synthase", "reductase", "oxidase", "transferase",
    "dehydrogenase", "hydrolase", "isomerase", "kinase", "mutase", "lyase",
    "ligase", "carrier", "binding", "subfamily", "enzyme", "terminal",
}

# 转运类注释要剔除：转运蛋白的存在说明该化合物**可被摄取**，而不是必须自合成，
# 把它算作「该生物能造它」的证据方向恰好相反。
_TRANSPORT_HINTS = (
    "transport", "permease", "symport", "antiport", "abc ", "abc-",
    "channel", "translocat", "efflux", "uptake", "porin", "pump",
    "atp-binding cassette",
)


def is_transport_ko(text):
    low = text.lower()
    return any(h in low for h in _TRANSPORT_HINTS)


def rel_keywords(name, role):
    """只取前体自身名字里的**有鉴别力**的词，不做 role 级扩展。"""
    kw = set()
    n = re.sub(r"-[cpe]\d+$", "", (name or "").lower())
    n = re.sub(r"\b\d+\b", " ", n)           # 去掉 "8"、"n" 等
    for w in re.split(r"[^a-z]+", n):
        if len(w) >= 4 and w not in _STOP and w not in _GENERIC_CHEM:
            kw.add(w)
    return kw


def related_kos(model_genes, kofam, keywords, cap=8):
    """弱参考信号：模型基因的 KOfam best-hit 注释文本里出现这些关键词的条目。

    这是**文本关键词命中**，不是通路归属；转运类注释已剔除。
    """
    hits = []
    for gid in model_genes:
        for ko in kofam.get(gid, []):
            low = ko.lower()
            if is_transport_ko(low):
                continue
            if any(k in low for k in keywords):
                hits.append((gid, ko))
    # 去重(同 gene 多条 KO 取第一条)
    seen, uniq = set(), []
    for gid, ko in hits:
        if gid in seen:
            continue
        seen.add(gid)
        uniq.append({"gene": gid, "ko": ko})
    return uniq[:cap]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model_xml", required=True, help="draft SBML(未 gapfill)")
    ap.add_argument("--biomass_json", required=True)
    ap.add_argument("--biomass_rxn", default="bio1")
    ap.add_argument("--medium_csv", default="", help="预测培养基({genome_id}-medium.csv)")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--kofam", default="", help="KOFamScan tsv，注释相关基因")
    ap.add_argument("--template_dir", default="", help="biomass 模板目录(跨模板对比)")
    ap.add_argument("--groups", default="", help="只保留这些组(调试)")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    tmpl, by_cpd, by_name = load_template(args.biomass_json)
    medium = load_medium(args.medium_csv)
    other_tmpl = scan_other_templates(args.template_dir or os.path.dirname(args.biomass_json),
                                      args.biomass_json)
    kofam = load_kofam(args.kofam)
    print(f"模板 {tmpl.get('id')} | 培养基 {len(medium)} 化合物 | 其它模板 {len(other_tmpl)} cpd "
          f"| kofam {len(kofam)} 基因")

    import warnings
    warnings.filterwarnings("ignore")
    import cobra
    from cobra import Configuration
    Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(args.model_xml)
    if args.biomass_rxn not in model.reactions:
        sys.exit(f"模型找不到 {args.biomass_rxn}")
    n_applied = apply_medium(model, medium)
    bio = model.reactions.get_by_id(args.biomass_rxn)
    wt = model.slim_optimize()
    print(f"模型 {len(model.reactions)} rxn | {args.biomass_rxn} 前体 {len(bio.metabolites)} "
          f"| 施加 medium exchange {n_applied} | WT={wt:.4f}")

    # 全部 gapseq biomass 模板目录（供 prompt 计算跨模板保守性）。
    # 域(Bacteria/Archaea)由模板自身的 name 字段推导 —— 模板名不在此处硬编码，
    # 新增模板时自动纳入。
    _bdir = args.template_dir or os.path.dirname(args.biomass_json)
    tmpl_ids = []
    template_catalog = []
    if os.path.isdir(_bdir):
        for _fn in sorted(os.listdir(_bdir)):
            if not _fn.endswith(".json"):
                continue
            try:
                _d = json.load(open(os.path.join(_bdir, _fn), encoding="utf-8"))
            except Exception:
                continue
            _tid = _d.get("id") or _fn
            if _tid in tmpl_ids:
                continue
            _nm = str(_d.get("name") or "")
            _low = _nm.lower()
            if "archae" in _low:
                _dom = "Archaea"
            elif "bacteri" in _low:
                _dom = "Bacteria"
            else:
                _dom = "Unknown"
            tmpl_ids.append(_tid)
            template_catalog.append({"id": _tid, "name": _nm, "domain": _dom})

    model_gene_ids = {g.id for g in model.genes}
    class_idx = build_class_index(model)
    recs = []
    for met, coef in sorted(bio.metabolites.items(), key=lambda x: x[0].id):
        base = re.sub(r"_[a-z]\d+$", "", met.id)
        t = (by_cpd.get(norm(base)) or by_cpd.get(norm(met.id))
             or by_name.get(norm(met.name))
             or by_name.get(norm(re.sub(r"-c\d+$", "", str(met.name)))))
        role_matches = classify_matches((t or {}).get("group"), met.name or met.id,
                                        (t or {}).get("cpd") or base)
        role = role_matches[0]
        hard_guard = role in HARD_GUARD_ROLES
        producible, pval = precursor_gap(model, met)
        from_medium = (norm(base) in medium) or (norm(met.name) in medium) \
            or (t is not None and (norm(t["cpd"]) in medium or norm(t["name"]) in medium))
        gap = (producible is False) and (not from_medium)
        kws = rel_keywords(met.name, role)
        relkos = related_kos(model_gene_ids, kofam, kws) if kws else []
        others = other_tmpl.get(norm(base)) or other_tmpl.get(norm(met.id)) or []
        analogues = class_analogues(model, class_idx, met)
        recs.append({
            "met_id": met.id, "met_name": met.name or met.id,
            "cpd": (t or {}).get("cpd") or base, "coef": coef,
            "template_group": (t or {}).get("group"),
            "template_group_mass": (t or {}).get("group_mass"),
            "in_template": t is not None,
            "role": role, "role_matches": role_matches,
            "audit_tier": audit_tier(role), "hard_guard": hard_guard,
            "producible": producible,
            "producible_flux": pval,
            "from_medium": from_medium,
            "d3_usable": bool(from_medium),
            "gap": gap,
            "in_other_templates": others,
            "analogue_candidates": analogues,
            "d2_usable": False,   # 框架决定：prompt 不给反应清单，D2 关闭
            "present_in_templates": sorted(
                set(others + ([tmpl.get("id")] if tmpl.get("id") else []))),
            "related_kos": relkos,
        })

    if args.groups:
        want = {g.strip() for g in args.groups.split(",")}
        recs = [r for r in recs if r["template_group"] in want]

    out = {"biomass_rxn": args.biomass_rxn, "model_xml": os.path.abspath(args.model_xml),
           "template_json": os.path.abspath(args.biomass_json),
           "template_id": tmpl.get("id"), "template_name": tmpl.get("name"),
           "wt_growth": wt, "medium_csv": os.path.abspath(args.medium_csv) if args.medium_csv else "",
           "n_medium_compounds": len(medium),
           "n_precursors": len(recs),
           "hard_guard_roles": sorted(HARD_GUARD_ROLES),
           "all_template_ids": tmpl_ids, "n_templates_total": len(tmpl_ids),
           "template_catalog": template_catalog,
           "precursors": recs}
    p = os.path.join(args.out_dir, "biomass_evidence.json")
    json.dump(out, open(p, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"\n输出 {p}")

    tc = defaultdict(int)
    for r in recs:
        tc[r["audit_tier"]] += 1
    print("\naudit tier 分布:")
    for k in sorted(tc):
        print(f"  {k:18s} {tc[k]:3d}")

    n_gap = sum(1 for r in recs if r["gap"])
    n_gap_audit = sum(1 for r in recs if r["gap"] and r["audit_tier"].startswith("AUDIT"))
    print(f"\n总前体 {len(recs)} | 不可生产(gap) {n_gap} | 其中 AUDIT 类 gap {n_gap_audit}")
    print("\ngap 前体清单(按 tier 分组):")
    for r in sorted(recs, key=lambda x: x["audit_tier"]):
        if r["gap"]:
            print(f"  [{r['audit_tier']:16s}] {r['cpd']:12s} {str(r['met_name'])[:38]:38s} "
                  f"coef={r['coef']:9.5f}")


if __name__ == "__main__":
    main()
