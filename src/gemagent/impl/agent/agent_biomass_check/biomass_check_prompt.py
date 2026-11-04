#!/usr/bin/env python3
"""biomass_check_prompt.py — 构建 biomass 逐前体审查 prompt（v4：全量前体、role 降为属性）。

设计要点（v4 相对 v3 的变化）：

  v3 的问题：先按 role 把前体分成 [AUTO_KEEP] / [AUDIT] 两段，AUTO_KEEP 只作
  上下文、明令不出决策；而 AUDIT 的四类 role 各自套一条决策规则，其中
  REDOX_CATALYTIC_POOL 那条（规则 3）是唯一带宽容 DROP 出口的，三个子条件里
  (a) 数值上不可判别、(b) 免证、(c) 可被幻觉满足。实测 bsub 的 Menaquinone-8
  正是从这条出口被误删，直接导致 F1 从 0.4348 掉到 0.3689。

  v4 的变化：
    1. role 降级为「属性」—— 仍显示，但不再决定谁被审、走哪条规则。因此关键词
       命中的先后顺序、以及「既是醌又是辅因子」这类本体嵌套，都不再影响判定。
    2. 全量前体都要求一个判定（KEEP | DROP | UNCERTAIN），不再有 AUDIT 分段。
    3. 规则收敛为一条统一的 DROP 判据：必须同时满足
         (i)   是周转性的催化/氧化还原/辅因子池；
         (ii)  D1 / D2 / D3 至少一条成立；
         (iii) reason 里具名给出该条要求的证据。
       D3（培养基可清除）由证据预计算并显式标 YES/NO，LLM 不得自行推理培养基。
    4. related_kos（本生物自身基因的 KOfam 注释）提升为「最强模型内证据」，
       非空时 D2 附加反驳义务。
    5. AUTO_KEEP 那三类 role 的 DROP 不在此处拦截，而由 apply_biomass_decisions.py
       的硬闸负责：记红旗、不落地。prompt 里如实告知，鼓励 LLM 给出真实判断
       作为诊断信号。

证据全部参考-free: role + 化学计量系数 + 跨模板保守性 + 培养基可得性 +
draft 自身基因的 KOfam 注释。

System prompt 复用 code/prompt/agent_prompts.py 的 BASE_BIOLOGY_CONTEXT
（其中 {mimedb_profile} 注入 MIMEdb 生理/栖息地/氧需求信息）。

输出: <out_dir>/prompt/<group>.txt + manifest.tsv
用法:
  python3 biomass_check_prompt.py \
      --evidence <out_dir>/biomass_evidence.json \
      --output_dir <pipeline_output_dir> \
      --out_dir <pipeline_output_dir>/biomass_check \
      --org_name "Escherichia coli"
"""
import argparse, csv, json, os, re, sys
from collections import OrderedDict

_HERE = os.path.dirname(os.path.abspath(__file__))
_CODE_DIR = os.path.normpath(os.path.join(_HERE, "..", ".."))
_PROMPT_DIR = os.path.normpath(os.path.join(_HERE, "..", "..", "prompt"))
for _p in (_PROMPT_DIR, _CODE_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import taxonomy  # noqa: E402  (code/taxonomy.py owns the naming tables)

BIOMASS_ROLE = """
[ROLE]
You are a Biomass Requirement Auditor. You audit the biomass reaction of a
draft genome-scale metabolic model (BEFORE gapfilling) and decide, PER
PRECURSOR, whether the model should REQUIRE that precursor for growth.

[THE REAL PROBLEM]
Gapfilling reads the biomass reaction and adds reactions to make every required
precursor producible. An over-requiring biomass reaction therefore has TWO
harmful effects:
  (1) gapfilling inserts spurious biosynthetic pathways, and those pathways can
      act as bypasses that mask true essential genes;
  (2) it manufactures false essential genes — the genes of such a pathway look
      essential only because the biomass demands their product.

[OBSERVABLES WE COMPUTE FOR YOU]
  PRODUCIBLE = can the draft model synthesize this precursor (FBA demand probe)?
  GAP        = NOT producible AND NOT supplied by the medium.
A GAP is exactly what gapfilling will try to fix.

[HOW TO READ THE EVIDENCE]
- role is INFORMATIONAL ONLY. It describes what kind of component this is. It
  does NOT select a rule: every precursor is judged by the SAME framework below.
  Two roles may both fit one metabolite (e.g. menaquinone is both a quinone pool
  and a cofactor); that overlap is chemistry, not a decision for you to resolve.
- cross-template conservation: which gapseq biomass templates contain this
  precursor. Presence in BOTH bacterial templates is the strongest conservation
  signal available here and supports KEEP. Absence from the Archaeal template
  says nothing about a bacterium.
- organism-own genes (related_kos): a WEAK reference signal. These are KOfam
  best hits whose annotation TEXT happens to contain a word from this precursor's
  name. They are keyword hits, NOT a curated biosynthetic pathway assignment: a
  hit can be an unrelated enzyme that merely shares a generic word (e.g.
  "phosphate"), and transport annotations are excluded on purpose because a
  transporter implies UPTAKE, which argues against self-synthesis. Treat a
  non-empty list as a hint worth checking, never as proof of a requirement.
- same-class analogues: other members of the same compound CLASS that already
  exist in the model (e.g. menaquinone-7 when the template asks for
  menaquinone-8). INFORMATIONAL ONLY: no substitution is possible in this
  run, so their presence simply shows the class is part of the model's
  chemistry — which argues for KEEP.

[DECISION FRAMEWORK]
Emit exactly ONE verdict per precursor: KEEP | DROP | UNCERTAIN.
Default is KEEP. Do not invent evidence: never cite an id or a fact that is not
in this prompt. UNCERTAIN means "I cannot establish this" and is treated as KEEP
downstream.

A DROP must satisfy ALL THREE of:
  (i)   the precursor is a catalytic / redox / cofactor pool that TURNS OVER,
        not accumulated cell mass; AND
  (ii)  criterion D1 or D3 below holds (D2 is CLOSED — see below); AND
  (iii) your "reason" names the evidence that criterion requires.

  D1 — LINEAGE MARKER
       The precursor is a diagnostic marker of a lineage this organism is not in.
       REQUIRED in "reason": the marker, AND the lineage it diagnoses.
         valid  : "teichoic acid is a Gram-positive wall marker; this organism
                   is Gram-negative."
         invalid: "this looks like a Gram-positive lipid."   (no marker named)

  D2 — IN-MODEL ALTERNATIVE   >> CLOSED IN THIS RUN
       This criterion would let you cite a functionally equivalent alternative
       already present in the model. We do NOT supply the model's reaction list
       in this prompt, so no citable in-model ID can be produced and any ID you
       wrote would be fabricated. D2 is therefore CLOSED: you may not DROP on D2
       grounds. If you believe a genuine in-model alternative exists, answer
       UNCERTAIN and name the suspected alternative in "reason" so a human can
       check it.

  D3 — SCAVENGABLE FROM THE MEDIUM
       The precursor is available in the predicted medium at non-zero flux.
       PRE-EVALUATED FOR YOU — see the field "D3 usable" on each precursor. If
       it says NO, D3 is closed and you may NOT invoke it. Do not reason about
       the medium yourself.

[NO SUBSTITUTION — there is no "swap in a different form" verdict]
A template may name a different FORM of a compound class than this organism makes
(for instance a different chain length of the same family). That is NOT a reason
to remove the requirement: keeping it only costs gapfilling a few extra reactions,
whereas dropping it deletes a real requirement and can erase genuinely essential
genes. So when you can see the class IS required, answer KEEP even if the specific
form looks like another organism's variant; when you cannot tell whether the class
is required, answer UNCERTAIN. Never invent a substitute compound.

[HARD GUARD — enforced downstream, do not rely on it]
Precursors whose role is MACROMOLECULE_MONOMER / INORGANIC /
ENERGY_REDOX_CURRENCY are NEVER removed, whatever verdict you emit:
such a verdict is logged as a red flag and ignored. Emit your honest verdict
anyway — we want the signal, not agreement.

[STRICT OUTPUT FORMAT]
Output ONLY valid JSON, no markdown, no extra text. Every string must be valid
JSON: never put an unescaped double quote or a raw newline inside "reason".
{"group": "<group name>", "decisions": [
  {"cpd": "<cpdXXXXX>", "name": "<name>",
   "verdict": "KEEP|DROP|UNCERTAIN",
   "reason": "<= 40 words naming the decisive evidence"}]}
Emit exactly one entry for EVERY precursor listed below. Do not skip any.
"""


def load_base_context():
    try:
        from agent_prompts import BASE_BIOLOGY_CONTEXT
        return BASE_BIOLOGY_CONTEXT
    except Exception as e:
        print(f"  [warn] 无法导入 BASE_BIOLOGY_CONTEXT ({e})")
        return ("You are an expert genome-scale metabolic model reconstruction "
                "specialist with deep biological and metabolic knowledge.\n\n"
                "[TARGET ORGANISM CONTEXT]\nThe organism annotated by GTDB-Tk is: "
                "{org_name}\n{mimedb_profile}\n"
                "--------------------------------------------------")


def load_taxonomy(out_root, override=""):
    """GTDB lineage for the prompt context.

    Explicit ``--taxonomy_tree`` wins; otherwise read the run's
    ``taxonomy/taxonomy.tsv``, which code/taxonomy.py is the sole writer of.
    """
    if override:
        return override
    return taxonomy.lineage_of(out_root)


def load_mimedb(out_root):
    p = os.path.join(out_root, "prompts", "pathway_verification_system.txt")
    if not os.path.exists(p):
        return ""
    text = open(p, encoding="utf-8").read()
    s = text.find("[TARGET ORGANISM CONTEXT]")
    if s < 0:
        return ""
    # 终止符必须是下一个 section 头，**不能**用 "\n---"：MIMEdb 块内部自带
    # "---" 分隔（Taxonomic Lineage / Aggregated Physiological Traits /
    # Representative Description），用 "---" 会在第一段标题后就把块截断，
    # 革兰氏染色、氧需求、栖息地等表型全部丢失（实测块长 226 -> 2575 字符）。
    e = -1
    for marker in ("\n[ROLE]", "\n[STRICT OUTPUT", "\n[INSTRUCTION",
                   "\n[OUTPUT", "\n[YOUR TASK]", "\n[TASK]"):
        k = text.find(marker, s + 30)
        if k > s and (e < 0 or k < e):
            e = k
    if e < 0:
        e = len(text)
    block = text[s:e].strip()
    # 上游 MIMEdb 的 strain 字段缺失时会被 format 成字面 "nan"（见
    # code/prompt/org_prompt_gen.py 第 42 行），这里做一次防御性清理。
    block = re.sub(r"using\s+nan\s+as the primary representative",
                   "consolidated across strains", block, flags=re.I)
    block = re.sub(r"^\[TARGET ORGANISM CONTEXT\]\s*", "", block)
    block = re.sub(r"\n?The organism annotated by GTDB-Tk is:.*?\n", "\n", block, count=1)
    return block.strip()


def organism_context(org_name, taxonomy, mimedb):
    parts = ["=" * 72, "ORGANISM CONTEXT", "=" * 72, f"Organism: {org_name}"]
    if taxonomy:
        parts.append(f"GTDB-Tk Taxonomy: {taxonomy}")
    if mimedb:
        parts += ["", mimedb]
    return "\n".join(parts)


def yn(v):
    return {True: "YES", False: "NO", None: "?"}.get(v, "?")


def normalise_catalog(cat):
    """归一化模板目录。catalog 来自 evidence 的 template_catalog；
    为兼容旧证据，也接受纯 id 列表（非 Archaea 视为 Bacteria）。"""
    out = []
    for t in cat or []:
        if isinstance(t, dict):
            out.append({"id": t.get("id"),
                        "name": t.get("name", ""),
                        "domain": t.get("domain") or "Unknown"})
        else:
            out.append({"id": t, "name": "",
                        "domain": "Archaea" if t == "Archaea" else "Bacteria"})
    return out


def conservation_line(catalog, present):
    """把跨模板分布渲染成保守性信号（而非排除信号）。

    域(Bacteria/Archaea)取自模板自身的 name 字段（见 normalise_catalog），
    因此新增/删除模板无需改这里。措辞刻意避免「conserved across ALL bacterial
    templates」这类读起来自相矛盾的表达（细菌模板只有 2 个时该说法尤其费解），
    改成「在 N 个细菌模板中的 M 个里出现」+ 明确列出缺失者。
    """
    catalog = normalise_catalog(catalog)
    present = list(present or [])
    all_ids = [t["id"] for t in catalog]
    bact_all = [t["id"] for t in catalog if t["domain"] == "Bacteria"]
    bact_here = [t for t in present if t in bact_all]
    nb = len(bact_all)
    if nb == 0:
        label = "no bacterial template in catalog"
    elif nb == 1:
        label = ("present in the only bacterial template (weak signal)"
                 if bact_here else "absent from the only bacterial template")
    elif len(bact_here) == nb:
        label = "broad conservation across bacterial templates (supports KEEP)"
    elif len(bact_here) >= 2:
        label = "shared by most bacterial templates"
    elif bact_here:
        label = "restricted: only one bacterial template contains it"
    else:
        label = "absent from every bacterial template"
    n_tot = len(all_ids)
    here = ", ".join(present) if present else "unmatched"
    out = [f"- cross-template conservation: present in {len(present)} of "
           f"{n_tot if n_tot else '?'} gapseq biomass templates ({here})"]
    if bact_all:
        out.append(f"    - among the {nb} BACTERIAL template(s) "
                   f"({', '.join(bact_all)}): {len(bact_here)} of {nb} contain it"
                   f"  -> {label}")
    miss = [t for t in all_ids if t not in present]
    out.append(f"    - absent from: {', '.join(miss) if miss else '(none)'}")
    return out


def precursor_block(r, catalog=None):
    L = []
    L.append(f"### {r['cpd']} — {r['met_name']}")
    L.append(f"- role (INFORMATIONAL only — does not select a rule): {r['role']}")
    rm = r.get("role_matches") or []
    if len(rm) > 1:
        L.append(f"    - role candidates (overlapping by chemistry; the first is the "
                 f"historical label, none of them selects a rule): {', '.join(rm)}")
    gm = r.get("template_group_mass")
    gmass = f"{gm}" if gm is not None else "n/a (not a template component)"
    L.append(f"- template group: {r.get('template_group') or '_not_in_template'}  "
             f"(group mass fraction {gmass})")
    coef = r.get("coef")
    if coef is None:
        L.append("- stoichiometric coefficient in biomass: n/a")
    else:
        L.append(f"- stoichiometric coefficient in biomass: {coef:.6f}")
    L += conservation_line(catalog or [], r.get("present_in_templates"))
    L.append(f"- producible by draft model (with predicted medium): {yn(r['producible'])}")
    gap = "YES  GAPFILLING WILL FILL THIS" if r["gap"] else "NO"
    L.append(f"- is this a GAP (not producible & not in medium): {gap}")
    L.append(f"- already supplied by the medium: {yn(r['from_medium'])}")
    if r.get("d3_usable"):
        L.append("- D3 usable: YES  (present in the predicted medium)")
    else:
        L.append("- D3 usable: NO  (not in the predicted medium — D3 is CLOSED for "
                 "this precursor; you may NOT invoke D3)")
    an = r.get("analogue_candidates") or []
    if an:
        L.append("- same-class analogues already in this model (INFORMATIONAL "
                 "ONLY — no substitution is possible in this run; their "
                 "presence argues the class is part of this chemistry):")
        for a in an:
            pr = {True: "producible", False: "not producible",
                  None: "producibility unknown"}.get(a.get("producible"), "?")
            L.append(f"    - {a['cpd']}  {a['name']}  ({pr})")
    rk = r.get("related_kos") or []
    if rk:
        L.append(f"- organism-own genes (related_kos, {len(rk)}) — WEAK reference "
                 f"signal only (keyword hits on KOfam text, not a curated pathway; "
                 f"transport annotations already excluded):")
        for g in rk[:8]:
            L.append(f"    - {g['gene']}  [{g['ko']}]")
    else:
        L.append("- organism-own genes (related_kos): none — no signal either way")
    if r.get("hard_guard"):
        L.append("- HARD GUARD: this precursor is NEVER removed "
                 "downstream; such a verdict is logged and ignored. Emit your "
                 "honest verdict.")
    L.append("")
    return "\n".join(L)

def build_user_prompt(group, items, ev):
    n_gap = sum(1 for r in items if r["gap"])
    n_hard = sum(1 for r in items if r.get("hard_guard"))
    L = [f"## Biomass template group: {group}",
         f"- template: {ev.get('template_id')} ({ev.get('template_name')})",
         f"- draft model wild-type growth rate (with predicted medium): "
         f"{ev.get('wt_growth'):.4f}",
         f"- predicted medium has {ev.get('n_medium_compounds')} compounds",
         f"- precursors in this group: {len(items)}  (of which GAP: {n_gap})",
         f"- {n_hard} precursor(s) in this group are under HARD GUARD "
         f"(never removed or replaced downstream; such a verdict is logged as a "
         f"red flag and ignored)",
         "",
         f"## Precursors — emit exactly {len(items)} decisions",
         ""]
    for r in items:
        L.append(precursor_block(
            r, ev.get("template_catalog") or ev.get("all_template_ids")))
    L += ["## Task",
          f"Decide, for each of the {len(items)} precursor(s) above, whether the "
          "biomass reaction should REQUIRE it for this organism.", "",
          "The GAP flag marks precursors gapfilling will fill. For a GAP precursor, "
          "the real question is: is this requirement real (the organism genuinely "
          "needs it, so gapfilling SHOULD fill it) or an artifact (over-requiring, so "
          "DROP it and prevent gapfilling from inserting a spurious pathway)?", "",
          "Apply the [DECISION FRAMEWORK] from your instructions. Default is KEEP. "
          "A DROP requires (i) a catalytic/redox/cofactor pool that turns over, "
          "(ii) D1 or D3 (D2 is CLOSED, D3 only where marked usable), and "
          "(iii) a reason naming the required evidence. There is no substitution "
          "verdict: if the class is required, keep it even when the template names "
          "the wrong form; if you cannot tell, answer UNCERTAIN.", "",
          "Output ONLY valid JSON:",
          '{"group": "%s", "decisions": [{"cpd": "...", "name": "...", '
          '"verdict": "KEEP|DROP|UNCERTAIN", '
          '"reason": "..."}]}' % group]
    return "\n".join(L)

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--org_name", default="")
    ap.add_argument("--taxonomy_tree", default="")
    ap.add_argument("--groups", default="")
    ap.add_argument("--audit_only", action="store_true",
                    help="只输出含 AUDIT 前体的组")
    args = ap.parse_args()

    ev = json.load(open(args.evidence, encoding="utf-8"))
    taxonomy = load_taxonomy(args.output_dir, args.taxonomy_tree)
    mimedb = load_mimedb(args.output_dir)
    print(f"taxonomy: {'OK' if taxonomy else 'EMPTY'}   mimedb: {'OK' if mimedb else 'EMPTY'}")

    system = load_base_context().format(org_name=args.org_name, mimedb_profile=mimedb) + BIOMASS_ROLE
    gc = organism_context(args.org_name, taxonomy, mimedb)

    groups = OrderedDict()
    for r in ev["precursors"]:
        groups.setdefault(r.get("template_group") or "_not_in_template", []).append(r)
    if args.groups:
        want = {g.strip() for g in args.groups.split(",")}
        groups = OrderedDict((k, v) for k, v in groups.items() if k in want)

    order = ["DNA", "RNA", "Protein", "Inorganics", "Cofactors", "Lipids", "Cellwall"]
    groups = OrderedDict(sorted(groups.items(),
                                key=lambda kv: (order.index(kv[0]) if kv[0] in order else 99, kv[0])))

    pdir = os.path.join(args.out_dir, "prompt")
    os.makedirs(pdir, exist_ok=True)
    manifest = []
    for g, items in groups.items():
        n_audit = sum(1 for r in items if r["audit_tier"].startswith("AUDIT"))
        if args.audit_only and n_audit == 0:
            continue
        user = build_user_prompt(g, items, ev)
        prompt = f"{system}\n\n---\n\n{gc}\n\n---\n\n{user}"
        fname = re.sub(r"[^\w.\-]", "_", g) + ".txt"
        open(os.path.join(pdir, fname), "w", encoding="utf-8").write(prompt)
        manifest.append((g, len(items), n_audit, fname))
        print(f"  {g:18s} 前体 {len(items):3d} (AUDIT {n_audit:3d}) -> {fname}")

    with open(os.path.join(pdir, "manifest.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["group", "n_precursors", "n_audit", "prompt_file"])
        w.writerows(manifest)

    tc = {}
    for r in ev["precursors"]:
        tc[r["audit_tier"]] = tc.get(r["audit_tier"], 0) + 1
    print(f"\n输出 {pdir} | prompt {len(manifest)} | 模板 {ev.get('template_id')} "
          f"| 前体 {ev.get('n_precursors')}")
    print("tier 分布: " + ", ".join(f"{k}={v}" for k, v in sorted(tc.items())))
    print(f"需 LLM 判定的前体 (AUDIT): {sum(v for k, v in tc.items() if k.startswith('AUDIT'))}")


if __name__ == "__main__":
    main()
