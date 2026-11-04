#!/usr/bin/env python3
"""
apply_biomass_decisions.py — 把 LLM 的 biomass 判定落到模型里。

（pipeline 步骤：biomass check 的第 3 步）

做的事:
  1. 读 biomass 模板 JSON + LLM 判定（biomass_decisions.json）
  2. 生成“裁剪版”模板 JSON：删掉判为 EXCLUDE 的 component（含 link 引用）
  3. 调 code/draft_construction/replace_biomass.R 把新 bio1 写回模型 ——
     该脚本内部经 parse_BMjson.R 按质量分数重新归一化化学计量数
     （Scoef := group_mass * coef / sum(mass*coef) * 1000），
     因此不需要手工配平。
  4. 输出 diff 报告 + 新模型 xml

用法:
  python3 apply_biomass_decisions.py \
      --draft_rds     <genome_id>-draft.RDS \
      --biomass_json  database/biomass/biomass_Gram_neg.json \
      --decisions     <out_dir>/biomass_decisions.json \
      --out_dir       <out_dir> \
      --out_rds       <genome_id>-draft-biomassadj.RDS \
      [--dry_run] [--keep_uncertain]
"""
import argparse, json, os, re, subprocess, sys

# <root>/code 提供共享的 paths（本文件在 <root>/code/agent/agent_biomass_check/）。
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import paths  # noqa: E402


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--draft_rds", required=True)
    ap.add_argument("--biomass_json", required=True)
    ap.add_argument("--decisions", required=True, help="biomass_decisions.json")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--out_rds", default="", help="默认 <draft> 同目录 -draft-biomassadj.RDS")
    ap.add_argument("--replace_biomass_r", default="", help="replace_biomass.R 路径（默认自动定位）")
    ap.add_argument("--exclude", default="", help="可选：直接给 cpd 黑名单（逗号分隔），覆盖 decisions")
    ap.add_argument("--evidence", default="",
                    help="biomass_evidence.json（硬闸 role 判定用）")
    ap.add_argument("--dry_run", action="store_true", help="只生成裁剪 JSON，不调 R")
    ap.add_argument("--keep_uncertain", action="store_true", default=True,
                    help="UNCERTAIN 视为保留（默认）")
    return ap.parse_args()


HARD_GUARD_ROLES = {"MACROMOLECULE_MONOMER", "INORGANIC", "ENERGY_REDOX_CURRENCY",
                    "NOT_A_PRECURSOR"}

# role 缺失时的退化代理：这几组按 classify() 的前置分支必然是硬闸类
HARD_GUARD_GROUPS = {"DNA", "RNA", "Protein", "Inorganics"}


def norm(s):
    return re.sub(r"[\s\-_]+", "", str(s or "").lower())


def lookup_met_names(seed_tsv):
    """{norm(cpd): name}，用于给 REPLACE 的替换物取名。"""
    out = {}
    if not seed_tsv or not os.path.exists(seed_tsv):
        return out
    with open(seed_tsv, encoding="utf-8", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            i_id, i_nm = hdr.index("id"), hdr.index("name")
        except ValueError:
            return out
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) > max(i_id, i_nm) and c[i_id]:
                out[norm(c[i_id])] = c[i_nm]
    return out


def load_roles(evidence_path, out_dir):
    """从 biomass_evidence.json 读 {norm(cpd): role} 和可用的类似物候选。"""
    cands = [evidence_path,
             os.path.join(out_dir, "biomass_check", "biomass_evidence.json"),
             os.path.join(out_dir, "biomass_evidence.json")]
    for c in cands:
        if not c or not os.path.exists(c):
            continue
        try:
            j = json.load(open(c, encoding="utf-8"))
        except Exception as e:
            print(f"  [warn] 读不了 {c}: {e}")
            continue
        roles, allowed = {}, {}
        for r in j.get("precursors", []):
            if not r.get("cpd"):
                continue
            k = norm(r["cpd"])
            roles[k] = r.get("role") or ""
            allowed[k] = {norm(a["cpd"]) for a in (r.get("analogue_candidates") or [])}
        print(f"硬闸证据: {os.path.basename(c)} ({len(roles)} 个 cpd 的 role)")
        return roles, allowed
    print("[warn] 未找到 biomass_evidence.json —— 硬闸退化为按 template_group 判断")
    return {}, {}


def collect_excludes(decisions_path, explicit):
    """返回 {norm(cpd): (cpd, name, group, reason)}"""
    ex = {}
    if explicit:
        for c in explicit.split(","):
            c = c.strip()
            if c:
                ex[norm(c)] = (c, "", "", "explicit --exclude")
        return ex
    d = json.load(open(decisions_path, encoding="utf-8"))
    for g, v in d.items():
        # 跳过 manifest 之类的非判定条目；_not_in_template 是真实组名，必须保留，
        # 否则该组的判定会被静默丢弃（它的 DROP 需要进红旗通道）。
        if not isinstance(v, dict):
            continue
        for it in v.get("decisions", []):
            if str(it.get("verdict", "")).upper() == "DROP":
                cpd = it.get("cpd", "")
                if cpd:
                    ex[norm(cpd)] = (cpd, it.get("name", ""), g, it.get("reason", ""))
    return ex


def strip_links(comp, ex):
    """把 comp['link'] 里指向被删 cpd 的项去掉；返回新 link 或 None。"""
    link = comp.get("link")
    if not link:
        return link
    keep = []
    for tok in str(link).split("|"):
        met = tok.split(":")[0].strip()
        if norm(met) in ex:
            continue
        keep.append(tok)
    return "|".join(keep) if keep else None


# 说明：本流水线**不做替换**（只有 KEEP / DROP / UNCERTAIN）。曾经实现过一条
# D4/REPLACE 通路（把模板里的错误形态换成同类正确形态），在 bsub 上实测有害：
# 把 Phosphatidylglycerol-dioctadecanoyl 换成 cpd23600 后 biomass 要求了一个
# gapfilling 无法制造的化合物，模型直接不可生长（gapfill "no feasible solution"）。
# 而且静态判据无法区分两者（DB 生产反应 5 vs 5、该生物候选集内 0 vs 0）。
# 因此移除替换：只保留 DROP，DROP 只会放松约束，不可能造成不可行。


def main():
    a = parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    tmpl = json.load(open(a.biomass_json, encoding="utf-8"))
    ex = collect_excludes(a.decisions, a.exclude)
    print(f"模板: {a.biomass_json} ({tmpl.get('id')})")
    print(f"待排除 cpd: {len(ex)}")
    for k, (cpd, nm, g, rs) in sorted(ex.items()):
        print(f"  - {cpd:10s} {nm[:38]:38s} [{g}] {rs[:60]}")
    replaced = []   # 替换通路已移除，恒为空；保留 `_audit.replaced` 以兼容下游

    # ── 硬闸：核心 role 的 DROP 不落地，只记红旗 ─────────────────────────
    # 理由：这几类是判错代价最大、合法该删概率最低的一组。LLM 的判定原样写进
    # biomass_red_flags.tsv 作为诊断信号。
    roles, _allowed = load_roles(a.evidence, a.out_dir)
    guarded = []
    for k in list(ex):
        cpd, nm, g, rs = ex[k]
        role = roles.get(k, "")
        if role in HARD_GUARD_ROLES or (not role and g in HARD_GUARD_GROUPS):
            guarded.append((cpd, nm, g, role or f"<{g}>", "DROP", rs))
            del ex[k]

    if guarded:
        print(f"\n[红旗] {len(guarded)} 条判定落在硬闸 role 上 -> 已忽略，不落地:")
        for cpd, nm, g, role, kind, rs in guarded:
            print(f"  ! {cpd:10s} {nm[:30]:30s} role={role:24s} {kind:16s} {rs[:44]}")
        rf = os.path.join(a.out_dir, "biomass_red_flags.tsv")
        with open(rf, "w", encoding="utf-8") as f:
            f.write("cpd\tname\tgroup\trole\tkind\treason\n")
            for cpd, nm, g, role, kind, rs in guarded:
                f.write(f"{cpd}\t{nm}\t{g}\t{role}\t{kind}\t{rs}\n")
        print(f"  已写 {rf}")

    print(f"落地排除 cpd: {len(ex)}")

    removed, kept_links = [], []
    dropped_groups = []
    new_groups = []
    for g in tmpl.get("met_groups", []):
        gname = g.get("group_name")
        new_comps = []
        group_removed = []
        for c in g.get("components", []):
            cpd = c.get("id", "")
            if norm(cpd) in ex:
                group_removed.append((gname, cpd, c.get("name")))
                continue
            c2 = dict(c)
            nl = strip_links(c2, ex)
            if nl != c2.get("link"):
                if nl is None:
                    c2.pop("link", None)
                    kept_links.append((gname, cpd))
                else:
                    c2["link"] = nl
            new_comps.append(c2)
        if not new_comps:
            # 整组的组分全部被判 DROP -> 接受该判定，整组移除。
            #
            # 旧行为是「保留原样以免破坏 biomass 定义」，这让 LLM 的判定静默失效：
            # mgen(M. genitalium，无细胞壁)被误分到 Gram_pos 模板后，LLM 正确地
            # 把 18 个 Cellwall 组分全判 DROP，但这里整组恢复，于是 gapfill 又补进
            # 21 个细胞壁反应，模型不生长、整条 GE 评测作废。
            #
            # 「避免空 biomass」由下面的核心组断言保证，而不是靠恢复被删组分。
            print(f"  [drop-group] 组 {gname} 的 {len(group_removed)} 个组分"
                  f"全部被判 DROP → 整组移除")
            dropped_groups.append(gname)
            removed.extend(group_removed)
            continue
        g2 = dict(g)
        g2["components"] = new_comps
        # MOLFRACTION 组若被删成员，提示 parse_BMjson 会自动 rescale
        new_groups.append(g2)
        removed.extend(group_removed)

    # 安全阀：整组移除不得把 biomass 掏空到没有核心大分子组。
    CORE_GROUPS = {"DNA", "RNA", "Protein"}
    kept_names = {g.get("group_name") for g in new_groups}
    missing = CORE_GROUPS - kept_names
    if missing:
        sys.exit(f"错误: 整组移除后 biomass 缺少核心组 {sorted(missing)}。"
                 f"被整组移除的是 {dropped_groups}，共删除 {len(removed)} 个组分。"
                 f"请检查 {a.decisions} 是否把核心组误判为 DROP。")


    out = dict(tmpl)
    out["met_groups"] = new_groups
    out["_audit"] = {
        "source_template": os.path.abspath(a.biomass_json),
        "removed": [{"group": g, "cpd": c, "name": n} for g, c, n in removed],
        "n_removed": len(removed),
        "dropped_groups": dropped_groups,
        "n_dropped_groups": len(dropped_groups),
        "replaced": [{"group": g, "from": a_, "from_name": b_,
                      "to": c_, "to_name": d_} for g, a_, b_, c_, d_ in replaced],
        "n_replaced": len(replaced),
        "note": "stoichiometry is renormalised by parse_BMjson.R on render",
    }
    adj = os.path.join(a.out_dir, f"biomass_{tmpl.get('id', 'template')}_adjusted.json")
    json.dump(out, open(adj, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print(f"\n裁剪模板: {adj}  (删除 {len(removed)} 个 component)")
    if dropped_groups:
        print(f"  整组移除: {dropped_groups}")
    for g, c, n in removed:
        print(f"   - [{g}] {c} {n}")
    if kept_links:
        print(f"   另有 {len(kept_links)} 个 link 引用被清理: {kept_links}")

    # 变化摘要（component 级）—— removed 现在是「真正写进裁剪模板」的集合，
    # 不再包含被整组恢复的组分，因此报告与模板始终一致。
    with open(os.path.join(a.out_dir, "biomass_apply_report.tsv"), "w",
              encoding="utf-8") as f:
        f.write("group\tcpd\tname\treason\n")
        for g, c, n in removed:
            f.write(f"{g}\t{c}\t{n}\t{ex.get(norm(c), ('', '', '', ''))[3]}\n")

    # 变化摘要（group 级）
    with open(os.path.join(a.out_dir, "biomass_apply_groups.tsv"), "w",
              encoding="utf-8") as f:
        f.write("group\tn_template\tn_kept\tstatus\n")
        for g in tmpl.get("met_groups", []):
            gname = g.get("group_name")
            n0 = len(g.get("components", []))
            if gname in dropped_groups:
                f.write(f"{gname}\t{n0}\t0\tDROPPED\n")
            else:
                kept = next((len(x.get("components", []))
                             for x in new_groups if x.get("group_name") == gname), 0)
                f.write(f"{gname}\t{n0}\t{kept}\t"
                        f"{'KEPT' if kept == n0 else 'TRIMMED'}\n")

    if a.dry_run:
        print("\n--dry_run: 不调用 R，仅生成裁剪模板。")
        return

    # 定位 replace_biomass.R
    rscript_bin = paths.RSCRIPT
    rep = a.replace_biomass_r
    if not rep:
        cand = os.path.join(paths.CODE_DIR, "draft_construction",
                            "replace_biomass.R")
        rep = cand if os.path.exists(cand) else ""
    if not rep or not os.path.exists(rep):
        sys.exit(f"错误: 找不到 replace_biomass.R (--replace_biomass_r): {rep!r}")

    out_rds = a.out_rds or re.sub(r"\.RDS$", "", a.draft_rds) + "-biomassadj.RDS"
    cmd = [rscript_bin, rep, a.draft_rds, adj, out_rds]
    print(f"\n调用: {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout[-4000:] if r.stdout else "")
    if r.returncode != 0:
        print(r.stderr[-4000:], file=sys.stderr)
        sys.exit(f"replace_biomass.R 失败 (exit {r.returncode})")
    print(f"\n完成。新模型: {out_rds}")
    out_stem = re.sub(r"\.RDS$", "", out_rds)
    print(f"  SBML: {out_stem}.xml")


if __name__ == "__main__":
    main()
