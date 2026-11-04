#!/usr/bin/env python3
"""
Translate a foreign GEM's opaque gene IDs into gapseq-style coordinate names.

（cross-tool 接入的第 0 步）

背景
----
GEMAgent 的三套 review 与官方 gapseqEval 都以**基因组坐标**为基因的身份：
gapseq 把坐标直接写进基因名（`<contig>_<start>:<end>`），官方
`predGeneEssentiality.GS` 再从基因名里解析回来。

而 ModelSEED / CarveMe 用的是不透明的 locus tag（`6666666.524380.peg.77`、
`NP_391348_1`）。若原样喂进管线：
  - GPR review 的 KOFam 坐标桥拿不到坐标（KO 覆盖率 0 → review 空转）；
  - 官方 GE 评测器无法把 locus 映回模型基因。

本脚本做**一次性的、无损的坐标重标注**：用每个工具自带的基因位置表
（ModelSEED `*.gtbl` / CarveMe `*_genes_pos.tsv`）把 locus tag 换成
`<contig>_<start>:<end>`。这正是官方 `predGeneEssentiality.MS/.CM`
内部所用的同一份坐标映射，所以基因身份没有任何丢失——重标注后
`predGeneEssentiality.GS` 得到的结果与 `.MS`/`.CM` 逐 locus 一致。

坐标朝向
--------
两个位置表都**保留链方向**：负链基因 start > end（例如
`BG10944 3790436 3789189`），loci 表同样如此。官方评测的
`calc_seq_overlap` 用 `start:end` 生成序列，若把负链统一成 min:max，
交集会变空 → essential 全判否。因此这里**原样保留** (start,end) 朝向。

用法
----
  relabel_model.py --model_xml <in.xml> --kind modelseed|carveme \
      --pos_table <table> --out_xml <out.xml> \
      [--out_genepos <gene,contig,start,end tsv>] [--out_map <old,new tsv>]
"""
import argparse
import csv
import json
import os
import re
import sys

import cobra
import warnings

warnings.filterwarnings("ignore")


def _as_num(x):
    try:
        return int(float(str(x).strip()))
    except (TypeError, ValueError):
        return None


def load_positions_modelseed(gtbl_path):
    """ModelSEED RAST .gtbl -> {gene: (contig, start, end)}。

    `location` 每项是 [contig, sstart, sdir, slength]；结束坐标按
    `prepare_MS_gene_table.R` 的口径：'+' 用 sstart+slength，
    '-' 用 sstart-slengen（保持负链 start>end 的降序朝向）。

    一个 feature 可能有多个 location（断裂基因），取第一个，与官方一致。
    """
    with open(gtbl_path, encoding="utf-8", errors="replace") as fh:
        data = json.load(fh)
    out = {}
    for feat in data.get("features", []):
        gid = str(feat.get("id", "")).replace("fig|", "")
        loc = feat.get("location") or []
        if not gid or not loc:
            continue
        first = loc[0]
        if len(first) < 4:
            continue
        contig, sstart, sdir, slength = first[0], first[1], first[2], first[3]
        s = _as_num(sstart)
        ln = _as_num(slength)
        if s is None or ln is None:
            continue
        e = s + ln if str(sdir).strip() == "+" else s - ln
        out[gid] = (str(contig).strip(), s, e)
    return out


def load_positions_carveme(pos_path):
    """CarveMe `*_genes_pos.tsv` -> {gene: (contig, start, end)}。

    文件无表头，列序 gene,contig,start,end；朝向原样保留（负链 start>end）。
    容忍带表头/空白/逗号分隔。
    """
    out = {}
    with open(pos_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = re.split(r"[\t,]", line)
            if len(cols) < 4:
                cols = line.split()
            if len(cols) < 4:
                continue
            g, contig, s, e = cols[0], cols[1], cols[2], cols[3]
            si, ei = _as_num(s), _as_num(e)
            if si is None or ei is None:
                continue  # skips a header row too
            out[g] = (str(contig).strip(), si, ei)
    return out


def coord_name(contig, start, end):
    """gapseq 风格基因名：<contig>_<start>:<end>（朝向保留）。"""
    return f"{contig}_{start}:{end}"


GENE_TOKEN = re.compile(r"[A-Za-z0-9_.:'\-]+")


def relabel(model, pos):
    """把模型里能查到坐标的基因名换掉；返回 old->new 映射。

    cobra 把 GPR 存成一个**字符串缓存**（`Reaction._gpr`，经 `gene_reaction_rule`
    属性暴露）。改 `gene.id` 不会更新它。若不显式重写，写出 SBML 时 geneProduct
    列表用的是新名、而 `fbc:geneProductAssociation` 仍是旧名 → 读回来时基因数
    翻倍，且 GPR 指向已不存在的基因（`cobrar::readSBMLmod` 直接 validObject 失败）。

    因此：先改名，再按 old->new 重写每条 GPR 字符串并重新链接。
    """
    mapping = {}
    used = {}
    for g in model.genes:
        old = g.id
        p = pos.get(old)
        if p is None:
            continue
        new = coord_name(*p)
        if new in used:
            # 两个 locus 落到同一坐标区间：保持唯一，加后缀（评测按坐标匹配，
            # 后缀不影响 .GS 解析出的 start:end）
            used[new] += 1
            new = f"{new}_dup{used[new]}"
        else:
            used[new] = 0
        mapping[old] = new

    # 1. 改名
    for g in model.genes:
        if g.id in mapping:
            g.id = mapping[g.id]

    # cobra 的 DictList 内部维护 id->列表下标 的索引；原地改 id 后必须重建，
    # 否则 gene_reaction_rule setter 走 get_by_id() 会拿到错误对象。
    model.genes._dict = {g.id: i for i, g in enumerate(model.genes)}

    # 2. 重写 GPR 字符串并重新链接
    for rxn in model.reactions:
        old_rule = rxn.gene_reaction_rule or ""
        if not old_rule:
            continue
        new_rule = GENE_TOKEN.sub(
            lambda mo: mapping.get(mo.group(0), mo.group(0)), old_rule)
        if new_rule == old_rule:
            continue
        rxn.gene_reaction_rule = new_rule   # 公开 setter 会按新名重新链接基因

    # 3. 自检
    referenced = set()
    for rxn in model.reactions:
        referenced.update(g.id for g in rxn.genes)
    declared = {g.id for g in model.genes}
    missing = sorted(referenced - declared)
    print(f"GPR 引用的基因: {len(referenced)}；model.genes: {len(declared)}；"
          f"引用但未声明: {len(missing)}"
          f"{(' ' + str(missing[:3])) if missing else ''}")
    n_unmapped = sum(1 for g in model.genes if g.id not in referenced)
    print(f"未被任何 GPR 引用的基因（多为去掉的 Unknown/伪基因）: {n_unmapped}")
    return mapping


FORMULA_RE = re.compile(r"^([A-Z][a-z]?[0-9]*)+$")


def sanitize_formulas(model):
    """把不合法的 fbc:chemicalFormula 清成合法字符串。

    CarveMe 会写出 `HPO4;HO4P` 这类分号分隔的备选式；libSBML 的 fbc 校验器
    （`cobrar::readSBMLmod` 走的是严格校验）会直接报错退出。取第一个备选；
    仍不合法就清空（空式在 cobrar 里是允许的）。
    """
    n_fixed = n_cleared = 0
    for met in model.metabolites:
        f = (met.formula or "").strip()
        if not f:
            continue
        first = f.split(";")[0].strip()
        if FORMULA_RE.match(first):
            if first != f:
                n_fixed += 1
            met.formula = first
        else:
            n_cleared += 1
            met.formula = ""
    if n_fixed or n_cleared:
        print(f"化学式清理: 取首个备选 {n_fixed} 个，清空非法 {n_cleared} 个")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model_xml", required=True)
    ap.add_argument("--kind", required=True, choices=["modelseed", "carveme"])
    ap.add_argument("--pos_table", required=True)
    ap.add_argument("--out_xml", required=True)
    ap.add_argument("--out_genepos", default="",
                    help="可选：输出 gene,contig,start,end（gene 为重标注后的名字）")
    ap.add_argument("--out_genepos_orig", default="",
                    help="可选：输出 gene,contig,start,end（gene 为**原始** locus tag）。"
                         "供评测端以同一张坐标表驱动 before/after 两行，避免依赖 "
                         "predGeneEssentiality.MS/.CM 里写死的表文件名")
    ap.add_argument("--out_map", default="",
                    help="可选：输出 old_gene,new_gene 对照表")
    ap.add_argument("--genepos_only", action="store_true",
                    help="只导出坐标表，不读/不写模型。此时 --model_xml/--out_xml "
                         "可给占位值；--out_genepos_orig 会拿到**全量**位置表"
                         "（与官方 prepare_MS_gene_table 的数据源一致，评测端再按 "
                         "@allGenes 过滤）")
    args = ap.parse_args()

    if args.kind == "modelseed":
        pos = load_positions_modelseed(args.pos_table)
    else:
        pos = load_positions_carveme(args.pos_table)
    print(f"位置表: {args.pos_table}  ({len(pos)} 个基因)")

    if args.genepos_only:
        if args.out_genepos_orig:
            with open(args.out_genepos_orig, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["gene", "contig", "start", "end"])
                for g, (contig, s, e) in pos.items():
                    w.writerow([g, contig, s, e])
            print(f"写出 gene_pos_table（全量原始 ID）: {args.out_genepos_orig}")
        if args.out_genepos:
            with open(args.out_genepos, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["gene", "contig", "start", "end"])
                for g, (contig, s, e) in pos.items():
                    w.writerow([coord_name(contig, s, e), contig, s, e])
            print(f"写出 gene_pos_table（全量坐标名）: {args.out_genepos}")
        return

    model = cobra.io.read_sbml_model(args.model_xml)
    n0 = len(model.genes)
    sanitize_formulas(model)
    mapping = relabel(model, pos)
    print(f"模型基因: {n0}  成功重标注: {len(mapping)}  "
          f"({100.0 * len(mapping) / max(1, n0):.1f}%)")
    if len(mapping) < n0:
        miss = [g.id for g in model.genes if g.id not in mapping.values()][:5]
        print(f"  未映射（保留原名，评测时不可定位）示例: {miss}")

    cobra.io.write_sbml_model(model, args.out_xml)
    print(f"写出: {args.out_xml}")

    # new name -> coords（供 --gene_pos_table 使用）
    if args.out_genepos:
        with open(args.out_genepos, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene", "contig", "start", "end"])
            for old, new in mapping.items():
                contig, s, e = pos[old]
                w.writerow([new, contig, s, e])
        print(f"写出 gene_pos_table: {args.out_genepos}")

    # 同一张坐标表、但 gene 列保留原始 locus tag —— 评测端 before 行的输入。
    # 只包含「模型里确实出现过的基因」，与 prepare_MS_gene_table 的
    # `mod.gene %in% mod@allGenes` 过滤等价（这里提前滤掉，评测端也会再滤一次）。
    if args.out_genepos_orig:
        with open(args.out_genepos_orig, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene", "contig", "start", "end"])
            for old in mapping:
                contig, s, e = pos[old]
                w.writerow([old, contig, s, e])
        print(f"写出 gene_pos_table（原始 ID）: {args.out_genepos_orig}")

    if args.out_map:
        with open(args.out_map, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["old_gene", "new_gene"])
            for old, new in sorted(mapping.items()):
                w.writerow([old, new])
        print(f"写出对照表: {args.out_map}")


if __name__ == "__main__":
    main()
