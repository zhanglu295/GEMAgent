#!/usr/bin/env python3
"""
CarveMe SBML 的无损重标注（坐标基因名）。

为什么不用 relabel_model.py（cobra 路线）
-----------------------------------------
CarveMe 1.2.2 写出的 bsub_protein.xml 是**纯 notes 式 GPR**（实测全文
`fbc:geneProduct` 出现 0 次、无 `listOfGeneProducts`、无
`fbc:geneProductAssociation`），基因关联只写在每个 <reaction> 的

    <notes><html:p>GENE_ASSOCIATION: G_NP_391348_1</html:p></notes>

里；化学式同样只写在 <species> 的 `<p>FORMULA: ...</p>` 里（无
`fbc:chemicalFormula`）。

cobra 的 read+write 往返会：
  1. 把 notes 里的 `G_<locus>` 拆成 gene 对象，写出时再生成
     `fbc:geneProduct`（id = "G_" + 转义后的 gene.id），
  2. **不更新 notes**，于是同一个反应的 GPR 在文件里出现两份：
     新的 fbc 关联 + 旧的 notes 关联 —— cobrar 读回来基因数从 1100
     涨到 1422（实测）；
  3. 转义 `"."`/`":"` 为 `__46__`/`__58__`，并把 `<p>FORMULA>` 转成
     `fbc:chemicalFormula`（CarveMe 的备选式 `HPO4;HO4P` 会让 cobrar
     的严格校验直接报错，所以 relabel_model.py 还得额外清洗化学式）。
  4. 实测最要命的是**基因丢失**：重标注的 1099 个基因里只有 853 个真正
     带坐标名出现在输出 SBML 中，其余 246 个在 sybil/cobrar 两种读法下
     都查不到 —— 评测端 postable 覆盖率因此只有 853/1099。

本脚本只做**纯文本替换**：把 GENE_ASSOCIATION notes 行里的 `G_<locus>`
逐字换成坐标名（去掉 `G_` 前缀，与评测端 postable 表同名），文件其余字节
（化学式 / 上下界 / objective / notes 结构）完全不变。因此：
  - 1100 个基因全部保真（无损性由 gene_map / postable 覆盖率闸门验证）；
  - sybil 与 cobrar 都在 notes 里读到**字面**坐标名，不再有 `__46__` 转义，
    cobrar 的 `gp_...__DOT__` 还原由 restore_gene_names.R 兜底；
  - cobra 读该 XML（gpr_check_prompt.py）时基因名逐字保留，
    `parse_model_coord()` 直接可用，KOFam 桥不再依赖转义还原。

用法:
  python relabel_carveme_text.py --model_xml bsub_protein.xml \
      --pos_table bsub_genes_pos.tsv \
      --out_xml bsub-relabeled.xml --out_genepos gene_pos_carveme.tsv \
      --out_map gene_map_carveme.tsv
"""
import argparse
import csv
import re
import sys

from relabel_model import load_positions_carveme, coord_name

# notes 里的基因 token：G_<locus>。locus 允许字母数字与 . _ -（CarveMe 的
# NP_391348_1 / YP_... / 少数带点的 RefSeq accession）。
TOKEN = re.compile(r"G_([A-Za-z0-9_.\-]+)")
# 只改 GENE_ASSOCIATION 行（行内到 "<" 为止，即 </html:p> 之前）
LINE = re.compile(r"GENE_ASSOCIATION:[^\n<]*")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model_xml", required=True)
    ap.add_argument("--pos_table", required=True, help="CarveMe *genes_pos.tsv")
    ap.add_argument("--out_xml", required=True)
    ap.add_argument("--out_genepos", default="")
    ap.add_argument("--out_map", default="")
    args = ap.parse_args()

    pos = load_positions_carveme(args.pos_table)
    newname = {g: coord_name(c, s, e) for g, (c, s, e) in pos.items()}
    print(f"位置表: {args.pos_table}  ({len(pos)} 个基因)")

    with open(args.model_xml, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    n_before = len(text)

    hits = {}
    unmapped = set()

    def repl_token(mo):
        old = mo.group(1)
        new = newname.get(old)
        if new is None:
            unmapped.add(old)
            return mo.group(0)          # 查不到就原样保留
        hits[old] = hits.get(old, 0) + 1
        return new

    def repl_line(mo):
        return TOKEN.sub(repl_token, mo.group(0))

    out, n_lines = LINE.subn(repl_line, text)
    print(f"GENE_ASSOCIATION 行: {n_lines}")
    print(f"命中基因 token 出现次数: {sum(hits.values())}（不同基因 {len(hits)}）")

    # 长度校验：输出与输入的长度差必须恰好等于基因名替换的长度差。
    # 这是「没有误改其它字节」的强证据（任何意外的正则吞噬都会破坏它）。
    delta = sum(c * (len(newname[g]) - len("G_" + g)) for g, c in hits.items())
    assert len(out) == len(text) + delta, \
        f"长度校验失败: out={len(out)} expect={len(text) + delta}"
    print(f"  长度校验通过（仅基因名长度变化 {delta:+d} bytes）")

    # 自检：GENE_ASSOCIATION 行里不应再残留「位置表里有坐标」的旧名
    left = set()
    for mo in LINE.finditer(out):
        for tm in TOKEN.finditer(mo.group(0)):
            if tm.group(1) in pos:
                left.add(tm.group(1))
    assert not left, f"仍有未替换的基因: {sorted(left)[:5]}"
    if unmapped:
        print(f"  未映射（保留原名，位置表里没有）: {sorted(unmapped)[:8]}")

    print(f"写出: {args.out_xml}  ({n_before} -> {len(out)} bytes)")
    with open(args.out_xml, "w", encoding="utf-8") as fh:
        fh.write(out)

    # 坐标表：只含模型里真正出现过的基因（与评测端 @allGenes 过滤等价）
    if args.out_genepos:
        with open(args.out_genepos, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene", "contig", "start", "end"])
            for g in sorted(hits):
                c, s, e = pos[g]
                w.writerow([newname[g], c, s, e])
        print(f"写出 gene_pos_table（坐标名，仅模型内基因）: {args.out_genepos}")

    if args.out_map:
        with open(args.out_map, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["old_gene", "new_gene"])
            for g in sorted(hits):
                w.writerow([g, newname[g]])
        print(f"写出对照表: {args.out_map}")

    # 无损性闸门：模型内基因数应等于位置表里与之匹配的数量
    print(f"\n闸门: 模型内被重标注的基因 = {len(hits)} / 位置表 {len(pos)}；"
          f"未映射保留原名 = {len(unmapped)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
