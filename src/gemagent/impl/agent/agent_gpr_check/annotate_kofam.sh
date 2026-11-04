#!/bin/bash
#
# GEMAgent — 基因功能注释前置工具（prodigal + KOFamScan）
#
# 为 GPR review（Step 11c）生成基因级 KO 功能注释，供 prompt 内联使用。
# KOFamScan 基于 KEGG Orthology 的 HMM 谱比对，不依赖任何参考模型，
# 对未知/新物种同样适用，是 LLM 判断"基因是否误挂"的功能证据来源。
#
# 用法:
#   annotate_kofam.sh <genome.fna> <out_dir> <genome_id>
#
# 输出:
#   <out_dir>/<genome_id>_proteins.faa     prodigal 预测的蛋白序列
#   <out_dir>/<genome_id>_kofam.tsv        KOFamScan 注释（detail-tsv 格式）
#
# 可通过环境变量覆盖工具/DB 路径:
#   PRODIGAL          prodigal 可执行文件
#   EXEC_ANNOTATION   KOFamScan 的 exec_annotation（ruby）
#   KOFAM_PROFILE     合并的 HMM 谱（prokaryote.hal）
#   KOFAM_KO_LIST     ko_list 文件（KO -> threshold -> definition）
#   KOFAM_DIR         profile/ko_list 所在目录（默认 ~/kofam_local）
#   CPU               HMMER 线程数（默认 8）
#
# gemagent.py 会用 code/paths.py 解析出的值注入上面这些变量，因此这些默认值
# 只在单独手动运行时兜底。
#
set -euo pipefail

if [ $# -lt 3 ]; then
    echo "用法: $0 <genome.fna> <out_dir> <genome_id>" >&2
    exit 1
fi

GENOME=$1
OUT_DIR=$2
GENOME_ID=$3

KOFAM_DIR="${KOFAM_DIR:-$HOME/kofam_local}"
PRODIGAL="${PRODIGAL:-$(command -v prodigal || true)}"
EXEC_ANNOTATION="${EXEC_ANNOTATION:-$(command -v exec_annotation || true)}"
KOFAM_PROFILE="${KOFAM_PROFILE:-$KOFAM_DIR/profiles/prokaryote.hal}"
KOFAM_KO_LIST="${KOFAM_KO_LIST:-$KOFAM_DIR/ko_list}"
CPU="${CPU:-8}"

# exec_annotation is a ruby script (#!/usr/bin/env ruby) and needs hmmer
# (hmmsearch/hmmscan) on PATH. Make the script self-contained.
# 注意：EXEC_ANNOTATION 为空时 dirname "" = "."，会把当前目录塞进 PATH。
if [ -n "$EXEC_ANNOTATION" ]; then
    export PATH="$(dirname "$EXEC_ANNOTATION"):$PATH"
fi

mkdir -p "$OUT_DIR"

# prodigal 无法直接读取 gzip 压缩的基因组，先解压到输出目录
if [[ "$GENOME" == *.gz ]]; then
    UNZIPPED="$OUT_DIR/${GENOME_ID}.fna"
    zcat "$GENOME" > "$UNZIPPED"
    GENOME="$UNZIPPED"
fi

FAA="$OUT_DIR/${GENOME_ID}_proteins.faa"
KO_TSV="$OUT_DIR/${GENOME_ID}_kofam.tsv"

for f in "$GENOME" "$KOFAM_PROFILE" "$KOFAM_KO_LIST"; do
    if [ ! -f "$f" ]; then
        echo "错误: 缺少输入文件 $f" >&2
        exit 1
    fi
done
for exe in "$PRODIGAL" "$EXEC_ANNOTATION"; do
    if [ ! -x "$exe" ]; then
        echo "错误: 缺少可执行文件 $exe" >&2
        exit 1
    fi
done

echo "==> [1/2] prodigal 预测蛋白 ($GENOME_ID)"
"$PRODIGAL" -i "$GENOME" -a "$FAA" \
    -o "$OUT_DIR/${GENOME_ID}_prodigal.gff" \
    -p single -q 2>&1 | tee "$OUT_DIR/${GENOME_ID}_prodigal.log"

KOFAM_TMP="$OUT_DIR/kofam_tmp"
rm -rf "$KOFAM_TMP"
mkdir -p "$KOFAM_TMP"
# per-genome kofam tmp dir isolates runs and avoids cross-species residue
# (exec_annotation defaults to ./tmp and Dir.glob("K*") picks up stale files)
echo "==> [2/2] KOFamScan 注释 ($GENOME_ID)"
"$EXEC_ANNOTATION" \
    -o "$KO_TSV" \
    -f detail-tsv \
    --cpu "$CPU" \
    --tmp-dir "$KOFAM_TMP" \
    --profile "$KOFAM_PROFILE" \
    --ko "$KOFAM_KO_LIST" \
    -E 1e-5 \
    "$FAA" 2>&1 | tee "$OUT_DIR/${GENOME_ID}_kofam.log"

echo "==> 完成: $KO_TSV"
grep -c '^\*' "$KO_TSV" | xargs echo "最佳命中基因数:"
