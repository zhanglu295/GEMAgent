#!/usr/bin/env bash
# =============================================================================
# run_biomass_check.sh — biomass 逐前体 LLM 审查流水线（三步）
#
#   1) biomass_evidence.py      收集每个前体的证据（FBA 可生产性 / 培养基 / 生产反应）
#   2) biomass_check_prompt.py  建 prompt（复用 BASE_BIOLOGY_CONTEXT 作 system）
#   3) review_biomass_check.py  LLM 逐组判定 KEEP/EXCLUDE/UNCERTAIN
#   4) apply_biomass_decisions.py 生成裁剪模板 + 经 replace_biomass.R 重归一出新 bio1
#      （--evidence 必传：硬闸 role 判定与 REPLACE 候选校验都依赖它）
#
# 用法:
#   bash run_biomass_check.sh \
#       --output_dir <pipeline_out> --genome_id ecol \
#       --org_name "Escherichia coli" \
#       --model_xml  <...>-draft.xml \
#       --draft_rds  <...>-draft.RDS \
#       --biomass_json <database>/biomass/biomass_Gram_neg.json \
#       [--medium_csv <...>-medium.csv] [--steps 1,2,3,4] [--dry_run]
#       [--kofam <...>_kofam.tsv] [--all_groups]
#   --kofam       KOfamScan 注释表；缺省自动查找 <output_dir>/kofam/<id>_kofam.tsv
#   --all_groups  连不含 AUDIT 前体的组也出 prompt（生产默认只出 AUDIT 组）
# =============================================================================
set -euo pipefail

# 生产入口是 gemagent.py（它按步骤调用下面 4 个脚本）；本脚本只是单独调试/
# 重跑 biomass review 用的便利封装，因此要求调用前已 activate gemagent 环境。
# 解释器与共享库都从当前环境推导，不写死绝对路径。
PY="${PY:-$(command -v python)}"
ENV_ROOT="$(cd "$(dirname "$PY")/.." && pwd)"
# 本机 libstdc++ 过旧，需优先用 conda env 自带的库（否则 cobra 的 sqlite3 导入失败）
if [ -d "$ENV_ROOT/lib" ]; then
    export LD_LIBRARY_PATH="$ENV_ROOT/lib:${LD_LIBRARY_PATH:-}"
fi
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

OUT_DIR=""; GENOME_ID=""; ORG_NAME=""; MODEL_XML=""; DRAFT_RDS=""
BM_JSON=""; MEDIUM=""; STEPS="1,2,3,4"; DRY=""; WORKERS=4
KOFAM=""; AUDIT_FLAG="--audit_only"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output_dir) OUT_DIR="$2"; shift 2;;
    --genome_id)  GENOME_ID="$2"; shift 2;;
    --org_name)   ORG_NAME="$2"; shift 2;;
    --model_xml)  MODEL_XML="$2"; shift 2;;
    --draft_rds)  DRAFT_RDS="$2"; shift 2;;
    --biomass_json) BM_JSON="$2"; shift 2;;
    --medium_csv) MEDIUM="$2"; shift 2;;
    --steps)      STEPS="$2"; shift 2;;
    --workers)    WORKERS="$2"; shift 2;;
    --kofam)      KOFAM="$2"; shift 2;;
    --all_groups) AUDIT_FLAG=""; shift;;
    --dry_run)    DRY="--dry_run"; shift;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

[[ -n "$OUT_DIR" && -n "$GENOME_ID" ]] || { echo "需要 --output_dir 与 --genome_id" >&2; exit 2; }
CK="$OUT_DIR/biomass_check"
mkdir -p "$CK"
echo "=================================================================="
echo " biomass check : $GENOME_ID"
echo " out_dir       : $CK"
echo " model_xml     : $MODEL_XML"
echo " draft_rds     : $DRAFT_RDS"
echo " biomass_json  : $BM_JSON"
echo " medium_csv    : $MEDIUM"
echo " kofam         : $KOFAM"
echo " audit_only    : ${AUDIT_FLAG:-no}"
echo " steps         : $STEPS"
echo "=================================================================="

has() { [[ ",$STEPS," == *",$1,"* ]]; }

# KOfamScan 注释表（供 related_kos 基因证据）。未显式指定时按生产命名约定自动查找。
if [[ -z "$KOFAM" ]]; then
  for c in "$OUT_DIR/kofam/${GENOME_ID}_kofam.tsv" \
           "$OUT_DIR/${GENOME_ID}_kofam.tsv" \
           "$(dirname "$OUT_DIR")/${GENOME_ID}_kofam/${GENOME_ID}_kofam.tsv"; do
    [[ -f "$c" ]] && { KOFAM="$c"; break; }
  done
fi
[[ -n "$KOFAM" ]] || echo "[note] 未找到 ${GENOME_ID}_kofam.tsv -> related_kos 为空（证据弱于生产路径）"

if has 1; then
  echo; echo "── [1/4] 证据收集 ──"
  $PY "$HERE/biomass_evidence.py" \
      --model_xml "$MODEL_XML" \
      --biomass_json "$BM_JSON" --medium_csv "$MEDIUM" \
      --out_dir "$CK" \
      --kofam "$KOFAM" \
      --template_dir "$(dirname "$BM_JSON")"
fi

if has 2; then
  echo; echo "── [2/4] 构建 prompt ──"
  $PY "$HERE/biomass_check_prompt.py" \
      --evidence "$CK/biomass_evidence.json" \
      --output_dir "$OUT_DIR" --out_dir "$CK" --org_name "$ORG_NAME" \
      $AUDIT_FLAG
fi

if has 3; then
  echo; echo "── [3/4] LLM 审查 ──"
  $PY "$HERE/review_biomass_check.py" \
      --prompt_dir "$CK/prompt" --out_dir "$CK" \
      --evidence "$CK/biomass_evidence.json" --workers "$WORKERS"
fi

if has 4; then
  echo; echo "── [4/4] 应用判定 → 裁剪模板 + 重归一出新 bio1 ──"
  $PY "$HERE/apply_biomass_decisions.py" \
      --draft_rds "$DRAFT_RDS" --biomass_json "$BM_JSON" \
      --decisions "$CK/biomass_decisions.json" --out_dir "$CK" \
      --evidence "$CK/biomass_evidence.json" $DRY
fi

echo; echo "完成。产物:"; ls -la "$CK" | sed 's/^/  /'
