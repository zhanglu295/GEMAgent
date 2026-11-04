#!/usr/bin/env python3
"""
review_biomass_check.py — 并发调用 LLM 审查每个 biomass 组分组的 prompt。

（pipeline 步骤：biomass check 的第 2 步）

复用项目统一 LLM 客户端（code/agent/llm_config.py 的 llm），并发方式与
review_gpr_check.py 一致：ThreadPoolExecutor + response_format=json_object +
空内容/解析失败重试。每个 prompt 文件自包含（SYSTEM + ORGANISM CONTEXT + 证据表
+ Task），本脚本按第一个 "---" 拆成 system / user 两段。

输入:
  <prompt_dir>/manifest.tsv + <prompt_dir>/<group>.txt
  <out_dir>/biomass_evidence.json   （用于校验每个 cpd 都有判定）
输出:
  <out_dir>/biomass_decisions.json  结构化（含 raw response）
  <out_dir>/biomass_decisions.tsv   扁平表: group, cpd, name, verdict, reason
  并汇总 <out_dir>/biomass_exclude.tsv  仅 EXCLUDE 的 cpd（供 apply 使用）

支持断点续跑：已存在于 biomass_decisions.json 的组会跳过。

用法:
  python3 review_biomass_check.py \
      --prompt_dir <out_dir>/prompt \
      --out_dir    <out_dir> \
      --workers 4
"""
import argparse
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# 本脚本位于 code/agent/agent_biomass_check/，llm_config 在其上一级 code/agent/
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
from llm_config import llm  # noqa: E402

_OUTPUT_CONSTRAINT = (
    "\n\n## FINAL OUTPUT FORMAT (STRICT - overrides the Task section above)\n"
    "Output ONLY a compact JSON object. No markdown fences, no commentary.\n"
    '{"group": "<group>", "decisions": [{"cpd": "<cpdXXXXX>", "name": "<name>", '
    '"verdict": "KEEP|DROP|UNCERTAIN", "reason": "<short>"}]}\n\n'
    "RULES:\n"
    "- One entry per precursor listed in the evidence table. Do NOT skip any.\n"
    "- verdict is exactly one of: KEEP, DROP, UNCERTAIN.\n"
    "- There is NO substitution verdict: a wrong-form template entry is kept.\n"
    "- reason <= 40 words; keep it valid JSON (no unescaped double quotes).\n"
    "- Output ONLY the JSON object."
)

VALID = ("KEEP", "DROP", "UNCERTAIN")


def split_prompt(text):
    marker = "\n---\n"
    i = text.find(marker)
    if i < 0:
        return "", text
    return text[:i].strip(), text[i + len(marker):].strip()


def load_manifest(prompt_dir):
    path = os.path.join(prompt_dir, "manifest.tsv")
    rows = []
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= len(header):
                rows.append(dict(zip(header, p)))
    return rows


def expected_from_evidence(evidence_path):
    """group -> [ {cpd,name}, ... ]，用于补齐 LLM 漏判的条目。

    只收 AUDIT_* 前体：AUTO_KEEP(大分子单体/无机离子/能量货币) 不需要 LLM 判定。
    """
    ev = json.load(open(evidence_path, encoding="utf-8"))
    out = {}
    for p in ev["precursors"]:
        if not str(p.get("audit_tier", "")).startswith("AUDIT"):
            continue
        g = p.get("template_group") or "_not_in_template"
        out.setdefault(g, []).append({"cpd": p["cpd"], "name": p["met_name"]})
    return out


def _regex_extract(text):
    out = []
    for m in re.finditer(
        r'"cpd"\s*:\s*"([^"]+)"[^}]*?"verdict"\s*:\s*"([^"]+)"', text
    ):
        v = m.group(2).strip().upper()
        if v in VALID:
            out.append({"cpd": m.group(1).strip(), "verdict": v, "reason": ""})
    return out


def parse_decisions(text, expected):
    """expected: [{"cpd","name"}]；缺失项补 UNCERTAIN（apply 时按 KEEP 处理）。"""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("empty LLM response")
    cleaned = text.strip()
    if cleaned.startswith("```"):
        cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned)
        cleaned = re.sub(r"\s*```$", "", cleaned)

    data = None
    try:
        data = json.loads(cleaned)
    except Exception:
        repaired = None
        for cut in [m.end() for m in re.finditer(r"\}\s*,", cleaned)][::-1]:
            try:
                obj = json.loads(cleaned[:cut] + "]}")
                if isinstance(obj, dict) and "decisions" in obj:
                    repaired = obj
                    break
            except Exception:
                continue
        if repaired is None and '"decisions"' in cleaned:
            reg = _regex_extract(cleaned)
            if reg:
                data = {"decisions": reg}
        else:
            data = repaired
        if data is None:
            raise ValueError(f"unparseable JSON: {cleaned[:140]}...")

    raw = data.get("decisions", [])
    if not isinstance(raw, list):
        raise ValueError("'decisions' is not a list")

    name_of = {e["cpd"]: e["name"] for e in expected}
    dec = []
    seen = set()
    for it in raw:
        if not isinstance(it, dict):
            continue
        cpd = str(it.get("cpd", "")).strip()
        if not cpd:
            continue
        verdict = str(it.get("verdict", "")).strip().upper()
        if verdict not in VALID:
            verdict = "UNCERTAIN"
        dec.append({"cpd": cpd, "name": it.get("name") or name_of.get(cpd, ""),
                    "verdict": verdict,
                    "reason": (it.get("reason") or "")[:300]})
        seen.add(cpd)
    for e in expected:
        if e["cpd"] not in seen:
            dec.append({"cpd": e["cpd"], "name": e["name"], "verdict": "UNCERTAIN",
                        "reason": "LLM did not decide; keep (conservative)."})
    return dec


def review_one(row, prompt_dir, max_retries, expected_by_group):
    group = row["group"]
    txt = open(os.path.join(prompt_dir, row["prompt_file"]), encoding="utf-8").read()
    system, user = split_prompt(txt)
    msgs = []
    if system:
        msgs.append({"role": "system", "content": system})
    msgs.append({"role": "user", "content": user + _OUTPUT_CONSTRAINT})

    expected = expected_by_group.get(group, [])
    bound = llm.bind(response_format={"type": "json_object"}, max_tokens=4096)
    last = ""
    for attempt in range(max_retries + 1):
        try:
            resp = bound.invoke(msgs)
            last = resp.content or ""
            return group, parse_decisions(last, expected), last
        except Exception as e:
            print(f"  [RETRY {attempt + 1}/{max_retries}] {group}: {e}", flush=True)
            if attempt < max_retries:
                time.sleep(2 ** attempt)
    dec = [{"cpd": e["cpd"], "name": e["name"], "verdict": "UNCERTAIN",
            "reason": "LLM error after retries; keep."} for e in expected]
    return group, dec, last


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prompt_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--evidence", default="", help="默认 <out_dir>/biomass_evidence.json")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max_retries", type=int, default=2)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    ev_path = args.evidence or os.path.join(args.out_dir, "biomass_evidence.json")
    expected_by_group = expected_from_evidence(ev_path)

    rows = load_manifest(args.prompt_dir)
    if args.limit:
        rows = rows[:args.limit]
    print(f"manifest 组数: {len(rows)}", flush=True)

    jpath = os.path.join(args.out_dir, "biomass_decisions.json")
    done = {}
    if os.path.exists(jpath):
        prev = json.load(open(jpath, encoding="utf-8"))
        done = {k: v for k, v in prev.items() if not k.startswith("_")}
        print(f"已有结果: {len(done)} 组（跳过）", flush=True)

    todo = [r for r in rows if r["group"] not in done]
    print(f"待处理: {len(todo)} 组", flush=True)
    if not todo:
        print("全部已处理。", flush=True)
    else:
        results = dict(done)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(review_one, r, args.prompt_dir, args.max_retries,
                                expected_by_group): r for r in todo}
            n = 0
            for fut in as_completed(futs):
                r = futs[fut]
                try:
                    g, dec, raw = fut.result()
                except Exception as e:
                    g, dec, raw = r["group"], [], f"worker error: {e}"
                results[g] = {"n_precursors": r.get("n_precursors", ""),
                              "decisions": dec, "raw_response": raw}
                n += 1
                print(f"  [progress] {n}/{len(todo)}  {g}", flush=True)
    if not todo:
        results = done

    cnt = {"KEEP": 0, "DROP": 0, "UNCERTAIN": 0}
    for g, v in results.items():
        for d in v["decisions"]:
            cnt[d["verdict"]] = cnt.get(d["verdict"], 0) + 1
    results["_summary"] = {"total_groups": len([k for k in results if not k.startswith("_")]),
                           "total_precursors": sum(cnt.values()), "verdict_counts": cnt}

    os.makedirs(args.out_dir, exist_ok=True)
    json.dump(results, open(jpath, "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    with open(os.path.join(args.out_dir, "biomass_decisions.tsv"), "w",
              encoding="utf-8") as f:
        f.write("group\tcpd\tname\tverdict\treason\n")
        for g in sorted(k for k in results if not k.startswith("_")):
            for d in results[g]["decisions"]:
                f.write(f"{g}\t{d['cpd']}\t{d['name']}\t{d['verdict']}\t{d['reason']}\n")

    with open(os.path.join(args.out_dir, "biomass_drop.tsv"), "w", encoding="utf-8") as f:
        f.write("cpd\tname\tgroup\treason\n")
        for g in sorted(k for k in results if not k.startswith("_")):
            for d in results[g]["decisions"]:
                if d["verdict"] == "DROP":
                    f.write(f"{d['cpd']}\t{d['name']}\t{g}\t{d['reason']}\n")

    print(f"\n完成。输出:\n  {jpath}\n  {os.path.join(args.out_dir, 'biomass_decisions.tsv')}"
          f"\n  {os.path.join(args.out_dir, 'biomass_drop.tsv')}")
    print(f"\nverdict 分布: {cnt}")
    for g in sorted(k for k in results if not k.startswith("_")):
        ex = [d["cpd"] for d in results[g]["decisions"] if d["verdict"] == "DROP"]
        if ex:
            print(f"  {g:18s} DROP {len(ex):2d}: {', '.join(ex)}")


if __name__ == "__main__":
    main()
