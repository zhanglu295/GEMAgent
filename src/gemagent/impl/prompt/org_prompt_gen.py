import os
import re
import sys
import argparse
import pandas as pd

# <root>/code 提供共享的 paths；本文件在 <root>/code/prompt/ 下。
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths  # noqa: E402

# 导入刚才单独创建的提示词模块
import agent_prompts

def load_data(file_path=None):
    if file_path is None:
        file_path = paths.db("mimedb", "mimedb_microbes_v2.csv")
    if not os.path.exists(file_path):
        # 降级处理防止本地运行时路径报错
        print(f"Warning: Database path {file_path} not found. Using empty data.")
        return pd.DataFrame()
    return pd.read_csv(file_path)

def clean_value(val):
    """Safely convert non-string/null values for prompt generation."""
    if pd.isna(val) or val is None or str(val).lower() == 'null':
        return "Unknown"
    return str(val).strip()

def consolidate_strains(results_df, note=""):
    """Deduplicates strains while preserving the essence and taxonomic diversity."""
    if results_df.empty:
        return "No comprehensive profile available."

    # 在副本上工作：调用方传进来的可能是 df 的切片，直接赋值会触发
    # SettingWithCopyWarning（每跑一次管线都会刷屏）。
    results_df = results_df.copy()

    # Identify the representative (longest description = most information)
    results_df['desc_len'] = results_df['description'].fillna("").str.len()
    representative = results_df.loc[results_df['desc_len'].idxmax()]
    
    # Aggregate unique metadata from all strains
    unique_metabolism = "; ".join(results_df['metabolism'].dropna().unique())
    unique_habitats = "; ".join(results_df['habitat'].dropna().unique())
    gramstaining = clean_value(representative['gram'])
    
    # Build the structured profile
    # strain 字段缺失时 pd.isna 为真，直接插进 f-string 会写出字面 "nan"。
    species_label = clean_value(representative['species'])
    strain_label = clean_value(representative['strain'])
    if strain_label in ("Unknown", "", "nan", "None"):
        provenance = (f"This profile integrates data from "
                      f"{len(results_df)} strain(s) of {species_label}.")
    else:
        provenance = (f"This profile integrates data from {len(results_df)} "
                      f"strain(s), using {strain_label} as the primary "
                      f"representative.")
    if note:
        provenance = f"{note}\n{provenance}"

    prompt = f"""
### Microbial Consolidated Profile: {species_label}
({provenance})

---
### Taxonomic Lineage
- **Species:** {species_label}
- **Kingdom/Phylum:** {clean_value(representative['kingdom'])} / {clean_value(representative['phylum'])}
- **Class/Order/Family/Genus:** {clean_value(representative['klass'])} / {clean_value(representative['order'])} / {clean_value(representative['family'])} / {clean_value(representative['genus'])}

### Aggregated Physiological Traits
- **Metabolism:** {unique_metabolism if unique_metabolism else "Unknown"}
- **Habitat(s):** {unique_habitats if unique_habitats else "Unknown"}
- **Gram Stain:** {gramstaining}
- **Oxygen Requirement:** {clean_value(representative['oxygen_requirement'])}

### Representative Description
{clean_value(representative['description'])}
"""
    return prompt.strip()


def _mimedb_queries(name):
    """Progressive MIMEdb queries, most specific first: given -> species -> genus.

    ``get_mimedb_info`` matches the query with ``str.contains``, so the query
    must be a *substring* of MIMEdb's bare-binomial ``species`` column. Every
    name we hold is richer than that:

        'Clostridium perfringens ATCC 13124'      -> 0 rows (strain qualifier)
        'Clostridium_S acetobutylicum'            -> 0 rows (GTDB decoration)
        '[Eubacterium] rectale ATCC 33656'        -> 0 rows (brackets)
        'Clostridium perfringens'                 -> 4 rows
        'Eubacterium rectale'                     -> 3 rows

    Verified on database/mimedb/mimedb_microbes_v2.csv. Without this ladder the
    lookup silently returns "No database records found" and every review loses
    its organism context.
    """
    out, seen = [], set()

    def push(value):
        value = value.strip()
        # GTDB suffixes: Clostridium_S -> Clostridium, Methanobrevibacter_A -> ...
        value = re.sub(r"_(?:[A-Z]|[A-Z]{2,3}\d*)$", "", value).strip()
        value = re.sub(r"\s+(?:sp|spp)\.?\d*$", "", value, flags=re.IGNORECASE).strip()
        if value and value.lower() not in seen:
            seen.add(value.lower())
            out.append(value)

    cleaned = re.sub(r"\s+", " ", (name or "").replace("[", "").replace("]", "")).strip()
    if not cleaned:
        return out
    push(cleaned)
    words = cleaned.split(" ")
    if len(words) >= 2:
        push(" ".join(words[:2]))
    if words:
        push(words[0])
    return out


def get_mimedb_info(query, dataframe):
    """Main pipeline: Search (progressively) -> Deduplicate -> Format."""
    if dataframe.empty or not query:
        return "No records found (Empty Database or Query)"

    candidates = _mimedb_queries(query)
    for idx, cand in enumerate(candidates):
        # Fuzzy matching on species
        pattern = re.sub(r'\s+', r'\\s+', cand)
        results = dataframe[dataframe['species'].str.contains(pattern, flags=re.IGNORECASE, na=False)]
        if results.empty:
            continue
        note = ""
        if idx:
            note = (f"NOTE: MimeDB has no entry for the exact name '{query}'. "
                    f"This is the closest available match (query '{cand}'), so "
                    f"the traits below may describe a relative, not this strain.")
        return consolidate_strains(results, note=note)

    return f"No database records found in MimeDB for: {query}"

def _render(template, org_name, mimedb_profile):
    """Fill the two placeholders without str.format.

    The templates contain JSON examples full of braces; ``str.format`` would treat
    them as fields and require them escaped. v1 escaped them as ``{{{{``, which
    ``format`` turned back into ``{{`` in the saved file -- so the LLM was shown
    a malformed JSON example. ``replace`` sidesteps the whole problem.
    """
    return (template
            .replace("{org_name}", org_name)
            .replace("{mimedb_profile}", mimedb_profile))


def main():
    parser = argparse.ArgumentParser(
        description="GemAgent: Automated Genome-Scale Metabolic Model Reconstruction",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Directory to save all output files')
    parser.add_argument('--org_name', type=str, default="",
                        help='Organism name')
    
    args = parser.parse_args()
    org_name = args.org_name if args.org_name else "Unknown Microbe"
    
    # 1. 创建输出目录及子 prompt 目录
    prompt_dir = os.path.join(args.output_dir, "prompts")
    os.makedirs(prompt_dir, exist_ok=True)
    
    # 2. 从 MimeDB 加载并提取特定物种的信息
    df = load_data()
    mimedb_profile = get_mimedb_info(org_name, df)
    
    # 3. 渲染 v2 真正消费的四个 system prompt
    #
    # 消费者（详见各读取处）：
    #   pathway_verification_system.txt -> biomass_check_prompt / gpr_check_prompt
    #                                      / review_gapfill_added（三者共用同一份
    #                                      organism context 作为 system prompt）
    #   dr_extraction_system.txt        -> agent_deepresearch/extract.py
    #   dr_selection_system.txt         -> agent_deepresearch/selectRxn.py
    # v1 前端遗留的 candidate_review / paper_extraction / paper_reading / validation
    # 四个 prompt 在 v2 里没有任何消费者（candidate review 阶段已随 v1 移除），
    # mimedb_gramstaining.txt 也从不被读取（Gram 染色实际由
    # step_read_biomass_template 从 draft 的 biomass 反应名反推），故一并删去。
    prompts_to_save = {
        "pathway_verification_system.txt": agent_prompts.PATHWAY_VERIFICATION_SYSTEM,
        "dr_extraction_system.txt": agent_prompts.DR_EXTRACTION_SYSTEM,
        "dr_selection_system.txt": agent_prompts.DR_SELECTION_SYSTEM,
    }
    
    print(f"\n[INFO] Generating and saving specific system prompts for: {org_name}")
    for filename, template in prompts_to_save.items():
        file_path = os.path.join(prompt_dir, filename)
        content = _render(template, org_name, mimedb_profile)
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        print(f" -> Saved: {file_path} ({len(content)} chars)")
        
    print("\n[SUCCESS] All system prompts successfully generated and archived.")

if __name__ == "__main__":
    main()