import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT / "src"
GPR_DIR = SRC / "gemagent" / "impl" / "agent" / "agent_gpr_check"
sys.path.insert(0, str(GPR_DIR))

from audit_or2_gpr_tsv import parse_pure_or  # noqa: E402
from review_or2_with_llm import validate_decisions  # noqa: E402


FBC_XML = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2">
  <model id="m">
    <listOfReactions>
      <reaction id="r1">
        <fbc:geneProductAssociation>
          <fbc:or>
            <fbc:geneProductRef fbc:geneProduct="g1"/>
            <fbc:geneProductRef fbc:geneProduct="g2"/>
          </fbc:or>
        </fbc:geneProductAssociation>
      </reaction>
      <reaction id="r2">
        <fbc:geneProductAssociation>
          <fbc:or>
            <fbc:geneProductRef fbc:geneProduct="g3"/>
            <fbc:geneProductRef fbc:geneProduct="g4"/>
          </fbc:or>
        </fbc:geneProductAssociation>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
"""


class OR2ReviewTests(unittest.TestCase):
    def test_all_remove_is_downgraded_to_uncertain(self):
        decisions = validate_decisions(
            [
                {"gene": "g1", "verdict": "REMOVE"},
                {"gene": "g2", "verdict": "REMOVE"},
            ],
            ["g1", "g2"],
        )
        self.assertEqual([d["verdict"] for d in decisions], ["UNCERTAIN", "UNCERTAIN"])

    def test_omitted_branch_is_uncertain(self):
        decisions = validate_decisions(
            [{"gene": "g1", "verdict": "KEEP"}],
            ["g1", "g2"],
        )
        self.assertEqual(decisions[1]["verdict"], "UNCERTAIN")
        self.assertIn("omitted", decisions[1]["reason"])

    def test_parse_pure_or_rejects_nested_rules(self):
        self.assertEqual(parse_pure_or("g1 or g2"), ["g1", "g2"])
        self.assertIsNone(parse_pure_or("(g1 or g2) and g3"))
        self.assertIsNone(parse_pure_or("g1 and g2"))

    def test_empty_cohort_does_not_require_llm_credentials(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            audit = root / "audit"
            out = root / "review"
            config = root / "config"
            audit.mkdir()
            config.mkdir()
            (audit / "or2_reactions.tsv").write_text(
                "reaction_id\tcategory\tn_genes\n", encoding="utf-8"
            )
            (audit / "or2_genes.tsv").write_text(
                "reaction_id\tgene\n", encoding="utf-8"
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "reaction_id\tname\tgpr\tec\n", encoding="utf-8"
            )
            kofam = root / "kofam.tsv"
            kofam.write_text("", encoding="utf-8")
            env = os.environ.copy()
            env["LLM_ENV_FILE"] = str(root / "missing.env")
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            result = subprocess.run(
                [
                    sys.executable,
                    str(GPR_DIR / "review_or2_with_llm.py"),
                    "--audit_dir",
                    str(audit),
                    "--legacy_manifest",
                    str(manifest),
                    "--kofam",
                    str(kofam),
                    "--out_dir",
                    str(out),
                    "--organism",
                    "test",
                    "--cohort",
                    "all",
                    "--limit",
                    "0",
                    "--workers",
                    "1",
                    "--llm_config_dir",
                    str(config),
                ],
                capture_output=True,
                text=True,
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            decision_rows = list(
                csv.DictReader(
                    (out / "llm_decisions.tsv").open(encoding="utf-8"),
                    delimiter="\t",
                )
            )
            self.assertEqual(decision_rows, [])
            manifest_data = json.loads(
                (out / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest_data["n_selected_reactions"], 0)
            self.assertEqual(manifest_data["llm_model"], "")

    def test_selected_reaction_requires_llm_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            audit = root / "audit"
            out = root / "review"
            config = root / "config"
            audit.mkdir()
            config.mkdir()
            (audit / "or2_reactions.tsv").write_text(
                "reaction_id\tcategory\tn_genes\tclusters\n"
                "r1\tambiguous_needs_review\t2\t1.1.1.1:2\n",
                encoding="utf-8",
            )
            (audit / "or2_genes.tsv").write_text(
                "reaction_id\tgene\tgene_label\tkofam_gene\tkos\tecs\n"
                "r1\tg1\tg1\t\t\t1.1.1.1\n"
                "r1\tg2\tg2\t\t\t1.1.1.1\n",
                encoding="utf-8",
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "reaction_id\tname\tgpr\tec\n"
                "r1\ttest\tg1 or g2\t1.1.1.1\n",
                encoding="utf-8",
            )
            kofam = root / "kofam.tsv"
            kofam.write_text("", encoding="utf-8")
            env = os.environ.copy()
            env["LLM_ENV_FILE"] = str(root / "missing.env")
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            result = subprocess.run(
                [
                    sys.executable,
                    str(GPR_DIR / "review_or2_with_llm.py"),
                    "--audit_dir",
                    str(audit),
                    "--legacy_manifest",
                    str(manifest),
                    "--kofam",
                    str(kofam),
                    "--out_dir",
                    str(out),
                    "--organism",
                    "test",
                    "--cohort",
                    "all",
                    "--limit",
                    "0",
                    "--workers",
                    "1",
                    "--llm_config_dir",
                    str(config),
                ],
                capture_output=True,
                text=True,
                env=env,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("LLM configuration error", result.stderr)


class OR2FBCApplierTests(unittest.TestCase):
    def test_remove_only_and_never_empty(self):
        decisions = (
            "reaction_id\tgene\tverdict\n"
            "r1\tg1\tKEEP\n"
            "r1\tg2\tREMOVE\n"
            "r2\tg3\tREMOVE\n"
            "r2\tg4\tREMOVE\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            model = root / "model.xml"
            dec = root / "decisions.tsv"
            out = root / "out.xml"
            model.write_text(FBC_XML, encoding="utf-8")
            dec.write_text(decisions, encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(GPR_DIR / "apply_or2_llm_xml.py"),
                    "--model",
                    str(model),
                    "--decisions",
                    str(dec),
                    "--out",
                    str(out),
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            text = out.read_text(encoding="utf-8")
            self.assertIn('fbc:geneProduct="g1"', text)
            self.assertNotIn('fbc:geneProduct="g2"', text)
            self.assertIn('fbc:geneProduct="g3"', text)
            self.assertIn('fbc:geneProduct="g4"', text)
            manifest = json.loads(
                (root / "out_apply_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["reactions_changed"], 1)
            self.assertEqual(manifest["gene_links_removed"], 1)
            self.assertEqual(manifest["skipped_all_removed"], 1)


if __name__ == "__main__":
    unittest.main()
