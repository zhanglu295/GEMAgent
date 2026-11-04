import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import cobra
import pandas as pd

from gemagent.impl.agent import agent_us_prune as prune
from gemagent.impl.agent import agent_us_expand as expand


def tiny_model(path, include_a_transport=True):
    model = cobra.Model("tiny")
    a_e = cobra.Metabolite("a_e0", compartment="e0")
    a_c = cobra.Metabolite("a_c0", compartment="c0")
    b_c = cobra.Metabolite("b_c0", compartment="c0")
    c_c = cobra.Metabolite("c_c0", compartment="c0")
    b_e = cobra.Metabolite("b_e0", compartment="e0")
    c_e = cobra.Metabolite("c_e0", compartment="e0")
    bio = cobra.Metabolite("bio_c0", compartment="c0")
    model.add_metabolites([a_e, a_c, b_c, c_c, b_e, c_e, bio])

    ex_a = cobra.Reaction("EX_a_e0")
    ex_a.bounds = (-1000, 1000)
    ex_a.add_metabolites({a_e: -1})
    model.add_reactions([ex_a])
    if include_a_transport:
        ta = cobra.Reaction("rxn_transport_a")
        ta.bounds = (-1000, 1000)
        ta.gene_reaction_rule = "geneA"
        ta.add_metabolites({a_e: -1, a_c: 1})
        model.add_reactions([ta])
    use = cobra.Reaction("rxn_use_a")
    use.bounds = (0, 1000)
    use.add_metabolites({a_c: -1})
    model.add_reactions([use])
    for rid, met in (("EX_b_e0", "b"), ("EX_c_e0", "c")):
        ex = cobra.Reaction(rid)
        ex.bounds = (-1000, 1000)
        ex.add_metabolites({model.metabolites.get_by_id(f"{met}_e0"): -1})
        model.add_reactions([ex])
        tr = cobra.Reaction(f"rxn_transport_{met}")
        tr.bounds = (-1000, 1000)
        tr.add_metabolites({model.metabolites.get_by_id(f"{met}_e0"): -1,
                            model.metabolites.get_by_id(f"{met}_c0"): 1})
        model.add_reactions([tr])
    conv = cobra.Reaction("rxn_bc")
    conv.bounds = (-1000, 1000)
    conv.add_metabolites({b_c: -1, c_c: 1})
    model.add_reactions([conv])
    model.objective = "rxn_use_a"
    cobra.io.write_sbml_model(model, str(path))
    return model


class TestPruneEvidence(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        model_path = self.root / "tiny-gapfilled-gprchecked-usreviewed.xml"
        full = tiny_model(model_path)
        draft = full.copy()
        draft.remove_reactions(["EX_b_e0", "rxn_transport_b", "EX_c_e0", "rxn_transport_c"])
        cobra.io.write_sbml_model(draft, str(self.root / "tiny-draft.xml"))
        cobra.io.write_sbml_model(full, str(self.root / "tiny-gapfilled.xml"))
        cobra.io.write_sbml_model(full, str(self.root / "tiny-gapfilled-gprchecked.xml"))

        (self.root / "tiny-gapfill-added.json").write_text(json.dumps([
            {"seed": "rxn_transport_b", "pathway": "|PW_B|", "pathway_status": "full"}
        ]), encoding="utf-8")
        (self.root / "tiny-all-Pathways.tbl").write_text("""ID\tName\tPrediction\tCompleteness\tVagueReactions\tKeyReactions\tKeyReactionsFound\tReactionsFound
|PW_A|\tA utilization\ttrue\t100\t0\t1\t1\trxn_use_a
|PW_B|\tB utilization\ttrue\t100\t0\t1\t1\trxn_transport_b
""", encoding="utf-8")
        (self.root / "tiny-all-Reactions.tbl").write_text("""rxn\tpathway\tstatus\tpathway.status\tpident\tbitscore\tqcovs
rxn_transport_a\t|PW_A|\tfull\tfull\t57\t400\t100
rxn_use_a\t|PW_A|\tfull\tfull\t45\t220\t95
""", encoding="utf-8")
        (self.root / "tiny-Transporter.tbl").write_text("""rea\tpident\tbitscore\tqcovs
rxn_transport_a\t57\t400\t100
rxn_transport_b\t55\t380\t99
""", encoding="utf-8")
        (self.root / "gpr_check").mkdir()
        (self.root / "gpr_check" / "decisions.tsv").write_text("""reaction_id\tgene\tverdict
rxn_transport_a\tgeneA\tKEEP
rxn_transport_b\tgeneB\tKEEP
""", encoding="utf-8")
        (self.root / "kofam").mkdir()
        (self.root / "kofam" / "tiny_kofam.tsv").write_text("""#\tgene name\tKO\tthrshld\tscore\tE-value\t"KO definition"
*\tgeneA\tK00001\t10\t100\t1e-20\t"A transporter"
""", encoding="utf-8")
        (self.root / "claims.tsv").write_text("""compound\tpolarity\taction\tevidence_form\tsubject_level\tattribution\tgapseq_exchange\tquote\tdoi
A\tpositive\tuptake\tphenotype_test\tstrain\tstrain\tEX_a_e0\tThe strain utilizes A.\t10.1/x
B\tnegative\tuptake\tphenotype_test\tstrain\tstrain\tEX_b_e0\tThe strain does not utilize B.\t10.1/y
""", encoding="utf-8")
        (self.root / "gtdb.tsv").write_text(
            "d__Bacteria;p__Firmicutes;c__Clostridia;o__Eubacteriales;f__Tiny;g__Tiny;s__Tiny testus\n",
            encoding="utf-8")
        (self.root / "mimedb.tsv").write_text("""gram\toxygen_requirement\tmetabolism\tenvironment
Positive\tAnaerobe\tSugar fermenter\tHost
""", encoding="utf-8")
        self.args = SimpleNamespace(
            workspace=str(self.root), genome_id="tiny", model=str(model_path),
            draft=str(self.root / "tiny-draft.xml"),
            gapfilled=str(self.root / "tiny-gapfilled.xml"),
            gprchecked=str(self.root / "tiny-gapfilled-gprchecked.xml"),
            gapfill_added=str(self.root / "tiny-gapfill-added.json"),
            pathways=str(self.root / "tiny-all-Pathways.tbl"),
            reactions=str(self.root / "tiny-all-Reactions.tbl"),
            transporter=str(self.root / "tiny-Transporter.tbl"),
            gpr_decisions=str(self.root / "gpr_check" / "decisions.tsv"),
            kofam=str(self.root / "kofam" / "tiny_kofam.tsv"),
            claims=str(self.root / "claims.tsv"),
            gtdb_summary=str(self.root / "gtdb.tsv"),
            mimedb_tsv=str(self.root / "mimedb.tsv"), org_name="Tiny testus")

    def tearDown(self):
        self.tmp.cleanup()

    def test_features_are_truth_free_and_trace_four_generations(self):
        _, features, gtdb, mimedb = prune.build_features(self.args)
        self.assertNotIn("observed", features.columns)
        self.assertIn("Tiny testus", gtdb)
        self.assertEqual(mimedb["gram"], "Positive")
        text = prune.prompt(features, "Tiny testus", gtdb, mimedb)
        self.assertIn("Tiny testus", text)
        self.assertIn("Anaerobe", text)
        self.assertNotIn("observed", text)
        row = features[(features.exchange_id == "EX_a_e0") &
                       (features.direction == "uptake")].iloc[0]
        self.assertTrue(row["present_after_gapfill"])
        self.assertTrue(row["positive_claim"])
        self.assertEqual(row["triage"], "PROTECTED")
        self.assertIn("PW_A", row["pathway_ids"])
        self.assertGreaterEqual(row["max_pident"], 30)

    def test_write_rds_symbol_is_available(self):
        self.assertTrue(callable(prune.write_rds))

    def test_direct_gtdb_and_single_row_mimedb(self):
        self.assertIn("d__Bacteria", prune.read_gtdb(
            str(self.root / "gtdb.tsv"), "tiny", str(self.root)))
        profile = prune.read_mimedb(str(self.root / "mimedb.tsv"), "unmatched")
        self.assertEqual(profile["oxygen_requirement"], "Anaerobe")

    def test_gtdb_lineage_column_and_optional_missing_inputs(self):
        tax = self.root / "taxonomy.tsv"
        tax.write_text(
            "genome_id\tgtdb_lineage\torg_name\n"
            "tiny\td__Bacteria;p__Firmicutes;g__Tiny;s__Tiny testus\tTiny testus\n",
            encoding="utf-8")
        self.assertIn("d__Bacteria", prune.read_gtdb(str(tax), "tiny", str(self.root)))
        self.assertEqual(prune.read_kofam("."), {})
        self.assertEqual(prune.read_gpr("."), {})
        self.assertTrue(prune.read_table(".").empty)

    def test_block_guard_rolls_back_one_item_and_keeps_another(self):
        features = pd.DataFrame([
            {"exchange_id": "EX_a_e0", "direction": "uptake", "positive_claim": False},
            {"exchange_id": "EX_b_e0", "direction": "uptake", "positive_claim": True},
            {"exchange_id": "EX_c_e0", "direction": "uptake", "positive_claim": False},
        ])
        decisions = [
            {"exchange_id": "EX_a_e0", "direction": "uptake", "decision": "BLOCK"},
            {"exchange_id": "EX_b_e0", "direction": "uptake", "decision": "BLOCK"},
            {"exchange_id": "EX_c_e0", "direction": "uptake", "decision": "BLOCK"},
        ]
        out = self.root / "pruned.xml"
        args = SimpleNamespace(model=str(self.root / "tiny-gapfilled-gprchecked-usreviewed.xml"),
                               out_xml=str(out))
        _, rows, _ = prune.apply_decisions(args, features, decisions)
        by_id = {r["exchange_id"]: r for r in rows}
        self.assertFalse(by_id["EX_a_e0"]["applied"])
        self.assertEqual(by_id["EX_a_e0"]["guard"], "growth_guard")
        self.assertFalse(by_id["EX_b_e0"]["applied"])
        self.assertEqual(by_id["EX_b_e0"]["guard"], "positive_claim_guard")
        self.assertTrue(by_id["EX_c_e0"]["applied"])
        self.assertTrue(out.exists())


class TestExpandWithoutPanel(unittest.TestCase):
    def test_blank_map_cell_does_not_become_nan(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            map_path = root / "exchange_map.tsv"
            map_path.write_text(
                "exchange\tabbr\tcompound_name\tcpd\tmatched_by\tgapseq_exchange\n"
                "EX_cellul(e)\tcellul\tCellulose\t\t\t\n", encoding="utf-8")
            crosswalk = expand.load_crosswalk(str(map_path))
            records = {"EX_cellul_e0": {"reaction_id": "EX_cellul_e0"}}
            ids, _ = expand.resolve_exchange_ids(
                "Cellulose", records, crosswalk, expand.aliases_dict(crosswalk))
            self.assertEqual(ids, ["EX_cellul_e0"])

    def test_direct_claims_emit_add_and_create_without_truth_panel(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            model = root / "tiny.xml"
            tiny_model(model, include_a_transport=False)
            exchange_map = root / "exchange_map.tsv"
            exchange_map.write_text(
                "exchange\tabbr\tcompound_name\tcpd\tmatched_by\tgapseq_exchange\n"
                "EX_a(e)\ta\tA\t\t\t\n"
                "EX_d(e)\td\tD\t\t\t\n", encoding="utf-8")
            claims = root / "claims.tsv"
            claims.write_text("""compound\tpolarity\taction\tevidence_form\tsubject_level\tattribution\tgapseq_exchange\tquote\tdoi
A\tpositive\tuptake\tphenotype_test\tstrain\tstrain\tEX_a_e0\tThe strain utilizes A as a carbon source.\t10.1/x
D\tpositive\tuptake\tphenotype_test\tstrain\tstrain\t\tThe strain utilizes D as a carbon source.\t10.1/z
""", encoding="utf-8")
            out = root / "out"
            xml = root / "expanded.xml"
            subprocess.run([
                sys.executable, str(Path(expand.__file__)), "--model", str(model),
                "--organism", "tiny", "--base", "test", "--claims", str(claims),
                "--exchange-map", str(exchange_map), "--allow-create",
                "--out-dir", str(out), "--out-xml", str(xml),
            ], check=True)
            decisions = pd.read_csv(out / "decisions.tsv", sep="\t")
            self.assertIn("ADD", set(decisions.action))
            self.assertIn("CREATE", set(decisions.action))
            self.assertNotIn("BLOCK", set(decisions.action))
            self.assertFalse(any(c.endswith("FN") or c.endswith("FP")
                                 for c in decisions.columns))
            self.assertTrue(xml.exists())

            loaded = cobra.io.read_sbml_model(str(xml))
            self.assertTrue(any(r.id.startswith("ADDTRAN_") for r in loaded.reactions))
            self.assertTrue(loaded.reactions.has_id("EX_d_e0"))


if __name__ == "__main__":
    unittest.main()
