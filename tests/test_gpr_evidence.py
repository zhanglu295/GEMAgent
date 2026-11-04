import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from gemagent.impl.agent.agent_gpr_check.gpr_evidence import (  # noqa: E402
    Consensus,
    classify_gene,
    decide_reaction,
    ec_consensus,
    gene_ec_map,
    is_full_ec,
    parse_ecs_from_definition,
)


class ECParsingTests(unittest.TestCase):
    def test_single_ec(self):
        self.assertEqual(parse_ecs_from_definition("A [EC:1.2.3.4]"), {"1.2.3.4"})

    def test_multiple_ec(self):
        self.assertEqual(
            parse_ecs_from_definition("A [EC:1.2.3.4 1.2.3.5]"),
            {"1.2.3.4", "1.2.3.5"},
        )

    def test_partial_ec(self):
        self.assertEqual(parse_ecs_from_definition("A [EC:2.7.1.-]"), {"2.7.1.-"})

    def test_no_ec(self):
        self.assertEqual(parse_ecs_from_definition("hypothetical protein"), set())

    def test_is_full(self):
        self.assertTrue(is_full_ec("1.2.3.4"))
        self.assertFalse(is_full_ec("1.2.3.-"))
        self.assertFalse(is_full_ec("1.2.3"))


class ConsensusTests(unittest.TestCase):
    def test_confident(self):
        ecs = {"g1": {"1.1.1.1"}, "g2": {"1.1.1.1"}, "g3": {"1.1.1.1"}}
        c = ec_consensus(ecs)
        self.assertEqual(c.status, "confident")
        self.assertEqual(c.ec, "1.1.1.1")
        self.assertEqual(c.n_support, 3)
        self.assertAlmostEqual(c.confidence, 1.0)

    def test_single_gene_insufficient(self):
        c = ec_consensus({"g1": {"1.1.1.1"}})
        self.assertEqual(c.status, "insufficient")

    def test_no_ec_insufficient(self):
        c = ec_consensus({"g1": set(), "g2": set()})
        self.assertEqual(c.status, "insufficient")
        self.assertIsNone(c.ec)

    def test_ambiguous_two_clusters(self):
        c = ec_consensus({"g1": {"1.1.1.1"}, "g2": {"2.2.2.2"}})
        self.assertIn(c.status, ("ambiguous", "insufficient"))
        self.assertFalse(c.confident)

    def test_promiscuous_gene_does_not_vote(self):
        c = ec_consensus(
            {
                "g1": {"1.1.1.1"},
                "g2": {"1.1.1.1"},
                "g3": {"2.2.2.2", "3.3.3.3", "4.4.4.4", "5.5.5.5"},
            }
        )
        self.assertEqual(c.status, "confident")
        self.assertEqual(c.ec, "1.1.1.1")


class ClassificationTests(unittest.TestCase):
    def test_states(self):
        c = Consensus(status="confident", ec="1.2.3.4", confidence=1.0, n_support=3)
        self.assertEqual(classify_gene({"1.2.3.4"}, c), "support")
        self.assertEqual(classify_gene({"1.2.3.-"}, c), "partial")
        self.assertEqual(classify_gene({"9.9.9.9"}, c), "conflict")
        self.assertEqual(classify_gene(set(), c), "unknown")

    def test_abstain_never_conflicts(self):
        c = Consensus(status="insufficient", ec=None, confidence=0.0, n_support=0)
        self.assertEqual(classify_gene({"9.9.9.9"}, c), "unknown")


class DecisionTests(unittest.TestCase):
    def test_remove_conflict_keeps_supporters(self):
        genes = ["g1", "g2", "g3"]
        ecs = {"g1": {"1.1.1.1"}, "g2": {"1.1.1.1"}, "g3": {"2.2.2.2"}}
        d = decide_reaction("r1", genes, ecs)
        self.assertEqual(d.removals, ["g3"])
        self.assertEqual(d.additions, [])
        self.assertFalse(d.needs_llm)

    def test_unknown_genes_are_kept(self):
        genes = ["g1", "g2", "g3"]
        ecs = {"g1": {"1.1.1.1"}, "g2": {"1.1.1.1"}, "g3": set()}
        d = decide_reaction("r1", genes, ecs)
        self.assertEqual(d.removals, [])
        self.assertEqual([x.action for x in d.decisions], ["KEEP", "KEEP", "KEEP"])

    def test_never_empty(self):
        genes = ["g1", "g2"]
        ecs = {"g1": {"1.1.1.1"}, "g2": {"2.2.2.2"}}
        d = decide_reaction("r1", genes, ecs)
        self.assertEqual(d.removals, [])
        self.assertFalse(all(x.action == "REMOVE" for x in d.decisions))

    def test_mixed_annotation_routes_to_llm(self):
        genes = ["g1", "g2"]
        ecs = {"g1": {"1.1.1.1"}, "g2": set()}
        d = decide_reaction("r1", genes, ecs)
        self.assertEqual(d.removals, [])
        self.assertTrue(d.needs_llm)

    def test_all_unknown_stays_noop(self):
        d = decide_reaction("r1", ["g1", "g2"], {"g1": set(), "g2": set()})
        self.assertEqual(d.removals, [])
        self.assertFalse(d.needs_llm)

    def test_gene_ec_map_union(self):
        mapping = gene_ec_map({
            "g1": [("K1", "a [EC:1.1.1.1]"), ("K2", "b [EC:2.2.2.2]")],
            "g2": [("K3", "c")],
        })
        self.assertEqual(mapping["g1"], {"1.1.1.1", "2.2.2.2"})
        self.assertEqual(mapping["g2"], set())


if __name__ == "__main__":
    unittest.main()
