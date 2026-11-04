"""Namespace adapter tests: one synthetic family per reconstruction tool."""
import unittest

from gemagent import namespace as ns


class TestBaseAliases(unittest.TestCase):
    def test_known_bases(self):
        self.assertEqual(ns.for_base("gapseq").family, "seed")
        self.assertEqual(ns.for_base("modelseed").family, "seed")
        self.assertEqual(ns.for_base("carveme").family, "bigg")
        self.assertEqual(ns.for_base("agora2").family, "agora")
        self.assertEqual(ns.for_base("AGORA2").family, "agora")

    def test_unknown_base_is_loud(self):
        with self.assertRaises(ns.UnknownNamespace):
            ns.for_base("mgen")


class TestSeedNamespace(unittest.TestCase):
    """gapseq and ModelSEED share cpdXXXXX / c0 / e0."""

    def setUp(self):
        self.adapter = ns.for_base("gapseq")

    def test_split_species(self):
        self.assertEqual(self.adapter.split_species("M_cpd00001_e0"), ("cpd00001", "e"))
        self.assertEqual(self.adapter.split_species("M_cpd00001_c0"), ("cpd00001", "c"))
        self.assertEqual(self.adapter.split_species("M_cpd00001_p0"), ("cpd00001", "p"))

    def test_cytosolic_partner(self):
        self.assertEqual(self.adapter.cytosolic_partner("M_cpd00001_e0"),
                         "M_cpd00001_c0")
        self.assertEqual(self.adapter.cytosolic_partner("M_cpd00001_c0"), "")

    def test_exchange_parsing(self):
        self.assertTrue(self.adapter.is_exchange("R_EX_cpd00009_e0"))
        self.assertTrue(self.adapter.is_exchange("EX_cpd00009_e0"))
        self.assertFalse(self.adapter.is_exchange("R_rxn00001_c0"))
        self.assertEqual(self.adapter.exchange_compound("R_EX_cpd00009_e0"), "cpd00009")
        self.assertEqual(self.adapter.new_exchange_id("cpd00009"), "EX_cpd00009_e0")


class TestBiggNamespace(unittest.TestCase):
    """CarveMe uses BiGG ids and a C_ compartment attribute."""

    def setUp(self):
        self.adapter = ns.for_base("carveme")

    def test_split_species(self):
        self.assertEqual(self.adapter.split_species("M_10fthf_c"), ("10fthf", "c"))
        self.assertEqual(self.adapter.split_species("M_10fthf_e"), ("10fthf", "e"))

    def test_cytosolic_partner(self):
        self.assertEqual(self.adapter.cytosolic_partner("M_10fthf_e"), "M_10fthf_c")

    def test_exchange_parsing(self):
        self.assertEqual(self.adapter.exchange_compound("R_EX_2ameph_e"), "2ameph")
        self.assertEqual(self.adapter.new_exchange_id("2ameph"), "EX_2ameph_e")


class TestAgoraNamespace(unittest.TestCase):
    """AGORA2 encodes brackets in both species and exchange ids."""

    def setUp(self):
        self.adapter = ns.for_base("agora2")

    def test_split_species(self):
        self.assertEqual(self.adapter.split_species("M_10fthf__91__c__93__"),
                         ("10fthf", "c"))
        self.assertEqual(self.adapter.split_species("M_10fthf__40__e__41__"),
                         ("10fthf", "e"))

    def test_cytosolic_partner(self):
        self.assertEqual(self.adapter.cytosolic_partner("M_10fthf__40__e__41__"),
                         "M_10fthf__91__c__93__")

    def test_exchange_parsing(self):
        self.assertEqual(
            self.adapter.exchange_compound("R_EX_12ppd_S__40__e__41__"), "12ppd_S")
        self.assertEqual(self.adapter.new_exchange_id("12ppd_S"), "EX_12ppd_S_e")


class TestDetection(unittest.TestCase):
    def test_detect_seed(self):
        text = '<species id="M_cpd00001_c0" compartment="c0"/>'
        self.assertEqual(ns.detect(text).family, "seed")

    def test_detect_bigg(self):
        text = '<species id="M_10fthf_c" compartment="C_c"/>'
        self.assertEqual(ns.detect(text).family, "bigg")

    def test_detect_agora(self):
        text = '<species id="M_10fthf__91__c__93__" compartment="c"/>'
        self.assertEqual(ns.detect(text).family, "agora")

    def test_detect_unknown_is_loud(self):
        with self.assertRaises(ns.UnknownNamespace):
            ns.detect('<species id="x_c" compartment="c"/>')


if __name__ == "__main__":
    unittest.main()
