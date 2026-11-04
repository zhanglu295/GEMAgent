import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from gemagent.impl.agent.agent_gpr_check.gpr_rewrite import apply_actions, parse_isozyme  # noqa: E402


class ParenthesizedTests(unittest.TestCase):
    def test_parenthesized_pure_or(self):
        self.assertEqual(parse_isozyme("(g1) or (g2)"), ["g1", "g2"])

    def test_nested_rejected(self):
        self.assertIsNone(parse_isozyme("(g1 or g2) and g3"))

    def test_parenthesized_rewrite(self):
        new, removed, added, skipped = apply_actions("(g1) or (g2)", {"g1": "REMOVE"})
        self.assertEqual(new, "g2")
        self.assertEqual(removed, ["g1"])


if __name__ == "__main__":
    unittest.main()
