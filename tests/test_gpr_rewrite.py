import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from gemagent.impl.agent.agent_gpr_check.gpr_rewrite import apply_actions, parse_isozyme  # noqa: E402


class RewriteTests(unittest.TestCase):
    def test_parse_pure_or(self):
        self.assertEqual(parse_isozyme("a or b"), ["a", "b"])
        self.assertIsNone(parse_isozyme("a and b"))
        self.assertIsNone(parse_isozyme("(a or b) and c"))

    def test_remove_keeps_order(self):
        new, removed, added, skipped = apply_actions("g1 or g2 or g3", {"g2": "REMOVE"})
        self.assertEqual(new, "g1 or g3")
        self.assertEqual(removed, ["g2"])
        self.assertEqual(added, [])
        self.assertFalse(skipped)

    def test_add_appends(self):
        new, removed, added, skipped = apply_actions("g1 or g2", {"g3": "ADD"})
        self.assertEqual(new, "g1 or g2 or g3")
        self.assertEqual(added, ["g3"])
        self.assertFalse(skipped)

    def test_never_empty(self):
        new, removed, added, skipped = apply_actions("g1 or g2", {"g1": "REMOVE", "g2": "REMOVE"})
        self.assertEqual(new, "g1 or g2")
        self.assertEqual(removed, [])
        self.assertTrue(skipped)

    def test_keep_and_uncertain_ignored(self):
        new, removed, added, skipped = apply_actions("g1 or g2", {"g1": "KEEP", "g2": "UNCERTAIN"})
        self.assertEqual(new, "g1 or g2")

    def test_non_isozyme_untouched(self):
        new, removed, added, skipped = apply_actions("g1 and g2", {"g1": "REMOVE"})
        self.assertEqual(new, "g1 and g2")
        self.assertEqual(removed, [])


if __name__ == "__main__":
    unittest.main()
