import unittest
from gemagent.skills import SKILLS


class TestRegistry(unittest.TestCase):
    def test_registry_contains_review_skills(self):
        self.assertEqual(set(SKILLS), {"biomass", "gapfill-review", "gpr", "us-review",
                                     "us-review-expand", "us-review-prune"})
        self.assertIn("expand", SKILLS["us-review"].description)
        for spec in SKILLS.values():
            self.assertTrue(spec.name)
            self.assertTrue(spec.requires)


if __name__ == "__main__":
    unittest.main()
