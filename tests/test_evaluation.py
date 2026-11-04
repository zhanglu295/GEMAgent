import unittest


class TestEvaluation(unittest.TestCase):
    def test_confusion_metrics(self):
        from gemagent.evaluation import confusion_metrics
        m = confusion_metrics(tp=8, fp=2, tn=6, fn=4)
        self.assertAlmostEqual(m["sensitivity"], 8 / 12)
        self.assertAlmostEqual(m["specificity"], 6 / 8)
        self.assertAlmostEqual(m["F2"], 5 * (8 / 10) * (8 / 12) / (4 * (8 / 10) + (8 / 12)))
        self.assertAlmostEqual(m["MCC"], 40 / (10 * 12 * 8 * 10) ** 0.5)

    def test_decision_mapping_applies_only_successful_blocks(self):
        import pandas as pd
        from gemagent.evaluation import decision_to_predicted
        frame = pd.DataFrame([
            {"decision": "BLOCK", "applied": True},
            {"decision": "BLOCK", "applied": False},
            {"decision": "KEEP", "applied": False},
            {"decision": "UNCERTAIN", "applied": False},
        ])
        mapped = decision_to_predicted(frame)
        self.assertEqual(mapped["predicted"].tolist(), [0, 1, 1, 1])


if __name__ == "__main__":
    unittest.main()


