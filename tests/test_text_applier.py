import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "src", "gemagent", "impl", "agent",
                      "agent_gpr_check", "apply_or2_llm_text.py")

XML = '''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core">
  <model id="m">
    <listOfReactions>
      <reaction id="r1" name="one">
        <notes><body><p>GENE_ASSOCIATION: g1 or g2 or g3</p></body></notes>
      </reaction>
      <reaction id="r2" name="two">
        <notes><body><p>GENE_ASSOCIATION: g4 and g5</p></body></notes>
      </reaction>
      <reaction id="r3" name="three">
        <notes><body><p>GENE_ASSOCIATION: g6 or g7</p></body></notes>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
'''

DECISIONS = '''reaction_id\tgene\tverdict\taction\treason\tevidence_ko\tevidence_ec\tsource
r1\tg3\tREMOVE\tREMOVE\t\tK00001\t1.1.1.1\trule
r2\tg4\tREMOVE\tREMOVE\t\tK00002\t1.1.1.1\trule
r3\tg6\tREMOVE\tREMOVE\t\tK00003\t1.1.1.1\trule
r3\tg7\tREMOVE\tREMOVE\t\tK00004\t1.1.1.1\trule
'''


class TextApplierTests(unittest.TestCase):
    def test_rewrite(self):
        with tempfile.TemporaryDirectory() as d:
            model = os.path.join(d, "m.xml")
            dec = os.path.join(d, "decisions.tsv")
            out = os.path.join(d, "out.xml")
            with open(model, "w", encoding="utf-8") as fh:
                fh.write(XML)
            with open(dec, "w", encoding="utf-8") as fh:
                fh.write(DECISIONS)
            res = subprocess.run(
                [sys.executable, SCRIPT, "--model", model, "--decisions", dec, "--out", out],
                capture_output=True, text=True,
            )
            self.assertEqual(res.returncode, 0, res.stderr + res.stdout)
            with open(out, encoding="utf-8") as fh:
                text = fh.read()
            self.assertIn("GENE_ASSOCIATION: g1 or g2", text)
            self.assertIn("GENE_ASSOCIATION: g4 and g5", text)
            self.assertIn("GENE_ASSOCIATION: g6 or g7", text)
            self.assertIn('id="r2"', text)


if __name__ == "__main__":
    unittest.main()
