import contextlib
import io
import unittest

from gemagent.cli import main


class TestCLI(unittest.TestCase):
    def test_skill_list(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = main(["skill", "list"])
        self.assertEqual(rc, 0)
        text = out.getvalue()
        for name in ("biomass", "gapfill-review", "gpr", "us-review",
                     "us-review-expand", "us-review-prune"):
            self.assertIn(name, text)


if __name__ == "__main__":
    unittest.main()
