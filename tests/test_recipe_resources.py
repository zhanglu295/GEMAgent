import contextlib
import io
import unittest

from gemagent.cli import main


class TestRecipeResources(unittest.TestCase):
    def test_packaged_recipes_are_listed(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = main(["recipe", "list"])
        self.assertEqual(status, 0)
        text = output.getvalue()
        self.assertIn("full", text)
        self.assertIn("base-only", text)


if __name__ == "__main__":
    unittest.main()
