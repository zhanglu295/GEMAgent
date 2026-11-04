import json
import tempfile
import unittest
from pathlib import Path

from gemagent.manifest import sha256, write_manifest


class TestManifest(unittest.TestCase):
    def test_manifest_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            f = root / "input.txt"
            f.write_text("abc", encoding="utf-8")
            self.assertEqual(len(sha256(f)), 64)
            out = write_manifest(root, status="completed", genome_id="x", outputs=[str(f)])
            data = json.loads(out.read_text(encoding="utf-8"))
            self.assertEqual(data["status"], "completed")
            self.assertEqual(data["genome_id"], "x")


if __name__ == "__main__":
    unittest.main()