import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
MAKE_INPUT = REPO_ROOT / "ampclasslib" / "make_input.sh"


def touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")


class MakeInputTests(unittest.TestCase):
    def test_make_input_pairs_cycles_and_graphs_and_writes_summary_map(self):
        with tempfile.TemporaryDirectory(prefix="ac-make-input-") as tmpdir:
            root = Path(tmpdir)
            aa_dir = root / "sampleA" / "sampleA_AA_results"
            touch(aa_dir / "sampleA_amplicon1_cycles.txt")
            touch(aa_dir / "sampleA_amplicon1_graph.txt")
            touch(aa_dir / "sampleA_amplicon2_graph.txt")
            touch(root / "sampleA" / "sampleA_summary.txt")
            outpre = root / "out" / "collection"

            result = subprocess.run(
                ["bash", str(MAKE_INPUT), str(root), str(outpre)],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Warning: Graph file has no matching cycles file", result.stderr)

            input_lines = (root / "out" / "collection.input").read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(input_lines), 1)
            fields = input_lines[0].split("\t")
            self.assertEqual(fields[0], "sampleA")
            self.assertEqual(fields[1], str(aa_dir / "sampleA_amplicon1_cycles.txt"))
            self.assertEqual(fields[2], str(aa_dir / "sampleA_amplicon1_graph.txt"))

            summary_lines = (root / "out" / "collection_summary_map.txt").read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertEqual(summary_lines, ["sampleA\t{}".format(root / "sampleA" / "sampleA_summary.txt")])

    def test_make_input_fails_when_cycles_file_has_no_matching_graph(self):
        with tempfile.TemporaryDirectory(prefix="ac-make-input-") as tmpdir:
            root = Path(tmpdir)
            aa_dir = root / "sampleB" / "sampleB_AA_results"
            touch(aa_dir / "sampleB_amplicon1_cycles.txt")
            outpre = root / "collection"

            result = subprocess.run(
                ["bash", str(MAKE_INPUT), str(root), str(outpre)],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Missing graph file for cycles file", result.stderr)
            self.assertIn("sampleB_amplicon1_graph.txt", result.stderr)

    def test_make_input_excludes_bfbarchitect_output_directories(self):
        with tempfile.TemporaryDirectory(prefix="ac-make-input-") as tmpdir:
            root = Path(tmpdir)
            aa_dir = root / "sampleC_AA_results"
            touch(aa_dir / "sampleC_amplicon1_cycles.txt")
            touch(aa_dir / "sampleC_amplicon1_graph.txt")

            for dirname in ("bfbarchitect_outputs", "rerun_BFBArchitect_outputs"):
                output_dir = root / dirname
                touch(output_dir / "sampleC_BFB_cycles.txt")
                touch(output_dir / "sampleC_BFB_graph.txt")
                touch(output_dir / "sampleC_summary.txt")

            outpre = root / "collection"
            result = subprocess.run(
                ["bash", str(MAKE_INPUT), str(root), str(outpre)],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            input_lines = (root / "collection.input").read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(input_lines), 1)
            self.assertIn(str(aa_dir / "sampleC_amplicon1_cycles.txt"), input_lines[0])
            self.assertNotIn("BFBArchitect", input_lines[0])
            self.assertNotIn("bfbarchitect", input_lines[0])
            self.assertEqual(
                (root / "collection_summary_map.txt").read_text(encoding="utf-8"),
                "",
            )


if __name__ == "__main__":
    unittest.main()
