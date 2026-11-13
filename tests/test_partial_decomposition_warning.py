import logging
import tempfile
import unittest
from pathlib import Path

import _pathfix  # noqa: F401
from ampclasslib.ac_io import (
    find_summary_beside_cycles,
    read_partial_decompositions,
    read_summary_paths,
    warn_partial_decompositions,
)
from ampclasslib.classification_records import AmpliconInput


SEP = "-----------------------------------------------"

# CoRAL >= 3.1.2 summary: amplicon1 stopped early by the time limit, amplicon2 finished.
CORAL_PARTIAL_SUMMARY = """CoRAL 3.1.2
1/2 amplicons solved.
Runtime Limit: 0 s
Profiling Enabled: False
{sep}
AmpliconID = 1
#Intervals = 1
Amplicon Intervals:
\tAmplicon 1 > chr6:39,444,248-41,645,280
Total Amplicon Size: 2,201,033
# Chromosomes: 1
# Sequence Edges: 4
# Concordant Edges: 3
# Discordant Edges: 2
# Non-Source Edges: 9
# Source Edges: 4
Cycle Decomposition Status: PARTIAL
\tModelMetadata: max_weight, k=1, alpha=0.01, total_weights=16730947.241071638, resolution=0.1, num_path_constraints=3
\tStopped early: global time limit (0 s) reached after 1 greedy iterations.
\tCN-weighted length explained: 75.8% (greedy normally stops at 90%)
\tPath constraints satisfied: 1/3 (greedy normally stops at 90%)
\tHeaviest graph walk solved was a cycle.
\tCycle=1;Copy_count=6.338641017762292;Segments=2+,3-
{sep}
AmpliconID = 2
#Intervals = 1
Amplicon Intervals:
\tAmplicon 2 > chr8:126,000,000-128,000,000
Total Amplicon Size: 2,000,001
# Chromosomes: 1
# Sequence Edges: 4
# Concordant Edges: 3
# Discordant Edges: 2
# Non-Source Edges: 9
# Source Edges: 4
Cycle Decomposition Status: SUCCESS
\tModelMetadata: max_weight, k=1, alpha=0.01, total_weights=1000.0, resolution=0.1, num_path_constraints=0
\tCN-weighted length explained: 95.4%
\tHeaviest graph walk solved was a cycle.
\tCycle=1;Copy_count=20.0;Segments=1+
{sep}
Solver Settings
Solver: SCIP
Threads: 1
Time Limit: 60 s
""".format(sep=SEP)

# CoRAL <= 3.1.1 summary: no PARTIAL status exists.
CORAL_OLD_SUMMARY = """CoRAL 3.1.1
1/3 amplicons solved.
Runtime Limit: 21600 s
Profiling Enabled: False
{sep}
AmpliconID = 1
#Intervals = 1
Cycle Decomposition Status: FAILURE
\tModelMetadata: max_weight, k=1, alpha=0.01, total_weights=1.0, resolution=0.1, num_path_constraints=0
{sep}
AmpliconID = 2
#Intervals = 1
Cycle Decomposition Status: SUCCESS
\tHeaviest graph walk solved was a cycle.
{sep}
AmpliconID = 3
#Intervals = 1
Cycle Decomposition Status: UNPERFORMED
""".format(sep=SEP)

# AmpliconArchitect summary (first two amplicons of the COLO320DM cell line run).
AA_SUMMARY = """#Amplicons = 2
-----------------------------------------------------------------------------------------
[amplicon1] AmpliconID = 1
[amplicon1] #Intervals = 1
[amplicon1] Intervals = chr1:146509560-148014548
[amplicon1] OncogenesAmplified = BCL9,CHD1L,
[amplicon1] TotalIntervalSize = 1504989
[amplicon1] AmplifiedIntervalSize = 1504984
[amplicon1] AverageAmplifiedCopyCount = 5.004494911219381
[amplicon1] #Chromosomes = 1
[amplicon1] #SequenceEdges = 5
[amplicon1] #BreakpointEdges = 1
[amplicon1] #CoverageShifts = 2
[amplicon1] #MeanshiftSegmentsCopyCount>5 = 0
[amplicon1] #Foldbacks = 1
[amplicon1] #CoverageShiftsWithBreakpointEdges = 2
[amplicon1] -----------------------------------------------------------------------------------------
[amplicon2] AmpliconID = 2
[amplicon2] #Intervals = 9
[amplicon2] OncogenesAmplified = MCL1,ARNT,
[amplicon2] AverageAmplifiedCopyCount = 10.028518816397451
[amplicon2] #Chromosomes = 7
[amplicon2] -----------------------------------------------------------------------------------------
"""


def amp_input(sample, ampN, root):
    return AmpliconInput(
        "{}_{}".format(sample, ampN), sample, ampN,
        str(root / "{}_{}_cycles.txt".format(sample, ampN)),
        str(root / "{}_{}_graph.txt".format(sample, ampN)),
    )


class ReadPartialDecompositionsTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="ac-partial-")
        self.root = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return str(path)

    def test_coral_partial_amplicon_is_reported_with_details(self):
        path = self.write("s_summary.txt", CORAL_PARTIAL_SUMMARY)
        self.assertEqual(read_partial_decompositions(path), {
            "amplicon1": [
                "Stopped early: global time limit (0 s) reached after 1 greedy iterations.",
                "CN-weighted length explained: 75.8% (greedy normally stops at 90%)",
                "Path constraints satisfied: 1/3 (greedy normally stops at 90%)",
            ]
        })

    def test_older_coral_summary_has_no_partial(self):
        path = self.write("s_summary.txt", CORAL_OLD_SUMMARY)
        self.assertEqual(read_partial_decompositions(path), {})

    def test_aa_summary_has_no_partial(self):
        path = self.write("s_summary.txt", AA_SUMMARY)
        self.assertEqual(read_partial_decompositions(path), {})

    def test_summary_without_status_lines_or_missing_file(self):
        path = self.write("s_summary.txt", "#Amplicons = 0\nsome other tool output\n")
        self.assertEqual(read_partial_decompositions(path), {})
        self.assertEqual(read_partial_decompositions(str(self.root / "missing_summary.txt")), {})

    def test_read_summary_paths(self):
        path = self.write("AC_summary_map.txt", "sA\t/x/sA_summary.txt\nsB\t\n\n")
        self.assertEqual(read_summary_paths(path), {"sA": "/x/sA_summary.txt"})

    def test_find_summary_beside_cycles(self):
        cycles = self.write("sA_amplicon1_cycles.txt", "")
        self.assertIsNone(find_summary_beside_cycles(cycles, "sA"))
        coral_summary = self.write("sA_amplicon_summary.txt", "")
        self.assertEqual(find_summary_beside_cycles(cycles, "sA"), coral_summary)
        aa_summary = self.write("sA_summary.txt", "")
        self.assertEqual(find_summary_beside_cycles(cycles, "sA"), aa_summary)


class WarnPartialDecompositionsTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="ac-partial-")
        self.root = Path(self._tmp.name)
        (self.root / "coral_summary.txt").write_text(CORAL_PARTIAL_SUMMARY, encoding="utf-8")
        (self.root / "aa_summary.txt").write_text(AA_SUMMARY, encoding="utf-8")

    def tearDown(self):
        self._tmp.cleanup()

    def test_warns_only_for_partial_coral_amplicon(self):
        inputs = [
            amp_input("coral", "amplicon1", self.root),
            amp_input("coral", "amplicon2", self.root),
            amp_input("aa", "amplicon1", self.root),
            amp_input("aa", "amplicon2", self.root),
            amp_input("no_summary", "amplicon1", self.root),
        ]
        summary_paths = {
            "coral": str(self.root / "coral_summary.txt"),
            "aa": str(self.root / "aa_summary.txt"),
        }
        with self.assertLogs("ampclasslib.ac_io", level="WARNING") as logs:
            n_warned = warn_partial_decompositions(inputs, summary_paths)

        self.assertEqual(n_warned, 1)
        self.assertEqual(len(logs.records), 1)
        message = logs.records[0].getMessage()
        self.assertTrue(message.startswith("coral amplicon1: CoRAL cycle decomposition is PARTIAL"), message)
        self.assertIn("global time limit (0 s)", message)
        self.assertIn("CN-weighted length explained: 75.8%", message)

    def test_no_warning_for_aa_only_input(self):
        inputs = [amp_input("aa", "amplicon1", self.root), amp_input("aa", "amplicon2", self.root)]
        # assertNoLogs needs Python 3.10; AC supports 3.8, so log a sentinel and check it is the only record.
        with self.assertLogs("ampclasslib.ac_io", level="WARNING") as logs:
            n_warned = warn_partial_decompositions(inputs, {"aa": str(self.root / "aa_summary.txt")})
            logging.getLogger("ampclasslib.ac_io").warning("sentinel")
        self.assertEqual(n_warned, 0)
        self.assertEqual([r.getMessage() for r in logs.records], ["sentinel"])


if __name__ == "__main__":
    unittest.main()
