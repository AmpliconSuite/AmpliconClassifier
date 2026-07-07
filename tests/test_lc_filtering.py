import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from intervaltree import IntervalTree

from ampclasslib.ac_util import ConfigVars, parseCycle


class LowComplexityCycleFilterTests(unittest.TestCase):
    def setUp(self):
        self.old_bp_fraction = ConfigVars.lc_cycle_max_bp_fraction
        self.old_breakend_fraction = ConfigVars.lc_cycle_max_breakend_fraction
        ConfigVars.lc_cycle_max_bp_fraction = 0.10
        ConfigVars.lc_cycle_max_breakend_fraction = 0.10

    def tearDown(self):
        ConfigVars.lc_cycle_max_bp_fraction = self.old_bp_fraction
        ConfigVars.lc_cycle_max_breakend_fraction = self.old_breakend_fraction

    def _write_graph_and_cycles(self, tmpdir, cycle_segments, discordant_edges=None):
        graph = Path(tmpdir) / "amp_graph.txt"
        cycles = Path(tmpdir) / "amp_cycles.txt"

        with graph.open("w", encoding="utf-8") as out:
            for ind in range(4):
                start = ind * 100
                end = start + 100
                out.write("sequence\tchr1:{}-\tchr1:{}+\t10\t100\t100\n".format(start, end))
            for v1, v2 in discordant_edges or []:
                out.write("discordant\t{}->{}\t1\t10\n".format(v1, v2))

        with cycles.open("w", encoding="utf-8") as out:
            for ind in range(4):
                start = ind * 100
                end = start + 100
                out.write("Segment\t{}\tchr1\t{}\t{}\n".format(ind + 1, start, end))
            out.write(
                "Cycle=1;Copy_count=5;Segments={};Path_constraints_satisfied=\n".format(cycle_segments)
            )

        return str(graph), str(cycles)

    def test_cycle_with_low_lc_bp_fraction_is_retained(self):
        lcD = defaultdict(IntervalTree)
        lcD["chr1"].addi(210, 230)

        with tempfile.TemporaryDirectory(prefix="ac-lc-filter-") as tmpdir:
            graph, cycles = self._write_graph_and_cycles(tmpdir, "1+,2+,3+,4+")

            _segSeqD, cycleList, cycleCNs, lc_filtered = parseCycle(cycles, graph, False, lcD, [])

        self.assertEqual(lc_filtered, 0)
        self.assertEqual(len(cycleList), 1)
        self.assertEqual(cycleCNs, [5.0])

    def test_cycle_with_high_lc_bp_fraction_is_filtered(self):
        lcD = defaultdict(IntervalTree)
        lcD["chr1"].addi(200, 250)

        with tempfile.TemporaryDirectory(prefix="ac-lc-filter-") as tmpdir:
            graph, cycles = self._write_graph_and_cycles(tmpdir, "1+,2+,3+,4+")

            _segSeqD, cycleList, _cycleCNs, lc_filtered = parseCycle(cycles, graph, False, lcD, [])

        self.assertEqual(lc_filtered, 1)
        self.assertEqual(cycleList, [])

    def test_cycle_with_high_lc_breakend_fraction_is_filtered(self):
        lcD = defaultdict(IntervalTree)
        lcD["chr1"].addi(100, 101)

        with tempfile.TemporaryDirectory(prefix="ac-lc-filter-") as tmpdir:
            graph, cycles = self._write_graph_and_cycles(
                tmpdir,
                "1+,3+,2+,4+",
                discordant_edges=[("chr1:100+", "chr1:200-")],
            )

            _segSeqD, cycleList, _cycleCNs, lc_filtered = parseCycle(cycles, graph, False, lcD, [])

        self.assertEqual(lc_filtered, 1)
        self.assertEqual(cycleList, [])


if __name__ == "__main__":
    unittest.main()
