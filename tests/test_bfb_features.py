import logging
import os
import unittest
from collections import defaultdict
from types import SimpleNamespace
from tempfile import NamedTemporaryFile, TemporaryDirectory

import _pathfix  # noqa: F401
import amplicon_classifier as ac
from ampclasslib.classification_records import Feature
from amplicon_classifier import (
    BFBARCHITECT_SOURCE,
    AC_SOURCE,
    bfbarchitect_whole_graph_preflight,
    cycle_indices_overlapping_bfb_intervals,
    is_foldback_qc_artifact,
    interval_dict_from_intervals,
    merge_bfb_features,
    snap_bfbarchitect_intervals,
)


class BFBFeatureTests(unittest.TestCase):
    def _run_with_fake_bfbarchitect(self, fake_reconstruct):
        with NamedTemporaryFile(mode="w", delete=True) as graph:
            graph.write("sequence\tchr1:100+\tchr1:200-\t6\t100\t60\n")
            graph.write("discordant\tchr1:150+->chr1:160+\t0\t5\t0\n")
            graph.flush()

            had_args = hasattr(ac, "args")
            old_args = getattr(ac, "args", None)
            old_reconstruct = ac.BFBARCHITECT_RECONSTRUCT
            old_centromeres = ac.BFBARCHITECT_CENTROMERES
            had_logger = hasattr(ac, "logger")
            old_logger = getattr(ac, "logger", None)
            try:
                ac.args = SimpleNamespace(no_bfbarchitect=False, bfb_threads=1, verbose_classification=False, o="/tmp/ac")
                ac.BFBARCHITECT_RECONSTRUCT = fake_reconstruct
                ac.BFBARCHITECT_CENTROMERES = {"chr1": 100000000}
                ac.logger = logging.getLogger("test_bfb_features")

                return ac.run_bfbarchitect(graph.name, fb_edges=2, output_prefix="test_amplicon")
            finally:
                if had_args:
                    ac.args = old_args
                else:
                    del ac.args
                ac.BFBARCHITECT_RECONSTRUCT = old_reconstruct
                ac.BFBARCHITECT_CENTROMERES = old_centromeres
                if had_logger:
                    ac.logger = old_logger
                else:
                    del ac.logger

    def test_foldback_qc_artifact_detects_extreme_edge_count(self):
        self.assertTrue(is_foldback_qc_artifact(101, 0.1))
        self.assertTrue(is_foldback_qc_artifact(16, 0.81))
        self.assertFalse(is_foldback_qc_artifact(100, 0.1))
        self.assertFalse(is_foldback_qc_artifact(16, 0.8))

    def test_bfbarchitect_version_uses_module_version(self):
        self.assertEqual(ac.get_bfbarchitect_version(SimpleNamespace(__version__="1.0.1")), "1.0.1")

    def test_bfbarchitect_version_falls_back_when_module_version_is_missing(self):
        self.assertEqual(ac.get_bfbarchitect_version(SimpleNamespace()), "unknown")

    def test_snap_bfbarchitect_interval_to_nearest_graph_segment_endpoints(self):
        seq_edges = [
            {"chrom": "chr1", "start": 100, "end": 200},
            {"chrom": "chr1", "start": 200, "end": 350},
            {"chrom": "chr1", "start": 350, "end": 500},
        ]

        snapped = snap_bfbarchitect_intervals([("chr1", 112, 339)], seq_edges, add_chr_tag=False)

        self.assertEqual(snapped, [("chr1", 100, 350)])

    def test_cycle_indices_use_high_overlap_threshold(self):
        segSeqD = {
            0: (None, 0, 0),
            1: ("chr1", 100, 200),
            2: ("chr1", 200, 300),
            3: ("chr1", 300, 500),
        }
        cycleList = [[1, 2], [2, 3]]
        bfb_intervals = interval_dict_from_intervals([("chr1", 100, 300)])

        hits = cycle_indices_overlapping_bfb_intervals(cycleList, segSeqD, bfb_intervals, 0.95)

        self.assertEqual(hits, {0})

    def test_overlapping_ac_and_bfbarchitect_features_merge(self):
        ac_feature = Feature(
            feature_id="BFB_1",
            feature_type="BFB",
            intervals=defaultdict(list, {"chr1": [(100, 300)]}),
            cycle_indices={0},
            sources={AC_SOURCE},
        )
        bfbarchitect_feature = Feature(
            feature_id="BFB_2",
            feature_type="BFB",
            intervals=defaultdict(list, {"chr1": [(250, 500)]}),
            cycle_indices={1},
            sources={BFBARCHITECT_SOURCE},
            score=1.5,
        )

        merged = merge_bfb_features([ac_feature, bfbarchitect_feature])

        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].feature_id, "BFB_1")
        self.assertEqual(merged[0].intervals["chr1"], [(100, 500)])
        self.assertEqual(merged[0].cycle_indices, {0, 1})
        self.assertEqual(merged[0].sources, {AC_SOURCE, BFBARCHITECT_SOURCE})
        self.assertEqual(merged[0].score, 1.5)

    def test_whole_graph_preflight_skips_complex_weak_foldback_graph(self):
        with NamedTemporaryFile(mode="w", delete=True) as graph:
            for i in range(251):
                graph.write("sequence\tchr1:{}+\tchr1:{}-\t6\t100\t60\n".format(i * 1000, i * 1000 + 999))
            for i in range(60):
                graph.write("discordant\tchr1:{}+->chr1:{}-\t0\t5\t0\n".format(i * 1000, i * 1000 + 100000))
            graph.flush()

            should_run, stats = bfbarchitect_whole_graph_preflight(graph.name)

        self.assertFalse(should_run)
        self.assertEqual(stats["sequence_edges"], 251)
        self.assertEqual(stats["discordant_edges"], 60)
        self.assertEqual(stats["foldback_edges"], 0)

    def test_whole_graph_preflight_skips_multichromosomal_graph(self):
        with NamedTemporaryFile(mode="w", delete=True) as graph:
            for chrom in ["chr1", "chr2", "chr3"]:
                graph.write("sequence\t{}:1+\t{}:300001-\t12\t100\t60\n".format(chrom, chrom))
            for i in range(6):
                graph.write("discordant\tchr1:{}+->chr1:{}+\t0\t5\t0\n".format(i * 1000, i * 1000 + 1000))
            graph.flush()

            should_run, stats = bfbarchitect_whole_graph_preflight(graph.name)

        self.assertFalse(should_run)
        self.assertEqual(stats["large_chrom_count"], 3)

    def test_whole_graph_preflight_allows_complex_foldback_enriched_graph(self):
        with NamedTemporaryFile(mode="w", delete=True) as graph:
            for i in range(251):
                graph.write("sequence\tchr1:{}+\tchr1:{}-\t6\t100\t60\n".format(i * 1000, i * 1000 + 999))
            for i in range(54):
                graph.write("discordant\tchr1:{}+->chr1:{}-\t0\t5\t0\n".format(i * 1000, i * 1000 + 100000))
            for i in range(6):
                graph.write("discordant\tchr1:{}+->chr1:{}+\t0\t5\t0\n".format(i * 1000, i * 1000 + 1000))
            graph.flush()

            should_run, stats = bfbarchitect_whole_graph_preflight(graph.name)

        self.assertTrue(should_run)
        self.assertEqual(stats["foldback_edges"], 6)

    def test_run_bfbarchitect_skips_whole_graph_retry_for_complex_weak_foldback_graph(self):
        calls = []

        def fake_reconstruct(_graphf, centromere_dict=None, solver=None, multiple=False, whole_graph=False,
                             threads=1, min_lp_bound=None, max_graph_segments=None, reverse_polarity=False):
            calls.append((whole_graph, reverse_polarity, min_lp_bound))
            return []

        with NamedTemporaryFile(mode="w", delete=True) as graph:
            for i in range(251):
                graph.write("sequence\tchr1:{}+\tchr1:{}-\t6\t100\t60\n".format(i * 1000, i * 1000 + 999))
            for i in range(60):
                graph.write("discordant\tchr1:{}+->chr1:{}-\t0\t5\t0\n".format(i * 1000, i * 1000 + 100000))
            graph.flush()

            had_args = hasattr(ac, "args")
            old_args = getattr(ac, "args", None)
            old_reconstruct = ac.BFBARCHITECT_RECONSTRUCT
            old_centromeres = ac.BFBARCHITECT_CENTROMERES
            had_logger = hasattr(ac, "logger")
            old_logger = getattr(ac, "logger", None)
            try:
                ac.args = SimpleNamespace(no_bfbarchitect=False, bfb_threads=1, verbose_classification=False, o="/tmp/ac")
                ac.BFBARCHITECT_RECONSTRUCT = fake_reconstruct
                ac.BFBARCHITECT_CENTROMERES = {"chr1": 100000000}
                ac.logger = logging.getLogger("test_bfb_features")

                summary = ac.run_bfbarchitect(graph.name, fb_edges=2, output_prefix="test_amplicon")
            finally:
                if had_args:
                    ac.args = old_args
                else:
                    del ac.args
                ac.BFBARCHITECT_RECONSTRUCT = old_reconstruct
                ac.BFBARCHITECT_CENTROMERES = old_centromeres
                if had_logger:
                    ac.logger = old_logger
                else:
                    del ac.logger

        self.assertEqual(calls, [(False, False, 25), (False, True, 25)])
        self.assertEqual(summary["whole_graph_used"], "False")
        self.assertEqual(summary["regions"], "whole_graph_skipped_complex_weak_foldback")

    def test_run_bfbarchitect_failed_summary_uses_best_whole_graph_reverse_score_mode(self):
        calls = []

        def fake_reconstruct(_graphf, centromere_dict=None, solver=None, multiple=False, whole_graph=False,
                             threads=1, min_lp_bound=None, max_graph_segments=None, reverse_polarity=False,
                             region=None):
            calls.append((whole_graph, reverse_polarity, region))
            if whole_graph and reverse_polarity:
                return [{"scores": [3.1], "region": "whole_graph_reverse", "multiplicity": 1}]
            if whole_graph:
                return [{"scores": [3.7], "region": "whole_graph", "multiplicity": 1}]
            if reverse_polarity:
                return [{"scores": [3.4], "region": region, "multiplicity": 1}]
            return [{"scores": [4.5], "region": ("chr1", 100, 200), "multiplicity": 1}]

        summary = self._run_with_fake_bfbarchitect(fake_reconstruct)

        self.assertEqual(calls, [
            (False, False, None),
            (True, False, None),
            (False, True, ("chr1", 100, 200)),
            (True, True, None),
        ])
        self.assertEqual(summary["min_score"], "3.10")
        self.assertEqual(summary["passing_region_count"], "0")
        self.assertEqual(summary["whole_graph_used"], "True")
        self.assertEqual(summary["reverse_polarity_used"], "True")

    def test_run_bfbarchitect_failed_summary_keeps_region_mode_when_region_score_is_best(self):
        def fake_reconstruct(_graphf, centromere_dict=None, solver=None, multiple=False, whole_graph=False,
                             threads=1, min_lp_bound=None, max_graph_segments=None, reverse_polarity=False,
                             region=None):
            if whole_graph and reverse_polarity:
                return [{"scores": [5.0], "region": "whole_graph_reverse", "multiplicity": 1}]
            if whole_graph:
                return [{"scores": [4.0], "region": "whole_graph", "multiplicity": 1}]
            if reverse_polarity:
                return [{"scores": [4.5], "region": region, "multiplicity": 1}]
            return [{"scores": [3.0], "region": ("chr1", 100, 200), "multiplicity": 1}]

        summary = self._run_with_fake_bfbarchitect(fake_reconstruct)

        self.assertEqual(summary["min_score"], "3.00")
        self.assertEqual(summary["passing_region_count"], "0")
        self.assertEqual(summary["whole_graph_used"], "False")
        self.assertEqual(summary["reverse_polarity_used"], "False")

    def test_run_bfbarchitect_stops_after_passing_region_result(self):
        calls = []

        def fake_reconstruct(_graphf, centromere_dict=None, solver=None, multiple=False, whole_graph=False,
                             threads=1, min_lp_bound=None, max_graph_segments=None, reverse_polarity=False,
                             region=None):
            calls.append((whole_graph, reverse_polarity, region))
            return [{"scores": [2.1], "region": ("chr1", 100, 200), "multiplicity": 1}]

        summary = self._run_with_fake_bfbarchitect(fake_reconstruct)

        self.assertEqual(calls, [(False, False, None)])
        self.assertEqual(summary["min_score"], "2.10")
        self.assertEqual(summary["passing_region_count"], "1")
        self.assertEqual(summary["whole_graph_used"], "False")
        self.assertEqual(summary["reverse_polarity_used"], "False")

    def test_visualize_uses_centromere_path_for_released_bfbarchitect(self):
        calls = []

        def fake_visualize(cycle_file, graph_file, cnr_file, output_prefix, multiple=False, centromere=None):
            calls.append((cycle_file, graph_file, cnr_file, output_prefix, multiple, centromere))

        old_visualize = ac.BFBARCHITECT_VISUALIZE
        old_centromere_path = ac.BFBARCHITECT_CENTROMERE_PATH
        try:
            ac.BFBARCHITECT_VISUALIZE = fake_visualize
            ac.BFBARCHITECT_CENTROMERE_PATH = "/data/GRCh38_centromere.bed"

            ac.run_bfbarchitect_visualize("cycles.txt", "graph.txt", "output")
        finally:
            ac.BFBARCHITECT_VISUALIZE = old_visualize
            ac.BFBARCHITECT_CENTROMERE_PATH = old_centromere_path

        self.assertEqual(calls, [
            ("cycles.txt", "graph.txt", None, "output", False, "/data/GRCh38_centromere.bed")
        ])

    def test_visualize_uses_centromere_dict_for_legacy_bfbarchitect(self):
        calls = []

        def fake_visualize(cycle_file, graph_file, cnr_file, output_prefix, multiple=False,
                           centromere_dict=None):
            calls.append((cycle_file, graph_file, cnr_file, output_prefix, multiple, centromere_dict))

        old_visualize = ac.BFBARCHITECT_VISUALIZE
        old_centromeres = ac.BFBARCHITECT_CENTROMERES
        try:
            ac.BFBARCHITECT_VISUALIZE = fake_visualize
            ac.BFBARCHITECT_CENTROMERES = {"chr1": 100000000}

            ac.run_bfbarchitect_visualize("cycles.txt", "graph.txt", "output")
        finally:
            ac.BFBARCHITECT_VISUALIZE = old_visualize
            ac.BFBARCHITECT_CENTROMERES = old_centromeres

        self.assertEqual(calls, [
            ("cycles.txt", "graph.txt", None, "output", False, {"chr1": 100000000})
        ])

    def test_bfbarchitect_outputs_are_written_without_verbose_classification(self):
        calls = []

        def fake_write_graph(path, new_segments, svs, sv_info):
            calls.append(("graph", path, new_segments, svs, sv_info))

        def fake_write_cycles(path, new_segments, bfb_strings, scores, multiplicity):
            calls.append(("cycles", path, new_segments, bfb_strings, scores, multiplicity))

        def fake_visualize(cycles_file, graph_file, output_prefix):
            calls.append(("visualize", cycles_file, graph_file, output_prefix))

        had_args = hasattr(ac, "args")
        old_args = getattr(ac, "args", None)
        old_write_graph = ac.BFBARCHITECT_WRITE_GRAPH
        old_write_cycles = ac.BFBARCHITECT_WRITE_CYCLES
        old_visualize_impl = ac.BFBARCHITECT_VISUALIZE
        old_visualize = ac.run_bfbarchitect_visualize
        try:
            with TemporaryDirectory() as tmpdir:
                ac.args = SimpleNamespace(
                    verbose_classification=False,
                    o=os.path.join(tmpdir, "classification", "sample")
                )
                ac.BFBARCHITECT_WRITE_GRAPH = fake_write_graph
                ac.BFBARCHITECT_WRITE_CYCLES = fake_write_cycles
                ac.BFBARCHITECT_VISUALIZE = fake_visualize
                ac.run_bfbarchitect_visualize = fake_visualize

                ac.write_bfbarchitect_outputs([{
                    "new_segments": ["segment"],
                    "svs": ["sv"],
                    "sv_info": {"sv": "info"},
                    "bfb_strings": ["+-"],
                    "scores": [1.0],
                    "multiplicity": 1,
                }], "sample_amplicon1", whole_graph_used=False)

                output_prefix = os.path.join(
                    tmpdir, "classification", "bfbarchitect_outputs",
                    "sample_amplicon1_region1_BFB"
                )
                self.assertEqual(calls, [
                    ("graph", output_prefix + "_graph.txt", ["segment"], ["sv"], {"sv": "info"}),
                    ("cycles", output_prefix + "_cycles.txt", ["segment"], ["+-"], [1.0], 1),
                    ("visualize", output_prefix + "_cycles.txt", output_prefix + "_graph.txt", output_prefix),
                ])
        finally:
            if had_args:
                ac.args = old_args
            else:
                del ac.args
            ac.BFBARCHITECT_WRITE_GRAPH = old_write_graph
            ac.BFBARCHITECT_WRITE_CYCLES = old_write_cycles
            ac.BFBARCHITECT_VISUALIZE = old_visualize_impl
            ac.run_bfbarchitect_visualize = old_visualize

    def test_run_bfbarchitect_reconstruct_warns_for_unsupported_expected_kwargs(self):
        calls = []

        def fake_reconstruct(_graphf, centromere_dict=None, solver=None, multiple=False, whole_graph=False, threads=1):
            calls.append(whole_graph)
            return []

        had_args = hasattr(ac, "args")
        old_args = getattr(ac, "args", None)
        old_reconstruct = ac.BFBARCHITECT_RECONSTRUCT
        old_centromeres = ac.BFBARCHITECT_CENTROMERES
        old_warned = set(ac.BFBARCHITECT_UNSUPPORTED_KWARGS_WARNED)
        had_logger = hasattr(ac, "logger")
        old_logger = getattr(ac, "logger", None)
        try:
            ac.args = SimpleNamespace(bfb_threads=1)
            ac.BFBARCHITECT_RECONSTRUCT = fake_reconstruct
            ac.BFBARCHITECT_CENTROMERES = {"chr1": 100000000}
            ac.BFBARCHITECT_UNSUPPORTED_KWARGS_WARNED.clear()
            ac.logger = logging.getLogger("test_bfb_features")

            with self.assertLogs("test_bfb_features", level="WARNING") as logs:
                result = ac.run_bfbarchitect_reconstruct("graph.txt", whole_graph=True)
        finally:
            if had_args:
                ac.args = old_args
            else:
                del ac.args
            ac.BFBARCHITECT_RECONSTRUCT = old_reconstruct
            ac.BFBARCHITECT_CENTROMERES = old_centromeres
            ac.BFBARCHITECT_UNSUPPORTED_KWARGS_WARNED.clear()
            ac.BFBARCHITECT_UNSUPPORTED_KWARGS_WARNED.update(old_warned)
            if had_logger:
                ac.logger = old_logger
            else:
                del ac.logger

        self.assertEqual(result, [])
        self.assertEqual(calls, [True])
        self.assertIn("min_lp_bound", logs.output[0])


if __name__ == "__main__":
    unittest.main()
