import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

from intervaltree import IntervalTree

import _pathfix  # noqa: F401
from ampclasslib.ac_annotation import amplicon_annotation
from ampclasslib.ac_io import bfbarchitect_profile_values, classification_profile_header
from amplicon_classifier import (
    INVALID_CLASS,
    NO_FSCNA_CLASS,
    infer_fan_decomposition_class,
)
from make_results_table import discover_feature_specs


def class_row(decomp, ecdna="None detected", bfb="None detected", chromo="None detected", ec_count="0"):
    return {
        "amplicon_decomposition_class": decomp,
        "ecDNA+": ecdna,
        "BFB+": bfb,
        "FAN+": chromo,
        "ecDNA_amplicons": ec_count,
    }


class ResultsReportingTests(unittest.TestCase):
    def test_default_profile_reports_bfbarchitect_score_and_source(self):
        args = SimpleNamespace(verbose_classification=False, bfbarchitect=True)
        header = classification_profile_header(args, ["raw_score"])
        values = bfbarchitect_profile_values(args, {
            "min_score": "1.09",
            "bfb_sources": "AC|BFBArchitect",
        })

        self.assertEqual(header[-2:], ["BFBArchitect_min_score", "BFB_source"])
        self.assertEqual(values, ["1.09", "AC|BFBArchitect"])

    def test_default_profile_reports_native_source_when_bfbarchitect_is_disabled(self):
        args = SimpleNamespace(verbose_classification=False, bfbarchitect=False)
        header = classification_profile_header(args, ["raw_score"])
        values = bfbarchitect_profile_values(args, {"bfb_sources": "AC"})

        self.assertEqual(header[-2:], ["BFBArchitect_min_score", "BFB_source"])
        self.assertEqual(values, ["NA", "AC"])

    def test_verbose_profile_preserves_extended_bfbarchitect_diagnostics(self):
        args = SimpleNamespace(verbose_classification=True, bfbarchitect=True)
        header = classification_profile_header(args, ["raw_score"])
        values = bfbarchitect_profile_values(args, {
            "min_score": "1.09",
            "passing_region_count": "1",
            "multiplicities": "1",
            "regions": "chr9:1-2:1.09:1",
            "whole_graph_used": "False",
            "reverse_polarity_used": "False",
            "bfb_sources": "BFBArchitect",
        })

        self.assertEqual(header, [
            "sample_name", "amplicon_number", "amplicon_decomposition_class", "ecDNA+", "BFB+", "FAN+",
            "ecDNA_amplicons", "contains_viral", "raw_score", "BFBArchitect_min_score",
            "BFBArchitect_passing_region_count", "BFBArchitect_multiplicities", "BFBArchitect_regions",
            "BFBArchitect_whole_graph_used", "BFBArchitect_reverse_polarity_used", "BFB_source",
        ])
        self.assertEqual(values, [
            "1.09", "1", "1", "chr9:1-2:1.09:1", "False", "False", "BFBArchitect",
        ])

    def test_fan_prevents_cyclic_fallback_result_row(self):
        specs = discover_feature_specs(
            class_row("Cyclic", chromo="Positive"),
            "sample_amplicon1",
            "/tmp/class_beds",
        )

        self.assertEqual(
            [(feature_id, feature) for feature_id, feature, _ in specs],
            [("sample_amplicon1_FAN_1", "FAN")],
        )

    def test_bare_cyclic_is_not_reported_as_feature(self):
        with self.assertLogs(level="WARNING") as logs:
            specs = discover_feature_specs(class_row("Cyclic"), "sample_amplicon1", "/tmp/class_beds")

        self.assertEqual(specs, [])
        self.assertIn("not reporting a Cyclic feature row", "\n".join(logs.output))

    def test_virus_is_reported_as_mechanism_feature(self):
        specs = discover_feature_specs(class_row("Virus"), "sample_amplicon1", "/tmp/class_beds")

        self.assertEqual(
            [(feature_id, feature) for feature_id, feature, _ in specs],
            [("sample_amplicon1_Virus_1", "Virus")],
        )

    def test_linear_and_complex_non_cyclic_remain_fallback_features(self):
        for fallback_class in ("Linear", "Complex-non-cyclic"):
            specs = discover_feature_specs(class_row(fallback_class), "sample_amplicon1", "/tmp/class_beds")
            self.assertEqual(
                [(feature_id, feature) for feature_id, feature, _ in specs],
                [("sample_amplicon1_{}_1".format(fallback_class), fallback_class)],
            )

    def test_annotation_does_not_create_cyclic_feature_bed_source(self):
        with tempfile.TemporaryDirectory(prefix="ac-reporting-") as tmpdir:
            graph = Path(tmpdir) / "sample_amplicon1_graph.txt"
            graph.write_text("sequence chr1:100+ chr1:200- 10\n", encoding="utf-8")

            with self.assertLogs(level="WARNING") as logs:
                _, _, _, feature_dict, _, amp_class = amplicon_annotation(
                    cycleList=[[1]],
                    segSeqD={0: ("", 0, 0), 1: ("chr1", 100, 200)},
                    bfb_cycle_inds=[],
                    ecIndexClusters=[],
                    invalidInds=[],
                    bfbStat=False,
                    ecStat=False,
                    ampClass="Cyclic",
                    graphf=str(graph),
                    add_chr_tag=False,
                    lcD=defaultdict(IntervalTree),
                    ref="GRCh38",
                )

        self.assertEqual(amp_class, "Invalid")
        self.assertNotIn("Cyclic_1", feature_dict)
        self.assertIn("Marking as Invalid", "\n".join(logs.output))

    def test_annotation_keeps_fan_without_cyclic_feature(self):
        with tempfile.TemporaryDirectory(prefix="ac-reporting-") as tmpdir:
            graph = Path(tmpdir) / "sample_amplicon1_graph.txt"
            graph.write_text("sequence chr1:100+ chr1:200- 10\n", encoding="utf-8")

            _, _, _, feature_dict, _, amp_class = amplicon_annotation(
                cycleList=[[1]],
                segSeqD={0: ("", 0, 0), 1: ("chr1", 100, 200)},
                bfb_cycle_inds=[],
                ecIndexClusters=[],
                invalidInds=[],
                bfbStat=False,
                ecStat=False,
                ampClass="Cyclic",
                graphf=str(graph),
                add_chr_tag=False,
                lcD=defaultdict(IntervalTree),
                ref="GRCh38",
                fan_intervals={"chr1": [(100, 200)]},
            )

        self.assertEqual(amp_class, "Cyclic")
        self.assertNotIn("Cyclic_1", feature_dict)
        self.assertIn("FAN_1", feature_dict)

    def test_annotation_keeps_fan_decomposition_when_fallback_intervals_are_filtered(self):
        with tempfile.TemporaryDirectory(prefix="ac-reporting-") as tmpdir:
            graph = Path(tmpdir) / "sample_amplicon1_graph.txt"
            graph.write_text("sequence chr1:100+ chr1:2200- 3\n", encoding="utf-8")

            _, _, _, feature_dict, _, amp_class = amplicon_annotation(
                cycleList=[[1]],
                segSeqD={0: ("", 0, 0), 1: ("chr1", 100, 2200)},
                bfb_cycle_inds=[],
                ecIndexClusters=[],
                invalidInds=[],
                bfbStat=False,
                ecStat=False,
                ampClass="Complex-non-cyclic",
                graphf=str(graph),
                add_chr_tag=False,
                lcD=defaultdict(IntervalTree),
                ref="GRCh38",
                fan_intervals={"chr1": [(100, 2200)]},
            )

        self.assertEqual(amp_class, "Complex-non-cyclic")
        self.assertIn("FAN_1", feature_dict)

    def test_fan_no_fscna_gets_structural_decomposition_class(self):
        amp_profile = defaultdict(float)
        amp_profile.update({
            "No amp/Invalid": 1.0,
            "Linear": 0.0,
            "Trivial cycle": 0.0,
            "Complex-non-cyclic": 0.0,
            "Complex-cyclic": 0.0,
        })

        self.assertEqual(
            infer_fan_decomposition_class(NO_FSCNA_CLASS, amp_profile, rearr_e=2),
            "Complex-non-cyclic",
        )
        self.assertEqual(
            infer_fan_decomposition_class(INVALID_CLASS, amp_profile, rearr_e=2),
            INVALID_CLASS,
        )


if __name__ == "__main__":
    unittest.main()
