import unittest
from collections import defaultdict

from ampclasslib.classification_records import Feature
from amplicon_classifier import (
    BFBARCHITECT_SOURCE,
    AC_SOURCE,
    cycle_indices_overlapping_bfb_intervals,
    interval_dict_from_intervals,
    merge_bfb_features,
    snap_bfbarchitect_intervals,
)


class BFBFeatureTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
