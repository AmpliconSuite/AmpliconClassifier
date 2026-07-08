import unittest

import _pathfix  # noqa: F401
import amplicon_classifier as ac


def registry_entry(sample_name, amplicon_number, feature_type, feature_number="1"):
    return {
        "sample_name": sample_name,
        "amplicon_number": amplicon_number,
        "feature_type": feature_type,
        "feature_number": feature_number,
    }


def similarity_row(feature_id_1, feature_id_2, pvalue=1.0, asym1=0.0, asym2=0.0):
    return [feature_id_1, feature_id_2, 1.0, 0.0, pvalue, asym1, asym2]


class SimilarityFilteringTests(unittest.TestCase):
    def test_same_sample_pairs_are_not_generated(self):
        feature_registry = {
            "sampleA_amplicon1_ecDNA_1": registry_entry("sampleA", "amplicon1", "ecDNA"),
            "sampleA_amplicon2_ecDNA_1": registry_entry("sampleA", "amplicon2", "ecDNA"),
            "sampleB_amplicon1_ecDNA_1": registry_entry("sampleB", "amplicon1", "ecDNA"),
        }
        feat_to_ivald = {
            feature_id: (None, None)
            for feature_id in feature_registry
        }

        pairs = list(ac.get_cross_sample_feature_pairs(feat_to_ivald, feature_registry, lambda _a, _b: True))
        self.assertNotIn(("sampleA_amplicon2_ecDNA_1", "sampleA_amplicon1_ecDNA_1"), pairs)
        self.assertEqual(len(pairs), 2)

    def test_fan_filters_amplicon_scope(self):
        feature_registry = {
            "sampleA_amplicon1_FAN_1": registry_entry("sampleA", "amplicon1", "FAN"),
            "sampleB_amplicon1_ecDNA_1": registry_entry("sampleB", "amplicon1", "ecDNA"),
        }
        rows = [similarity_row("sampleA_amplicon1_FAN_1", "sampleB_amplicon1_ecDNA_1", pvalue=1e-10)]

        feats_to_filter, amps_to_filter, filter_events = ac.choose_similarity_filter_targets(
            rows, feature_registry, pval=1e-5
        )

        self.assertEqual(feats_to_filter, {("sampleB", "amplicon1", "ecDNA", "1")})
        self.assertEqual(amps_to_filter, {("sampleA", "amplicon1")})
        self.assertEqual(len(filter_events), 1)

    def test_existing_feature_types_filter_feature_scope(self):
        feature_registry = {
            "sampleA_amplicon1_ecDNA_1": registry_entry("sampleA", "amplicon1", "ecDNA"),
            "sampleB_amplicon1_BFB_1": registry_entry("sampleB", "amplicon1", "BFB"),
        }
        rows = [similarity_row("sampleA_amplicon1_ecDNA_1", "sampleB_amplicon1_BFB_1", pvalue=1e-10)]

        feats_to_filter, amps_to_filter, _filter_events = ac.choose_similarity_filter_targets(
            rows, feature_registry, pval=1e-5
        )

        self.assertEqual(
            feats_to_filter,
            {("sampleA", "amplicon1", "ecDNA", "1"), ("sampleB", "amplicon1", "BFB", "1")}
        )
        self.assertEqual(amps_to_filter, set())

    def test_linear_requires_asymmetric_scores(self):
        feature_registry = {
            "sampleA_amplicon1_Linear_1": registry_entry("sampleA", "amplicon1", "Linear"),
            "sampleB_amplicon1_Linear_1": registry_entry("sampleB", "amplicon1", "Linear"),
        }
        rows = [
            similarity_row(
                "sampleA_amplicon1_Linear_1", "sampleB_amplicon1_Linear_1",
                pvalue=1e-10, asym1=0.95, asym2=0.5
            )
        ]

        feats_to_filter, amps_to_filter, filter_events = ac.choose_similarity_filter_targets(
            rows, feature_registry, pval=1e-5
        )

        self.assertEqual(feats_to_filter, set())
        self.assertEqual(amps_to_filter, set())
        self.assertEqual(filter_events, [])

    def test_nonpassing_fan_row_does_not_stop_legacy_filtering(self):
        feature_registry = {
            "sampleA_amplicon1_FAN_1": registry_entry("sampleA", "amplicon1", "FAN"),
            "sampleB_amplicon1_ecDNA_1": registry_entry("sampleB", "amplicon1", "ecDNA"),
            "sampleC_amplicon1_ecDNA_1": registry_entry("sampleC", "amplicon1", "ecDNA"),
        }
        rows = [
            similarity_row("sampleA_amplicon1_FAN_1", "sampleB_amplicon1_ecDNA_1", pvalue=0.5),
            similarity_row("sampleB_amplicon1_ecDNA_1", "sampleC_amplicon1_ecDNA_1", pvalue=1e-10),
        ]

        feats_to_filter, amps_to_filter, filter_events = ac.choose_similarity_filter_targets(
            rows, feature_registry, pval=1e-5
        )

        self.assertEqual(
            feats_to_filter,
            {("sampleB", "amplicon1", "ecDNA", "1"), ("sampleC", "amplicon1", "ecDNA", "1")}
        )
        self.assertEqual(amps_to_filter, set())
        self.assertEqual(len(filter_events), 1)

    def test_nonpassing_legacy_non_linear_row_preserves_break_behavior(self):
        feature_registry = {
            "sampleA_amplicon1_ecDNA_1": registry_entry("sampleA", "amplicon1", "ecDNA"),
            "sampleB_amplicon1_BFB_1": registry_entry("sampleB", "amplicon1", "BFB"),
            "sampleC_amplicon1_ecDNA_1": registry_entry("sampleC", "amplicon1", "ecDNA"),
        }
        rows = [
            similarity_row("sampleA_amplicon1_ecDNA_1", "sampleB_amplicon1_BFB_1", pvalue=0.5),
            similarity_row("sampleB_amplicon1_BFB_1", "sampleC_amplicon1_ecDNA_1", pvalue=1e-10),
        ]

        feats_to_filter, amps_to_filter, filter_events = ac.choose_similarity_filter_targets(
            rows, feature_registry, pval=1e-5
        )

        self.assertEqual(feats_to_filter, set())
        self.assertEqual(amps_to_filter, set())
        self.assertEqual(filter_events, [])

    def test_merge_interval_trees_unions_intervals(self):
        target = ac.interval_dict_to_tree({"chr1": [(10, 20)]})
        source = ac.interval_dict_to_tree({"chr1": [(30, 40)], "chr2": [(5, 8)]})

        ac.merge_interval_trees(target, source)

        self.assertTrue(target["chr1"][10:20])
        self.assertTrue(target["chr1"][30:40])
        self.assertTrue(target["chr2"][5:8])


if __name__ == "__main__":
    unittest.main()
