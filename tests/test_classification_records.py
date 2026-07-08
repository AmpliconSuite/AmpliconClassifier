import unittest

import _pathfix  # noqa: F401
from ampclasslib.classification_records import AmpliconRecord, collect_amplicon_records


def make_record(sample_name, amplicon_number, ec_count=0):
    feature_id = "{}_{}_ecDNA_1".format(sample_name, amplicon_number)
    return AmpliconRecord(
        sample_name=sample_name,
        amplicon_number=amplicon_number,
        cycles_file="{}_{}_cycles.txt".format(sample_name, amplicon_number),
        graph_file="{}_{}_graph.txt".format(sample_name, amplicon_number),
        feature_gene_entry=[sample_name, amplicon_number, {}, {}],
        annotated_cycles_entry=["{}_{}_annotated_cycles.txt".format(sample_name, amplicon_number), [], [], {}, [], [], [], set()],
        breakpoint_graph_lines=[[sample_name, amplicon_number]],
        feature_dict={"ecDNA_1": {"chr1": [(1, 10)]}},
        prop_dict={"ecDNA_1": [9, 5.0, 6.0]},
        classification=("Cyclic", ec_count > 0, False, ec_count),
        dvalues=[float(ec_count)],
        bfbarchitect_summary={"min_score": "NA"},
        fan_result={"decision": "not_FAN", "probability": 0.0, "features": {}},
        feature_complexity={(sample_name, amplicon_number, "ecDNA_1"): (1.0, 2.0, 3.0)},
        ec_count=ec_count,
        samp_amp_to_graph={"{}_{}".format(sample_name, amplicon_number): "{}_{}_graph.txt".format(sample_name, amplicon_number)},
        feature_registry={feature_id: {"sample_name": sample_name, "amplicon_number": amplicon_number}},
        full_featname_to_graph={feature_id: "{}_{}_graph.txt".format(sample_name, amplicon_number)},
        full_featname_to_intervals={feature_id: {"chr1": [(1, 10)]}},
    )


class ClassificationRecordsTests(unittest.TestCase):
    def test_collect_amplicon_records_preserves_order_and_merges_fields(self):
        records = [
            make_record("sampleA", "amplicon1", ec_count=1),
            make_record("sampleA", "amplicon2", ec_count=2),
        ]

        results = collect_amplicon_records(records)

        self.assertEqual(results.records, records)
        self.assertEqual(results.sampNames, ["sampleA", "sampleA"])
        self.assertEqual(results.cyclesFiles, ["sampleA_amplicon1_cycles.txt", "sampleA_amplicon2_cycles.txt"])
        self.assertEqual(results.AMP_classifications, [records[0].classification, records[1].classification])
        self.assertEqual(results.samp_to_ec_count["sampleA"], 3)
        self.assertIn(("sampleA", "amplicon1", "ecDNA_1"), results.featComplexityD)
        self.assertIn("sampleA_amplicon2_ecDNA_1", results.feature_registry)
        self.assertEqual(
            results.samp_amp_to_graph["sampleA_amplicon1"],
            "sampleA_amplicon1_graph.txt"
        )


if __name__ == "__main__":
    unittest.main()
