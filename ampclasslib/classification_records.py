from dataclasses import dataclass, field
from collections import defaultdict


@dataclass
class AmpliconInput:
    original_sample_name: str
    sample_name: str
    amplicon_number: str
    cycles_file: str
    graph_file: str


@dataclass
class ClassificationContext:
    args: object
    add_chr_tag: bool
    decomposition_strictness: float
    lcD: dict
    patch_links: list
    gene_lookup: dict
    ncRNA_lookup: dict


@dataclass
class AmpliconRecord:
    sample_name: str
    amplicon_number: str
    cycles_file: str
    graph_file: str
    feature_gene_entry: list
    annotated_cycles_entry: list
    breakpoint_graph_lines: list
    feature_dict: dict
    prop_dict: dict
    classification: tuple
    dvalues: list
    bfbarchitect_summary: dict
    fan_result: dict
    feature_complexity: dict = field(default_factory=dict)
    ec_count: int = 0
    samp_amp_to_graph: dict = field(default_factory=dict)
    feature_registry: dict = field(default_factory=dict)
    full_featname_to_graph: dict = field(default_factory=dict)
    full_featname_to_intervals: dict = field(default_factory=dict)


@dataclass
class Feature:
    feature_id: str
    feature_type: str
    intervals: dict
    cycle_indices: set = field(default_factory=set)
    sources: set = field(default_factory=set)
    score: object = None
    metadata: dict = field(default_factory=dict)


@dataclass
class ClassificationResults:
    records: list = field(default_factory=list)
    ftgd_list: list = field(default_factory=list)
    ftci_list: list = field(default_factory=list)
    bpgi_list: list = field(default_factory=list)
    fd_list: list = field(default_factory=list)
    prop_list: list = field(default_factory=list)
    AMP_dvaluesList: list = field(default_factory=list)
    AMP_classifications: list = field(default_factory=list)
    bfbarchitect_summaries: list = field(default_factory=list)
    fan_results: list = field(default_factory=list)
    sampNames: list = field(default_factory=list)
    cyclesFiles: list = field(default_factory=list)
    featComplexityD: dict = field(default_factory=dict)
    samp_to_ec_count: dict = field(default_factory=lambda: defaultdict(int))
    samp_amp_to_graph: dict = field(default_factory=dict)
    full_featname_to_graph: dict = field(default_factory=dict)
    full_featname_to_intervals: dict = field(default_factory=dict)
    feature_registry: dict = field(default_factory=dict)


def collect_amplicon_records(records):
    results = ClassificationResults(records=list(records))
    for record in records:
        results.ftgd_list.append(record.feature_gene_entry)
        results.ftci_list.append(record.annotated_cycles_entry)
        results.bpgi_list.append(record.breakpoint_graph_lines)
        results.fd_list.append(record.feature_dict)
        results.prop_list.append(record.prop_dict)
        results.AMP_dvaluesList.append(record.dvalues)
        results.AMP_classifications.append(record.classification)
        results.bfbarchitect_summaries.append(record.bfbarchitect_summary)
        results.fan_results.append(record.fan_result)
        results.sampNames.append(record.sample_name)
        results.cyclesFiles.append(record.cycles_file)
        results.featComplexityD.update(record.feature_complexity)
        results.samp_to_ec_count[record.sample_name] += record.ec_count
        results.samp_amp_to_graph.update(record.samp_amp_to_graph)
        results.full_featname_to_graph.update(record.full_featname_to_graph)
        results.full_featname_to_intervals.update(record.full_featname_to_intervals)
        results.feature_registry.update(record.feature_registry)

    return results
