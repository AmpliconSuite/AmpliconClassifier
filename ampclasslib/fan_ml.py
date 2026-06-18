"""
fan_ml — ML classifier for FAN (Focal amplification in neochromosome) detection.

Model: model_3_5_rank05_20260514_142539 (L2 logistic regression, 5 features)
Balanced accuracy: 0.961 (TCGA CV, grouped by patient), PR-AUC: 0.988

Public API
----------
    # AC integration — pass pre-parsed edges to avoid re-reading the graph file:
    result = classify_from_edges(seq_edges, sv_edges)

    # standalone / CLI — reads the file internally:
    result = classify_graph(graph_path)

    # result: {"decision": "FAN" | "not_FAN",
    #          "probability": float, "features": dict}

seq_edges format: list of dicts with keys chrom, start, end, cn, size_bp
sv_edges  format: list of dicts with keys chrom1, pos1, strand1,
                                          chrom2, pos2, strand2, cn, sv_type
"""

from __future__ import annotations

import math
import re
import sys
from collections import defaultdict

# ── Hardcoded model weights (model_3_5_rank05_20260514_142539) ────────────────
# Features in model order
_FEATURE_NAMES = [
    "amplified_span_mb",
    "sv_qualifying_inter_or_large",
    "chrom_dominant_frac",
    "sv_crossing_frac",
    "amplified_cn40_span_mb",
]
# Features log1p-transformed before scaling
_LOG1P_FEATS = {"amplified_span_mb", "sv_qualifying_inter_or_large", "amplified_cn40_span_mb"}
# Features zero-imputed when NaN (before log1p)
_ZERO_IMPUTE_FEATS = {"sv_crossing_frac", "amplified_cn40_span_mb"}
# Training-set medians for any remaining NaN imputation (after log1p)
_TRAINING_MEDIANS = {
    "amplified_span_mb":            3.9658,
    "sv_qualifying_inter_or_large": 2.8537,
    "chrom_dominant_frac":          0.9186,
    "sv_crossing_frac":             0.5833,
    "amplified_cn40_span_mb":       0.0,
}
# StandardScaler: fit on training data
_MEAN  = [3.965828078652812, 2.853699887324316, 0.8112811240721095,
          0.5730005302226935, 0.37467435011187045]
_SCALE = [0.8978377894837156, 1.220536575460846, 0.21568088034489855,
          0.2521635023323574, 1.017638668508902]
# Logistic regression coefficients (one per feature) and intercept
_COEF      = [0.32133598192120943, 5.3712848707745815, -1.3620222815703678,
              0.8605780334474445, -1.3488868597793138]
_INTERCEPT = -1.2625938392862028
_THRESHOLD = 0.5

# ── SV classification constants ───────────────────────────────────────────────
_PREREQ_SV_DEL_DUP_MIN_BP = 50_000
_PREREQ_SV_LARGE_BP       = 1_000_000


# ── Graph parsing (used by CLI / standalone classify_graph only) ──────────────

def _parse_vertex(tok: str) -> tuple:
    m = re.match(r'^(.+):(\d+)([+-])$', tok)
    if not m:
        raise ValueError("Cannot parse vertex token: {!r}".format(tok))
    return m.group(1), int(m.group(2)), m.group(3)


def _sv_type(chrom1, pos1, strand1, chrom2, pos2, strand2) -> str:
    if chrom1 != chrom2:
        return "interchromosomal"
    if strand1 == strand2:
        return "inversion"
    if strand1 == "+" and strand2 == "-":
        return "deletion" if pos2 > pos1 else "duplication"
    return "duplication" if pos2 > pos1 else "deletion"


def _parse_graph_and_sv_edges(graph_path: str) -> tuple:
    """Read a graph file once, returning (seq_edges, sv_edges)."""
    seq_edges = []
    sv_edges  = []
    with open(graph_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            if parts[0] == "sequence":
                try:
                    chrom, pos_start, _ = _parse_vertex(parts[1])
                    _,     pos_end,   _ = _parse_vertex(parts[2])
                    cn      = float(parts[3])
                    size_bp = int(parts[5])
                except (ValueError, IndexError):
                    continue
                seq_edges.append({"chrom": chrom, "start": pos_start,
                                   "end": pos_end, "cn": cn, "size_bp": size_bp})
            elif parts[0] == "discordant":
                try:
                    v1_tok, v2_tok = parts[1].split("->", 1)
                    chrom1, pos1, strand1 = _parse_vertex(v1_tok)
                    chrom2, pos2, strand2 = _parse_vertex(v2_tok)
                    cn  = float(parts[2])
                    svt = _sv_type(chrom1, pos1, strand1, chrom2, pos2, strand2)
                except (ValueError, IndexError):
                    continue
                sv_edges.append({"chrom1": chrom1, "pos1": pos1, "strand1": strand1,
                                  "chrom2": chrom2, "pos2": pos2, "strand2": strand2,
                                  "cn": cn, "sv_type": svt})
    return seq_edges, sv_edges


# ── Feature extraction ────────────────────────────────────────────────────────

def _is_short_intrachrom_del_dup(edge: dict) -> bool:
    return (
        edge.get("sv_type") in {"deletion", "duplication"}
        and edge["chrom1"] == edge["chrom2"]
        and abs(edge["pos2"] - edge["pos1"]) < _PREREQ_SV_DEL_DUP_MIN_BP
    )


def _is_inter_or_large(edge: dict) -> bool:
    return (
        edge["chrom1"] != edge["chrom2"]
        or abs(edge["pos2"] - edge["pos1"]) > _PREREQ_SV_LARGE_BP
    )


def _envelope_span_mb(edges: list) -> float:
    if not edges:
        return 0.0
    by_chrom: dict = defaultdict(list)
    for e in edges:
        by_chrom[e["chrom"]].append(e)
    span_bp = sum(
        max(e["end"] for e in ce) - min(e["start"] for e in ce)
        for ce in by_chrom.values()
    )
    return span_bp / 1e6


def extract_features_from_edges(seq_edges: list, sv_edges: list,
                                 cn_floor: float = 4.0) -> dict:
    """Compute the 5 model features from pre-parsed graph edges."""
    amp_edges    = [e for e in seq_edges if e["cn"] > cn_floor]
    amp_cn40     = [e for e in seq_edges if e["cn"] > 40.0]
    qualifying   = [e for e in sv_edges  if not _is_short_intrachrom_del_dup(e)]

    amplified_span_mb = round(_envelope_span_mb(amp_edges), 4)

    sv_inter_large = sum(1 for e in qualifying if _is_inter_or_large(e))

    if seq_edges:
        by_chrom: dict = defaultdict(int)
        for e in seq_edges:
            by_chrom[e["chrom"]] += e["size_bp"]
        total_bp = sum(by_chrom.values())
        chrom_dom = round(max(by_chrom.values()) / total_bp, 4) if total_bp else 0.0
    else:
        chrom_dom = 0.0

    n = len(qualifying)
    if n > 0:
        crossing = 0
        for i, sv in enumerate(qualifying):
            if sv["chrom1"] != sv["chrom2"]:
                continue
            lo    = min(sv["pos1"], sv["pos2"])
            hi    = max(sv["pos1"], sv["pos2"])
            chrom = sv["chrom1"]
            for j, other in enumerate(qualifying):
                if j == i:
                    continue
                if ((other["chrom1"] == chrom and lo < other["pos1"] < hi) or
                        (other["chrom2"] == chrom and lo < other["pos2"] < hi)):
                    crossing += 1
                    break
        sv_crossing = round(crossing / n, 4)
    else:
        sv_crossing = 0.0

    cn40_span_mb = round(_envelope_span_mb(amp_cn40), 4)

    return {
        "amplified_span_mb":            amplified_span_mb,
        "sv_qualifying_inter_or_large": sv_inter_large,
        "chrom_dominant_frac":          chrom_dom,
        "sv_crossing_frac":             sv_crossing,
        "amplified_cn40_span_mb":       cn40_span_mb,
    }


# ── Scoring ───────────────────────────────────────────────────────────────────

def _score(feats: dict) -> float:
    """Apply preprocessing pipeline and return P(FAN)."""
    row = []
    for i, name in enumerate(_FEATURE_NAMES):
        val = feats.get(name, float("nan"))
        if name in _ZERO_IMPUTE_FEATS and isinstance(val, float) and math.isnan(val):
            val = 0.0
        if name in _LOG1P_FEATS:
            val = math.log1p(float(val))
        if isinstance(val, float) and math.isnan(val):
            val = _TRAINING_MEDIANS.get(name, 0.0)
        row.append((float(val) - _MEAN[i]) / _SCALE[i])

    logit = sum(c * x for c, x in zip(_COEF, row)) + _INTERCEPT
    return 1.0 / (1.0 + math.exp(-logit))


# ── Public API ────────────────────────────────────────────────────────────────

def classify_from_edges(seq_edges: list, sv_edges: list,
                        cn_floor: float = 4.0) -> dict:
    """
    Classify an amplicon given pre-parsed graph edges.
    Preferred entry point for AmpliconClassifier integration —
    avoids re-reading the graph file.
    """
    feats = extract_features_from_edges(seq_edges, sv_edges, cn_floor)
    prob  = _score(feats)
    label = "FAN" if prob >= _THRESHOLD else "not_FAN"
    return {"decision": label, "probability": round(prob, 4), "features": feats}


def classify_graph(graph_path: str, cn_floor: float = 4.0) -> dict:
    """
    Classify one AmpliconArchitect graph file.
    Reads the file once internally. Prefer classify_from_edges when data
    is already parsed (e.g. inside AmpliconClassifier).
    """
    seq_edges, sv_edges = _parse_graph_and_sv_edges(str(graph_path))
    return classify_from_edges(seq_edges, sv_edges, cn_floor)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Classify an AA graph file as FAN (Focal amplification in neochromosome) or not_FAN."
    )
    parser.add_argument("graph", metavar="GRAPH_TXT",
                        help="Path to *_amplicon*_graph.txt")
    parser.add_argument("--cn-floor", type=float, default=4.0,
                        help="CN floor for amplified segments (default 4.0)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print feature values alongside the result")
    args = parser.parse_args()

    result = classify_graph(args.graph, cn_floor=args.cn_floor)
    print("{}\t{:.4f}".format(result["decision"], result["probability"]))
    if args.verbose:
        for k, v in result["features"].items():
            print("  {}: {}".format(k, v))


if __name__ == "__main__":
    _main()
