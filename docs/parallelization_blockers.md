# Parallelization Status in AmpliconClassifier

This note originally documented blockers to amplicon-level parallelization. The major blockers have now been addressed by moving per-amplicon state into record objects and collecting worker results in input order.

## Current State

- Amplicons can be classified in parallel with `--jobs N`.
- Output row order follows input manifest order.
- Per-amplicon classification returns structured records instead of appending directly to shared global accumulators.
- Output writing remains a sequential post-classification stage.
- BFBArchitect integration is default-on when installed and can be disabled with `--no_bfbarchitect`.
- BFBArchitect thread use is nested under AmpliconClassifier workers, so the approximate solver thread budget is `--jobs * --bfb_threads`.

## Resolved Blockers

| Issue | Current handling |
|---|---|
| Shared result accumulators | Worker results are collected as records and merged in the parent process. |
| Per-amplicon dictionary sprawl | Core per-amplicon inputs/results/features are represented with lightweight objects in `ampclasslib/classification_records.py`. |
| Input-order stability | Parallel results are collected against the original input list, preserving deterministic output order. |
| BFBArchitect subprocess/thread interaction | Each classification worker can run BFBArchitect, with `--bfb_threads` controlling per-amplicon solver threads. |
| Per-amplicon output collisions | Feature/annotation files are still emitted after classification from collected records. |

## Remaining Caveats

- `--jobs * --bfb_threads` can oversubscribe CPUs on large runs; choose both values together.
- Logging from workers is intentionally summarized rather than fully streaming every low-level BFBArchitect message.
- Spawn-style multiprocessing portability should be kept in mind if supporting non-Linux platforms; Linux fork behavior is the main exercised path.
- `--filter_similar` remains a post-classification operation and should be smoke-tested after changes to record collection or feature identity.

## Useful Validation

Run the focused unit tests:

```bash
python -m unittest discover -s tests
```

Run a representative parallel smoke test:

```bash
python amplicon_classifier.py --ref GRCh38 --input collection.input --jobs 3 --bfb_threads 2 --verbose_classification -o debug_prefix
```

When BFBArchitect is available, compare the same input with `--no_bfbarchitect` to confirm native AC BFB calls still behave as expected.
