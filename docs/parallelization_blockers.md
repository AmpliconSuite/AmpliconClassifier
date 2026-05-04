# Parallelization Blockers in AmpliconClassifier

AmpliconClassifier processes amplicons one at a time in a sequential loop. This document catalogs the concrete blockers to parallelizing at the amplicon level (e.g., using `multiprocessing` or `concurrent.futures`).

---

## Blocker 1 — Global mutable result accumulators
**Severity: Significant redesign**

The loop appends results into ~12 module-level lists and dicts declared at lines 1412–1427 of `amplicon_classifier.py`:

```python
ftgd_list, ftci_list, bpgi_list, fd_list, prop_list,
AMP_dvaluesList, AMP_classifications, bfbarchitect_summaries,
sampNames, cyclesFiles, samp_to_ec_count,
featEntropyD, full_featname_to_graph, full_featname_to_intervals
```

These are written in `run_classification()` (lines ~1131, 1143, 1148, 1153–1177) and `run_bfbarchitect()` (line 1095). They are consumed by `write_outputs()` and `filter_similar_amplicons()` after the loop ends.

In `multiprocessing`, forked workers get copies so writes are silently lost in the main process. With threads, concurrent appends corrupt the lists.

**Fix required:** `run_classification()` and `run_bfbarchitect()` must return their results as values rather than side-effecting globals. The main process collects them (e.g., via `Pool.map`) and populates the accumulators after all workers finish. This touches the core loop structure and both functions substantially.

---

## Blocker 2 — `args` namespace accessed deep in the call chain
**Severity: Moderate**

The argparse `args` object is a module-level global read inside many per-amplicon functions called from the loop:

| Function | Access | Location |
|---|---|---|
| `run_classification()` | `args.add_chr_tag`, `args.ref` | Lines ~1060, 1151 |
| `run_bfbarchitect()` | `args.bfbarchitect`, `args.bfb_threads` | Lines 730, 744, 748 |
| `get_raw_cycle_props()` | `args.ref`, `args.force` | Lines ~1020, 1027 |
| `nonbfb_cycles_are_ecdna()`, `clusterECCycles()`, `cycleIsNoAmpInvalid()` | `args.ref` | Lines ~264, 471, 529 |

`args` is read-only after line 1244, so it is safe to share — but it must be explicitly passed into worker processes. With `multiprocessing`, the entire namespace pickles fine (it contains only primitives and strings). The work is mechanical but widespread: every worker-callable function that reads `args` directly needs it threaded in as a parameter, or the worker initializer must re-establish it as a process-local global.

---

## Blocker 3 — `BFBARCHITECT_RECONSTRUCT` / `BFBARCHITECT_CENTROMERES` globals
**Severity: Moderate**

`BFBARCHITECT_RECONSTRUCT` is assigned at line 1385 to the imported function `reconstruct_bfb_from_graph`. `BFBARCHITECT_CENTROMERES` is a dict loaded at line 1386. Both are read inside `run_bfbarchitect()` (lines 736–748), which is called per-amplicon.

`multiprocessing` with the `fork` start method (Linux default) handles this transparently — both globals exist in forked workers. With `spawn` (required on macOS/Windows, also safer), the function reference pickles correctly only if BFBArchitect's module is importable in workers, and the centromere dict pickles as a plain dict.

A secondary concern: `--bfb_threads` controls ILP solver threads per call. Running N amplicons in parallel, each using M solver threads, gives N×M threads. These need to be budgeted together.

---

## Blocker 4 — Logger and shared log file
**Severity: Easy**

`logger` is a module-level global (line 1270) backed by a `FileHandler` writing to `{prefix}.log`. Multiple workers writing concurrently will interleave or corrupt log lines.

**Fix:** Use `logging.handlers.QueueHandler` in workers and a `QueueListener` in the main process. This is a standard pattern and self-contained.

---

## Blocker 5 — `ConfigVars` class-level mutable attributes
**Severity: Easy**

`ConfigVars` (`ac_util.py` lines 12–33) uses class attributes as config storage. `set_config_vars()` (`ac_util.py` ~line 549) mutates them at startup via `setattr(ConfigVars, key, value)`. This happens before the loop — during the loop, `ConfigVars` is read-only everywhere.

With `fork`, workers inherit the already-configured state. With `spawn`, `ConfigVars` needs to be initialized in each worker (straightforward: pass the config dict and call `set_config_vars` in the worker initializer). Not a blocker in practice.

---

## Blocker 6 — `filter_similar_amplicons()` mutates globals post-loop
**Severity: Moderate (conditional)**

When `--filter_similar` is passed, `filter_similar_amplicons()` (called at line 1458) directly mutates the global accumulators from Blocker 1 — deleting entries from `ftgd_list`, `fd_list`, `prop_list`, `AMP_classifications`, decrementing `samp_to_ec_count`, etc. It is a post-loop operation, so it does not block parallelizing the loop itself, but fixing Blocker 1 (returning results rather than appending to globals) requires this function to be adapted to operate on the collected results instead.

---

## Blocker 7 — Per-amplicon file writes inside `write_outputs()`
**Severity: Easy**

`write_outputs()` in `ac_io.py` (lines ~428–489) is called once after the loop and iterates over the accumulated lists to write SV summaries, BED files, and annotated cycles. This is sequential post-processing and is not itself a blocker. Individual per-amplicon subdirectory writes (`write_bpg_summary`, `write_interval_beds`) open separate files per amplicon. If those were moved inside the loop (to avoid holding all results in memory), they would need path-collision guards but are otherwise parallelism-safe.

---

## Summary

| # | Issue | Location | Severity |
|---|---|---|---|
| 1 | Global mutable result accumulators written during loop | `amplicon_classifier.py` lines 1412–1427, writes throughout loop | **Significant** |
| 2 | `args` global read inside many per-amplicon functions | Pervasive — `run_classification`, `run_bfbarchitect`, helpers | **Moderate** |
| 3 | `BFBARCHITECT_*` globals + nested solver threading | `amplicon_classifier.py` lines 26–27, 1383–1387, 729–755 | **Moderate** |
| 4 | Shared logger / log file handle | `amplicon_classifier.py` line 1270, `ac_io.py` lines 10–23 | **Easy** |
| 5 | `ConfigVars` class attribute mutation at startup | `ac_util.py` lines 12–33, 549–562 | **Easy** |
| 6 | `filter_similar_amplicons()` mutates globals post-loop | `amplicon_classifier.py` lines 1458–1461 | **Moderate** (conditional) |
| 7 | Post-loop file I/O in `write_outputs()` | `ac_io.py` lines 428–489 | **Easy** |

**The single hardest blocker is #1.** Everything else is routine. The cleanest path to parallelization is: define a `process_amplicon(fpair, args, ...)` function that returns a result dict, run it with `concurrent.futures.ProcessPoolExecutor`, and collect results in the main process before passing them to `write_outputs()`.
