# Changelog

This repository treats `ddca.py` as the authoritative implementation. Locally retained historical and current sources were compared directly. Repeated dataset and experiment copies were grouped by SHA-256 checksum rather than treated as separate versions.

## Private-repository commit anchors

The former private `Machine-learning-research/DDCA-Test/ddca.py` history and diffs supply the following commit-level anchors. Hashes are retained for provenance but are not linked because the legacy repository is not publicly accessible. The local numbered files additionally preserve experimental snapshots that were not committed on the main `ddca.py` path.

| Commit date | Commit | Recorded subject | Evidential use |
|---|---|---|---|
| 2021-06-27 | `12bcb7f` | `Create bd-ddca_v0.1.py` | Creates the first tracked Python scaffold. |
| 2021-06-29 | `ef52b7a` | `initial` | Adds Plaso dynamic-CSV discovery, time resampling, initial safe/danger derivation and graphing. |
| 2021-06-30 | `4fbf76a` | `Update bd-ddca_v0.1.py` | Adds the first working dense-array dDCA loop, migration/logging, MCAV/Kalpha calculation and output graphs. |
| 2022-02-12 | `0fa47fa` | `Updates to initial implementation` | Updates v0.1 and adds the tracked v0.2 and v0.3 files. |
| 2022-09-25 | `717c01e` | `Create bd-ddca_v1.0.py` | Adds the first object-based v1.0 file with `dDCA_cells` and `dDCA` classes. |
| 2022-09-28 | `2a3c689` | `2022/09/25 Implement object based code` | Renames v1.0 to `DDCA-Test/ddca.py`, moves reset behaviour into the cell class, tidies aggregate logging/statistics and removes the tracked v0.1–v0.3 files. |
| 2022-10-02 | `5ca03d8` | `Initial optimisation of logAntigen` | Replaces the inner per-occurrence loop with aggregate arithmetic and resets the dense per-cell profile in one operation. |
| 2024-07-01 | `3973695` | `Update ddca.py` | Implements the 2024-06-30 tidy: array arithmetic in `logAntigen()`, direct cell iteration in `doSignals()`, removal of the separate `updateDC()`, reproducible Docker notes and a CC licence marker. |
| 2024-07-06 | `f124f89` | `Update ddca.py` | Precomputes `csm`/`k` and changes `logAntigen()` to accept a cell object; the commit temporarily leaves the final flush passing integer indices. |
| 2024-07-09 | `b73a480` | `Update ddca.py` | Introduces dynamic per-cell antigen lists, fixed typed aggregate arrays, cycle-based assignment and the correct cell-object flush, while removing legacy per-cell metrics. |
| 2024-07-09 | `ecc3983` | `Update ddca.py` | Corrects two comments only: “historicl” and “Adcance”. |

## Current public release — 2026

- Sort each cell's dynamic antigen list before aggregate logging, matching the ascending antigen traversal of the shared C implementation.
- Correct the initial sorting expression: `list.sort()` mutates a list and returns `None`, so sorting and iteration must be separate operations.
- Clarify that cycle pre-advancement reproduces the C implementation's pre-incremented cell index; for populations larger than one, the first antigen is assigned to cell 1.
- Establish dDCA as a separately packaged, tested, Cython-buildable research component.
- Publish the authoritative module as `ddca.py`, with release information kept in package metadata rather than encoded in the filename.
- Retain lists for dynamic per-cell antigen collections and NumPy arrays for fixed-size aggregate profiles.
- Deliberately preserve floating-point migration thresholds across resets. The shared C implementation changes from a floating initial interval to an integer reset interval; this truncation is believed unintentional.

## v8 — 2025-12

- Handle a cell population of one without dividing by `TotCells - 1` while distributing migration lifespans.

## v7 — 2024-11 to 2024-12

- 2024-11-28: add `Tmcav`, track total danger and safe signals separately, and calculate the MCAV threshold as their ratio.
- Retain the v6 Kalpha threshold while expressing its weighted sum through `TotDanger - 2*TotSafe`.
- Remove the embedded experiment runner from the core module.
- 2024-12-04: restore `dc_stats()` with a lazy pandas import and per-cell identifiers.
- 2024-12-10: parenthesise the Kalpha threshold expression for clarity without changing precedence.

## v6 — 2024-11-25

- Restore the Kalpha anomaly threshold (`Tk`) and the signal-pair count used to calculate it.
- Restore per-cell iteration and incarnation statistics required by the threshold calculation.
- Add a pandas-backed `dc_stats()` helper in the development snapshot.

## v5 — 2024-07-09

- Replace each dense per-cell antigen profile with a dynamic Python list while retaining NumPy arrays for fixed-size aggregate profiles.
- Restore true multi-antigen collection after the v4 single-antigen experiment.
- Use an `itertools.cycle` object for deterministic round-robin assignment.
- Record a runtime reduction from 9.66 seconds to 0.64 seconds for 13,128 entries on an Apple M3. The original benchmark fixture and memory measurements have not yet been recovered.
- Commit `b73a480` introduces the list/array hybrid and corrects the final flush to pass cell objects; `ecc3983` is comment-only cleanup.

## v4 — 2024-07-08

- Trial a single-antigen-per-cell representation using one antigen index and counter.
- Record that the preceding threading experiment was ineffective for the observed workload.
- Use fixed integer/float aggregate arrays and remove several legacy metrics variables.
- This numbered source is a retained experimental snapshot. It does not appear as a committed state of `DDCA-Test/ddca.py`; commits `f124f89` and `b73a480` bracket it chronologically.

## v1.0 and early object implementation — 2022-09 to 2024-06

- 2022-09-25 (`717c01e`): create `bd-ddca_v1.0.py`, extracting object-based `dDCA_cells` and `dDCA` classes from the monolithic analysis script.
- 2022-09-28 (`2a3c689`): rename the file to `DDCA-Test/ddca.py`, add `dDCA_cells.reset()` and consolidate reset/statistics handling.
- 2022-10-02 (`5ca03d8`): perform the initial `logAntigen()` aggregate-arithmetic optimisation.
- Store a dense NumPy antigen-count array in every cell and vectorise aggregate logging.
- Preserve the explicit `doAntigen()` then `doSignals()` call sequence established in v0.7.
- Add a single-cell guard in one retained monolithic snapshot, although the numbered implementation did not carry it forward until v8.
- Replace suitable loops with vector operations and separate more of the algorithm from the surrounding analysis script.

## v0.10 — 2022-09

- Change safe-signal calculation from an inline `np.diff` pipeline expression to separately developed logic.
- Continue clustering experiments using KMeans and one-class SVM in the surrounding analysis script; these were not carried into the extracted dDCA engine.

## v0.9 clustering branch — 2022-08

- Add KMeans clustering experiments around the timeline-derived dDCA output.
- Treat this as an experimental pipeline branch, not the direct core-engine predecessor of v1.0.

## v0.8 — 2022-04 to 2022-07

- Restrict analysis to configured date ranges.
- Change danger and safe signal scaling to logarithmic ratios.
- Retain antigen-before-signal ordering.
- The v0.8 header says `2020/04/10` and `2020/04/11`; comparison with v0.7, the file's 2022 modification date, and all later headers indicate these are typographical errors for 2022.

## v0.7 — 2022-04-02

- Move antigen assignment before signal processing in the main loop. This ordering survives in the current implementation and is supported by Gu, Greensmith and Aickelin's segmented workflow.

## v0.5–v0.6 — 2022-03 to 2022-04

- Refactor the monolithic runner into a callable `main()` path and improve repeated-file processing.
- Simplify signal-loop data retained for the analysis and preserve cached derived `csm` and `k` values.
- Carry an explicit open question about antigen/signal order, resolved in v0.7.

## v0.4 — 2022-03-10

- Derive danger input from the number of YARA matches exported with the Plaso timeline.
- Account for invalid or unset timestamps during the surrounding timeline preprocessing.

## v0.3 and earlier pipeline — 2021-06 to 2022-02

- 2021-06-27 (`12bcb7f`): create the first tracked v0.1 Python scaffold.
- 2021-06-29 (`ef52b7a`): add Plaso dynamic-CSV discovery, pandas time resampling, initial activity-derived danger/safe values and graphs.
- 2021-06-30 (`4fbf76a`): implement the first full Python dDCA loop with dense per-cell antigen arrays, migration, mature/semi-mature aggregation, MCAV/Kalpha calculation and graphs.
- 2021-07-29 (source-header record): improve graphing in the monolithic pipeline script.
- 2021-07-30 (source-header record): add pickle-based caching around the pipeline.
- 2022-02-12 (`0fa47fa`): update v0.1 and add tracked v0.2/v0.3 implementations. v0.3 combines Plaso loading, time resampling, safe/danger derivation, dense cell/antigen arrays, dDCA execution, output persistence and graphing in one script.
- 2022-09-28 (`2a3c689`): remove tracked v0.1–v0.3 files when the object-based v1.0 implementation is promoted to `ddca.py`; separately retained copies preserve the experimental history.

The v0.x entries describe the former surrounding analysis pipeline as well as algorithm changes. The original repository now provides commit anchors from the first v0.1 scaffold through the July 2024 list-based implementation. Retained numbered sources fill experimental stages developed outside, between, or after those committed `ddca.py` states.
