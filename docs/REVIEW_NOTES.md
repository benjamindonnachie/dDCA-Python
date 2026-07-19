# Implementation review notes

The public source comments have been revised where the retained C implementation and local history provide direct evidence. The following behaviours still need a research decision before they should be changed; the current implementation preserves them for compatibility.

## Comments resolved for this release

- The old “start with first cell” comment was misleading. The cycle is pre-advanced so the first antigen goes to cell 1 when multiple cells exist, matching the C implementation's increment-before-assignment behaviour.
- The placeholder punctuation around the Gu, Greensmith and Aickelin citation was replaced with a verified bibliographic reference.
- The antigen-order comment now states exactly why the per-cell list is sorted.
- The licence identifier was corrected from the non-standard `CC-BY-SA-NC` order to `CC-BY-NC-SA-4.0`.
- Pipeline-era comments about Plaso, graphing, pickles and YARA were removed from the module because those operations are not performed by dDCA.
- The v0.8 source's `2020/04/10` and `2020/04/11` history dates were treated as typographical errors for 2022: v0.7 is dated 2022-04-02, v0.8 was modified in July 2022, and later headers use 2022.
- The shared C implementation initialises the migration interval as `float` but declares the reset interval as `int`. This is believed unintentional. dDCA deliberately preserves the floating-point maximum lifespan across resets; strict C-executable reproduction would require integer truncation.

## Behavioural assumptions to revisit

1. **Input validation.** `cells`, `antigenCount`, `maxMigration`, antigen IDs and signals are not validated. A negative antigen ID uses NumPy negative indexing; an ID above the configured range raises `IndexError`. Decide whether the surrounding pipeline remains responsible or a future release should reject invalid input explicitly.
2. **Zero denominators.** Calling `results()` without signal pairs divides by zero while calculating `Tk`; a zero total safe signal divides by zero while calculating `Tmcav`. Unobserved antigen types yield `NaN` for MCAV and Kalpha. Decide whether to raise, return `NaN`, or use a configured policy.
3. **Cell 0 migration threshold.** For multiple cells, the inherited expression gives cell 0 a zero maximum lifespan. This matches the shared C implementation, but its analytical purpose should be confirmed before refactoring.
4. **Incomplete active-incarnation statistics.** `TotIter` is updated only on cell reset, so `results()` does not include a cell's unfinished current incarnation when calculating the average used by `Tk`.
5. **Signal domain.** The implementation accepts negative and non-finite safe/danger values. The intended signal domain and normalisation contract belong in the wider pipeline specification.
6. **Output precision.** Aggregate Kalpha storage uses NumPy's `"f"` dtype (normally `float32`) while Python signal totals use native floating point. Confirm whether reproducibility requires explicit precision.
7. **Legacy API names.** `dDCA`, `dDCA_cells`, `TotDanger`, `Is`, and similar names are retained for compatibility. A future version could expose PEP 8 names through a compatibility layer.
8. **Unused compatibility parameter.** `dDCA_cells(..., antigenCount)` retains an unused argument from earlier versions. Remove only with a deliberate API change.
9. **Historical performance claim.** The 9.66-to-0.64-second observation survives in the source history, but its exact input, command, dependency versions, and memory measurements are not preserved. Treat it as historical context until benchmark fixtures are recovered.

These are review flags, not claims that the current behaviour is incorrect. Algorithm changes should be accompanied by a new version, fixed fixtures, and comparison against both pure-Python and Cython builds.
