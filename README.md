# dDCA

dDCA is a compact Python implementation of the deterministic Dendritic Cell Algorithm for correlating integer antigen identifiers with danger and safe signals in a research data pipeline.

> **Project status: work in progress.** This repository is the authoritative home of the Python implementation, but the code is a research pipeline component rather than a standalone product, detector, or verdict engine. Its interfaces and assumptions may change as the wider research develops.

The module retains dynamic per-cell antigen collections as Python lists while using NumPy arrays only for fixed-size aggregate output profiles. For the intended workload, lists avoid repeatedly growing dense arrays in the hot path and reduce RAM and CPU use. A recorded 2024 test of the wider implementation processed 13,128 entries in 0.64 seconds on an Apple M3, down from 9.66 seconds; that historical observation is not a general benchmark and should be reproduced for other environments.

## Place in the pipeline

dDCA does not ingest forensic images, Plaso output, network traffic, or raw logs. Upstream pipeline stages are responsible for acquiring and normalising evidence, assigning a bounded integer antigen ID to each entity of interest, and deriving numeric danger and safe signals. The dDCA stage then correlates those inputs across a deterministic population of cells. Downstream stages join the resulting metrics back to the source data for analysis, visualisation, ground-truth comparison, and human evidential review.

```mermaid
flowchart LR
    A["Evidence and source data"] --> B["Normalise and enrich"]
    B --> C["Derive antigen IDs and danger/safe signals"]
    C --> D["dDCA"]
    D --> E["MCAV, Kalpha and thresholds"]
    E --> F["Join, evaluate and review"]
```

The pipeline must present the antigen before the corresponding signal pair. This order follows the segmented dDCA workflow described by Gu, Greensmith and Aickelin.

## Inputs and outputs

Create the engine with:

- `cells`: number of deterministic cells;
- `antigenCount`: number of valid antigen types, addressed from `0` to `antigenCount - 1`; and
- `maxMigration`: maximum cell migration threshold.

For each input observation, call `doAntigen(integer_id)` followed by `doSignals(danger, safe)`. After all observations, call `results()` to populate:

- `mcav`: mature context antigen value by antigen type;
- `ka`: Kalpha by antigen type;
- `Tk`: calculated Kalpha anomaly threshold; and
- `Tmcav`: calculated MCAV threshold.

Unobserved antigen types produce `NaN` in `mcav` and `ka`. `results()` also assumes at least one signal pair and a non-zero total safe signal. These and other retained implementation assumptions are listed in [review notes](docs/REVIEW_NOTES.md).

## Deliberate divergence from the shared C implementation

The shared C implementation initialises the migration interval as a floating-point value, but declares the equivalent reset interval as an integer. At the commonly used configuration of 100 cells and maximum migration 100, the initial interval is `100 / 99`, while subsequent C resets truncate that interval to `1`.

This change of type between initialisation and reset is believed to be unintentional. dDCA therefore preserves each cell's floating-point maximum lifespan on every reset. This keeps the migration thresholds consistent across cell incarnations, but means results can diverge from the C executable after repeated migrations. Strict reproduction of the C executable would require intentionally restoring the integer truncation.

## Installation

dDCA requires Python 3.10 or newer and NumPy.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Install the optional `stats` dependency to use `dc_stats()`:

```bash
python -m pip install -e '.[stats]'
```

## Usage

```python
from ddca import dDCA

engine = dDCA(cells=100, antigenCount=500, maxMigration=100)

for antigen_id, danger, safe in derived_signal_stream:
    engine.doAntigen(int(antigen_id))
    engine.doSignals(float(danger), float(safe))

engine.results()

mcav = engine.mcav
kalpha = engine.ka
kalpha_threshold = engine.Tk
mcav_threshold = engine.Tmcav
```

The example deliberately begins with already-derived values: input feature engineering is part of the surrounding pipeline, not this module.

## Build with Cython

A C compiler and Python development headers are required. On macOS, install the Xcode Command Line Tools; on Debian/Ubuntu, install `build-essential` and the development package for the selected Python version.

From an activated virtual environment:

```bash
python -m pip install --upgrade pip
python -m pip install -e '.[build]'
python setup.py build_ext --inplace
```

Confirm that Python loads the compiled extension (`.so` on macOS/Linux or `.pyd` on Windows):

```bash
python -c 'import ddca; print(ddca.__version__, ddca.__file__)'
```

The generated `ddca.c`, extension binary, and `build/` directory are build artefacts and are excluded from version control. The authoritative source remains `ddca.py`; the version is exposed separately as `ddca.__version__`. Run the tests against either the pure-Python module or the compiled extension:

```bash
python -m unittest discover -s tests -v
```

## Provenance and acknowledgements

Benjamin Donnachie developed this Python research implementation. Particular thanks are due to **Dr Julie Greensmith for sharing her C implementation**, dated 27 March 2008. Its header records subsequent modification by Feng Gu on 11 July 2008. That implementation informed behavioural compatibility, including cell selection and ordered antigen logging. The C source is not redistributed here and is not relicensed by this repository.

The algorithmic basis should be cited as:

- Julie Greensmith and Uwe Aickelin, “[The Deterministic Dendritic Cell Algorithm](https://doi.org/10.1007/978-3-540-85072-4_26),” ICARIS 2008, LNCS 5132, pp. 291–302.
- Feng Gu, Julie Greensmith and Uwe Aickelin, “[Integrating Real-Time Analysis with the Dendritic Cell Algorithm through Segmentation](https://doi.org/10.1145/1569901.1570063),” GECCO 2009, pp. 1203–1210.

OpenAI Codex, Anthropic Claude, and [OpenWolf](https://github.com/cytostack/openwolf) were used only to prepare this public repository, including its documentation, packaging, provenance audit, and release checks. They were not involved in the original research or implementation.

## Version history

The concrete changes are documented in [CHANGELOG.md](CHANGELOG.md). Locally retained historical and current sources were compared directly; repeated experiment copies were grouped by checksum so they were not mistaken for separate releases. The original private repository was also inspected locally, providing commit anchors from the first v0.1 scaffold in June 2021 through the object-based implementation and four July 2024 updates. Retained snapshots document intermediate experiments that were not committed on the main `ddca.py` path.

## Licence

Copyright © 2021–2026 Benjamin Donnachie.

dDCA is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Licence](https://creativecommons.org/licenses/by-nc-sa/4.0/) (CC BY-NC-SA 4.0). See [LICENSE.md](LICENSE.md).

The cited papers, shared C implementation, third-party tools, and upstream pipeline data are not relicensed by this repository.
