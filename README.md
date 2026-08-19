# HVF — OCR extraction pipeline for paired OCT and Humphrey visual field reports

Research code accompanying a study on structure–function correlation in glaucoma.
The pipeline reads clinical report images (Cirrus OCT and HFA single-field
analysis printouts), extracts the numeric fields by OCR, validates them against a
manually transcribed reference, and produces the correlation analyses reported in
the manuscript.

**Research code, provided as-is.** It was written for one cohort on one pair of
devices and is not a general-purpose tool.

## What this repository contains

| Path | Contents |
|---|---|
| `src/hvf/` | Library modules — config loading, result registry, cohort filter, OCR routines for VF / RNFL / GCIPL |
| `scripts/` | Pipeline entry points (`stage_00` … `stage_30`) and the analysis scripts that produce every reported number |
| `config/params.yaml` | Single source of truth for constants: crop geometry, seeds, matching windows, plausibility ranges |
| `results/` | Result registry (JSON). Aggregate values only — the numbers cited in the manuscript |
| `paper/results_frozen/` | Frozen copies of the registry as used for the submitted version |
| `paper/analysis_scripts/` | Frozen copies of the analysis scripts matching those results |
| `docs/coordinate_frame_convention.md` | Laterality and mirroring convention (on-screen vs OD-normalized frames). Referenced by the scoring notes inside the result JSONs |
| `tests/` | Smoke tests for config and registry |
| `scripts/manuscript_tables/` | Cohort definition, manuscript tables, and regeneration scripts, with their text/JSON outputs. `scripts/manuscript_tables/_cohort_D_freeze.py` is the sole producer of `cohort_final_D.csv`, which defines the analysis cohort; `src/hvf/cohort.py` refuses to run without it |
| `scripts/extraction/` | Original OCR and laterality-flip implementations that the `src/hvf/` modules were derived from |
| `scripts/manual_review/` | Manual visual-field review pass (Phase B) that produced the corrected threshold table |

## What this repository does not contain

- **Raw imaging data.** No report images, no extracted crops, no per-eye tables.
  Withheld under the governing IRB approval. Everything in `results/` is
  aggregate.
- **Any patient identifier.** No patient IDs, names, dates of birth, or
  examination dates. Date columns in derived tables are month-granularity at
  most, and those tables are not distributed either.
- **The manuscript.** LaTeX sources, the compiled PDF, and the figure files are
  not here. The typeset article is distributed by the journal. This repository
  exists to make the *computation* reproducible, not to redistribute the paper.
- **Development history.** The repository begins at a single commit; the
  day-to-day working history is not part of the release.

Because the raw inputs are withheld, the extraction stages cannot be re-executed
from a clone. What can be checked is the analysis layer: the scripts, the
constants they read, and the frozen registry they wrote.

### Expected behaviour on a fresh clone

`make check` passes: the configuration validates and the smoke tests run. The
analysis scripts, however, stop immediately with `FileNotFoundError` on a path
under `data/` — for example `data/analysis_master.csv` for the correlation
scripts, or `cohort_final_D.csv` for the cohort definition. **This is expected.**
Those inputs hold per-eye rows and are withheld; the scripts deliberately fail
rather than fall back to a default, so no number can be produced from absent
data. To inspect what they would have written, read `results/` and
`paper/results_frozen/`.

## Reproduced values and their producers

| Reported quantity | Script | Output |
|---|---|---|
| VF / GCIPL / RNFL field-level extraction accuracy, error taxonomy | `scripts/phase3_extraction_accuracy.py` | `results/phase3_extraction_accuracy.json` |
| Spatial specificity of VF point errors (χ²) | `scripts/phase3_vf_spatial_chi2.py` | `results/phase3_vf_spatial_chi2.json` |
| Global structure–function correlations | `scripts/phase5_1_correlations.py` | `results/phase5_1_correlations.json` |
| Sector-level structure–function correlations | `scripts/phase5_2_sector_vf.py` | `results/phase5_2_sector_correlations.json` |
| Covariate-adjusted models | `scripts/phase5_3_covariate.py` | `results/phase5_3_covariate.json` |
| GCIPL vs RNFL joint model | `scripts/phase5_4_gcl_rnfl_compare.py` | `results/phase5_4_gcl_rnfl_compare.json` |
| SITA strategy sensitivity | `scripts/phase5_sita_sensitivity.py` | `results/sita_sensitivity.json` |
| Cohort definition (n = 187 eyes) | `scripts/manuscript_tables/_cohort_D_freeze.py` | `cohort_final_D.csv` (not distributed — derived from withheld data) |
| Manuscript tables and cohort curve | `scripts/manuscript_tables/_cohort_D_final_tables.py`, `scripts/manuscript_tables/_manuscript_numbers_freeze.py` | `results/cohortcurve.json` |
| Figures | `scripts/make_sf_matrix.py`, `scripts/make_autovsgt_forest.py`, `paper/analysis_scripts/make_oct_field_accuracy.py` | figure PDFs (not distributed) |

`results/README.md` documents the registry format.

## Requirements

Python 3.14.3 and Tesseract OCR 5.5.0.20241111 (leptonica 1.85.0). Pinned
versions are in `requirements.txt`; `ENVIRONMENT.md` records the environment the
reported numbers were produced in. Model training was run on a separate
machine and is not included in this repository; it is not part of the
reported analysis.

```bash
pip install -r requirements-dev.txt && pip install -e .
make check     # config validation + tests
```

Paths are derived from the repository location. Set `HVF_ROOT` to override, and
`HVF_ZEISS_ARCHIVE` / `HVF_RAW_OCT_ARCHIVE` / `HVF_RAW_VF_ARCHIVE` if you are
supplying your own source archives.

**Inline comments are in Korean.** They record why a particular coordinate
convention, crop window, cohort filter, or correction was chosen, and are kept in
the language they were written in rather than translated after the fact. The
parts a reader needs in order to use the repository are in English: this file,
`ENVIRONMENT.md`, the coordinate-frame convention under `docs/`, and the registry
format in `results/README.md`.

**The pins are deliberate and are not kept current.** They record the environment
the reported numbers were produced in, so dependency scanners will flag advisories
against them. Upgrading is expected to change results: OCR output is tied to the
Tesseract build, and the numerical libraries affect the reported statistics.
Anyone reusing this pipeline on new data should upgrade and re-validate; anyone
reproducing the reported values should install exactly these versions.

## Design rules the code follows

1. No identifying data in the repository; derived tables holding per-eye rows stay local.
2. Deterministic execution: seeds, library versions, and the Tesseract version are pinned.
3. No hand-copied numbers — every reported value is written to `results/` through `hvf.registry`.
4. No hard-coded constants — they live in `config/params.yaml`.

## Verification

`scripts/verify_no_phi.py` scans for identifiers across file names, the working
tree, binary containers (xlsx/docx/pdf), `HEAD`, and the full history, and
separately checks the prose for internal status notes. It proves its own
patterns fire on a synthetic probe before reporting anything, and exits
non-zero if a layer is not clean.

## License

MIT — see `LICENSE`.
