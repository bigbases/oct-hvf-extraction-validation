# Coordinate frame convention (OD / OS)

Two coordinate frames are used in this project, and which one a given file uses
was for a long time implicit in the code. Comparing data across frames silently
counts correct values as errors. This document is the reference for that
convention; the `scoring_note` fields inside `results/*.json` point here.

## Summary

1. Two frames: **on-screen** (as the device drew it) and **OD-normalized** (both
   eyes expressed in right-eye orientation).
2. Files whose name ends in `_flip` are OD-normalized; the rest are on-screen.
   The exception is the GCA `os_s_*` columns of the ground-truth table, which are
   row-eye dependent (§3).
3. OCR output is always on-screen. Scoring it directly against an OD-normalized
   reference makes every OS temporal/nasal and clock-hour value count as a
   mismatch, so extraction accuracy comes out far lower than it is.

## 1. The two frames

| Frame | Definition | Consequence |
|---|---|---|
| **on-screen** | Exactly as the device (Cirrus, HFA) printed it | On an OD report temporal is on the left; on an OS report temporal is on the right. The two eyes are mirror images |
| **OD-normalized** | Every eye expressed in right-eye orientation (OS mirrored) | The same feature index means the same anatomical location, so ISNT and Garway-Heath topography stay consistent across eyes |

The mirroring is horizontal: temporal and nasal swap, superior and inferior do
not.

## 2. Frame by artefact

| Artefact | Frame | Note |
|---|---|---|
| OCR output | on-screen | As printed. In the GCA six-sector pie chart the temporal label sits at the same screen position for both eyes, because that chart is a schematic overlay rather than a photograph — see `parse_sectors()` |
| `clockhours.csv`, `clockhours_180d.csv` | on-screen | RNFL reference, unflipped |
| `clockhours_flip.csv`, `clockhours_180d_flip.csv` | OD-normalized | Produced by `scripts/extraction/make_rnfl_flip.py` |
| `ml_final_*.csv` (unflipped) | on-screen | Derived from OCR output, **not** an independent reference — see §2.1 |
| `ml_final_*_flip.csv` | OD-normalized | Same lineage. The reported statistics do not use these files |
| `sfa_thresholds.csv` | on-screen | Visual-field thresholds, unflipped |
| `oct_tabular_90d.csv` (manual reference) | mixed, row-eye dependent | `od_s_*` is always on-screen; the four `os_s_{sup_t,sup_n,inf_t,inf_n}` columns depend on the row's eye — see §3 |
| `data/analysis_master.csv` | mixed, inherited from the reference | The file the reported correlations are computed from. It inherits the row-eye dependence above, but the analysis code only ever reads the columns matching each row's own eye, so the results are unaffected (§3) |

**Rule of thumb:** `_flip` in the filename means OD-normalized; no `_flip` means
on-screen; OCR output is always on-screen.

### 2.1 Why `ml_final_*` is not a validation reference

Compared column by column against the OCR output it was built from, the GCA
sector columns of `ml_final_*` agree 68–99% of the time. It is a derivative of
that output, not an independent transcription. Scoring the OCR against it would
be scoring the extraction against itself. Extraction accuracy is therefore
measured against the manually transcribed reference (`oct_tabular_90d.csv`, and
`data/analysis_master.csv` which inherits it), where agreement is markedly
lower.

## 3. The row-eye dependence in the reference table

In `oct_tabular_90d.csv` the same report, on the same date, is recorded
differently depending on whether the row's `eye` is OD or OS:

| Row eye | `os_s_sup_t` | `os_s_sup_n` |
|---|---|---|
| OD | 20 | 45 |
| OS | 45 | 20 |

Checked against the optic-disc position in the thickness map, the anatomically
correct value is 45: on an OS report the disc sits on the screen-left side, so
screen-right is temporal. In other words:

- on OS rows, `os_s_*` is on-screen (and matches the anatomy directly);
- on OD rows, `os_s_*` is OD-normalized (swapped);
- `od_s_*` shows no such dependence — only those four `os_s_*` columns are
  affected.

**Why this matters.** Scoring on-screen OCR output against this table without
accounting for the dependence mixes two conventions and lands near the average
of the two, around 50%. Applying the swap conditionally on the row's eye
recovers the true figure. The before/after comparison is recorded in
`results/gca_row_eye_scoring_comparison.json`; the check is reproduced by
`scripts/verify_gca_os_mirroring.py`.

**Why the reported analysis is unaffected.** `scripts/phase5_2_sector_vf.py` and
`scripts/phase5_1_correlations.py` build the column name from the row itself
(`e = row['eye'].lower(); col = f'{e}_s_{s}'`), so an OD row never reads an
`os_s_*` column. The mixed convention is never exercised.

## 4. The normalization rule

Applied to OS rows only when producing `*_flip` files; values are preserved and
positions swapped. OD rows are left untouched.

### 4.1 RNFL

- Quadrants: `rnfl_q_t ↔ rnfl_q_n`. `rnfl_q_s` and `rnfl_q_i` are fixed.
- Clock hours: `H12` and `H06` are fixed points on the vertical axis;
  `H01↔H11`, `H02↔H10`, `H03↔H09`, `H04↔H08`, `H05↔H07`.

### 4.2 GCA sectors

- `os_s_sup_t ↔ os_s_sup_n`, `os_s_inf_t ↔ os_s_inf_n`. Superior and inferior
  are fixed.
- The sector names encode direction, so a Methods description written in
  on-screen terms would contradict the stored data, which is OD-normalized.

### 4.3 Visual field, 52 points

- Reversed within each StatPac row (`p01–p04`, `p05–p10`, …, `p51–p54`).
- Row order is preserved, since the mirroring is horizontal.

## 5. Was the normalization applied to every eye?

Checked directly against the data — the flip scripts reverse every OS row
unconditionally:

| Data | OS eyes | Swap confirmed | Undetectable | OD rows changed |
|---|---:|---:|---:|---:|
| RNFL clock hours | 142 | 140 | 2 (missing values) | 0 |
| GCA sectors | 139 | 129 | remainder (t equals n, or missing) | few* |
| VF row 1 (`p01–p04`) | 139 | 110 | 22 (symmetric, so indistinguishable) | — |

Every OS eye carrying a value was normalized; partial normalization is ruled
out. The undetectable cases are those where the swap cannot be observed because
the two values are equal or missing, not evidence of a missed row.

\* The few changed OD rows come from individual manual corrections, unrelated to
the coordinate frame.

## 6. Measuring extraction accuracy correctly

Score on-screen OCR output against an on-screen reference. Either:

**A. Keep the reference on-screen** — compare against `clockhours.csv`
(unflipped) for RNFL, or `sfa_thresholds.csv` for the visual field.

**B. Normalize the OCR output first** — apply the §4 swap to OS rows, then
compare against the `*_flip` reference.

The two must agree. What must *not* be done is comparing on-screen OCR output
directly against an OD-normalized reference: every OS temporal/nasal and clock
hour is then structurally counted as wrong.

## 7. Checklist before a new comparison

- Which frame is each file in? (§2, and the `_flip` naming rule.)
- When merging or scoring two sources, are they in the same frame — or was one
  flipped first?
- Does the Methods description of the coordinate frame match the frame the data
  is actually in? (For GCA, the data is OD-normalized.)
- Before concluding that a reference is wrong, has a frame mismatch been ruled
  out?
