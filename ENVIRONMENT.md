# Execution environment

The reported numbers were produced with the versions below. OCR output is tied
to the Tesseract build, so a different version can change extracted values.

| Component | Version |
|---|---|
| Python | 3.14.3 |
| Tesseract OCR | 5.5.0.20241111 (leptonica 1.85.0) |
| Platform | Windows 11 |

Exact package pins are in `requirements.txt` (numpy 2.4.3, pandas 3.0.2,
scipy 1.17.1, scikit-learn 1.8.0, xgboost 3.2.0, pytesseract 0.3.13,
Pillow 12.1.1, matplotlib 3.10.8, PyYAML 6.0.3, openpyxl 3.1.5).

## Install

```bash
python -m venv .venv && source .venv/bin/activate   # .venv/Scripts/activate on Windows
pip install -r requirements.txt      # analysis and OCR
pip install -r requirements-dev.txt  # + tests
pip install -e .                     # the hvf package (src layout)
make check                           # config validation + tests
```

Tesseract is installed separately; point `pytesseract` at it through
`config/params.yaml` (`paths.tesseract`) or the `HVF_TESSERACT_CMD` environment
variable if it is not on `PATH`.

## Checking the versions

```bash
python --version
tesseract --version
pip freeze
```

Model training was run on a separate machine and is not included in this
repository; it is not part of the reported analysis.
