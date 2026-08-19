"""hvf — 녹내장 OCT–VF 재현 파이프라인 공용 라이브러리.

핵심 진입점:
  from hvf.config   import load_config, get, seed, set_seed, data_dir, results_dir, configure_tesseract
  from hvf.registry import record_result, load_result, sha256_file
"""

__version__ = "0.1.0"
