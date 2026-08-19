"""전역 설정 로더 + 결정론 유틸.

config/params.yaml 을 단일 출처로 읽는다. 상수·경로·seed 는 전부 여기 경유.
무거운 의존(numpy/torch/pytesseract)은 필요할 때만 지연 import.
"""
from __future__ import annotations

import os
import random
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml


def project_root() -> Path:
    """레포 루트 (src/hvf/config.py → parents[2])."""
    return Path(__file__).resolve().parents[2]


@lru_cache(maxsize=None)
def load_config(path: str | os.PathLike | None = None) -> dict:
    cfg_path = Path(path) if path else project_root() / "config" / "params.yaml"
    with open(cfg_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def get(*keys: str, default: Any = None) -> Any:
    """중첩 키 접근. 예: get('constants', 'vf', 'max_db')."""
    node: Any = load_config()
    for k in keys:
        if not isinstance(node, dict) or k not in node:
            return default
        node = node[k]
    return node


# ── 경로 헬퍼 ────────────────────────────────────────────────────────────────
def _resolve(p: str | os.PathLike) -> Path:
    p = Path(p)
    return p if p.is_absolute() else (project_root() / p)


def data_dir() -> Path:
    d = _resolve(get("paths", "data_dir", default="data"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def results_dir() -> Path:
    d = _resolve(get("paths", "results_dir", default="results"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def raw_path(which: str) -> Path:
    """which ∈ {vf, oct, zeiss}. 원본 PHI 경로 (존재 보장 안 함)."""
    return Path(get("paths", f"raw_{which}"))


# ── 결정론 ──────────────────────────────────────────────────────────────────
def seed() -> int:
    return int(get("seed", default=42))


def set_seed(s: int | None = None) -> int:
    """random / numpy / torch / PYTHONHASHSEED 전부 고정. 설치 안 된 건 조용히 건너뜀."""
    s = seed() if s is None else s
    os.environ["PYTHONHASHSEED"] = str(s)
    random.seed(s)
    try:
        import numpy as np
        np.random.seed(s)
    except ImportError:
        pass
    try:
        import torch
        torch.manual_seed(s)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(s)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except ImportError:
        pass
    return s


def configure_tesseract() -> str:
    """pytesseract 실행경로를 설정하고 그 경로를 반환.

    우선순위: 환경변수 HVF_TESSERACT_CMD > config paths.tesseract.
    (Docker/서버에서 Windows 경로를 덮어쓰기 위함)
    """
    import pytesseract
    cmd = os.environ.get("HVF_TESSERACT_CMD") or get("paths", "tesseract")
    if cmd:
        pytesseract.pytesseract.tesseract_cmd = cmd
    return cmd
