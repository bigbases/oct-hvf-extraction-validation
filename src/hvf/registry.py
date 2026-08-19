"""Results registry — 논문에 들어갈 모든 수치의 단일 기록소.

원칙: 논문/발표의 숫자는 손으로 적지 않는다. 계산 함수가 record_result 로
      {값, 계산 스크립트, 입력파일 해시, 생성시각, git commit} 과 함께 JSON 을 남긴다.
      results/<key>.json 하나 = 검증 가능한 수치 하나.

의존: 표준 라이브러리만 (numpy/pandas 불필요) → 테스트가 가볍다.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from .config import project_root, results_dir


def sha256_file(path: str | Path, _chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(_chunk), b""):
            h.update(block)
    return h.hexdigest()


def _git_commit() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=project_root(), capture_output=True, text=True, timeout=5,
        )
        return out.stdout.strip() or None
    except Exception:
        return None


def _describe_input(path: str | Path) -> dict:
    """입력 파일 1개를 {식별자, 해시, 크기}로 기록.

    PHI/사용자명 누출 방지: 레포 루트 밑이면 상대경로, 밖(원본 데이터)이면 파일명만.
    """
    p = Path(path)
    try:
        ident = str(p.resolve().relative_to(project_root())).replace("\\", "/")
    except ValueError:
        ident = p.name  # 레포 밖(raw) → 이름만, 전체 경로는 남기지 않음
    entry: dict[str, Any] = {"id": ident}
    if p.exists():
        entry["sha256"] = sha256_file(p)
        entry["bytes"] = p.stat().st_size
    else:
        entry["sha256"] = None
        entry["bytes"] = None
        entry["missing"] = True
    return entry


def record_result(
    key: str,
    value: Any,
    *,
    script: str,
    inputs: Iterable[str | Path] | None = None,
    extra: dict | None = None,
) -> Path:
    """수치 하나를 results/<key>.json 으로 기록하고 경로를 반환.

    Parameters
    ----------
    key    : 파일명이자 고유 식별자 (예: "pairs_count", "sf_pearson_od").
    value  : JSON 직렬화 가능한 값 (숫자/dict/list). 계산 결과 그 자체.
    script : 이 값을 계산한 스크립트 (레포 상대경로 권장).
    inputs : 계산에 쓴 입력 파일들 → 해시로 봉인. 입력이 바뀌면 해시가 달라짐.
    extra  : 부가 메타(포함/제외 기준, 보정변수, n 등)를 자유롭게.
    """
    entry = {
        "key": key,
        "value": value,
        "script": script,
        "inputs": [_describe_input(p) for p in (inputs or [])],
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": _git_commit(),
        "python": platform.python_version(),
        "extra": extra or {},
    }
    out = results_dir() / f"{key}.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(entry, f, ensure_ascii=False, indent=2, sort_keys=False)
    return out


def load_result(key: str) -> dict:
    with open(results_dir() / f"{key}.json", encoding="utf-8") as f:
        return json.load(f)


def list_results() -> list[str]:
    return sorted(p.stem for p in results_dir().glob("*.json"))
