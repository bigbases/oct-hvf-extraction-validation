"""분석 코호트(cohort D, 187눈) 필터 — 논문 수치의 단일 출처.

`data/analysis_master.csv` 는 매칭 성공한 275눈 전부를 담는 상위 집합이고,
논문에 보고되는 값은 그중 수동검수 참조표준이 확보된 187눈(= `cohort_final_D.csv`)
에서만 산출한다. 두 집합을 혼동하면 같은 입력 파일로도 전혀 다른 계수가 나온다
(예: rnfl_ch_mean β = 0.546(187) vs 0.508(275)).

2026-08-15 재현성 감사에서 `phase5_3_covariate.py` / `phase5_4_gcl_rnfl_compare.py`
가 이 필터 없이 돌아가 봉인된 results/*.json 을 재현하지 못하는 것이 확인되어
공용 헬퍼로 분리했다. 새 분석 스크립트는 반드시 이 모듈을 경유할 것.

주의: `scripts/make_autovsgt_forest.py` 는 동등한 필터를 인라인으로 갖고 있다
(Fig 4 산출물이 바이트 단위로 검증돼 있어 건드리지 않음). 로직을 바꿀 일이
생기면 그쪽도 같이 고쳐야 한다.
"""
from __future__ import annotations

import csv
from pathlib import Path

from hvf.config import project_root

COHORT_FILE = "cohort_final_D.csv"
_HEADER = ["patient_id", "eye", "vf_date"]


def load_cohort_keys(path: str | Path | None = None) -> set[tuple[str, str, str]]:
    """cohort_final_D.csv → {(patient_id, eye, vf_date)} 집합.

    파일이 없으면 예외를 던진다. 조용히 전체(275눈)로 넘어가면 논문 값과 다른
    숫자가 아무 경고 없이 나오므로, 실패는 반드시 시끄러워야 한다.
    """
    p = Path(path) if path else project_root() / COHORT_FILE
    if not p.exists():
        raise FileNotFoundError(
            f"분석 코호트 정의 파일을 찾을 수 없음: {p}\n"
            "논문 수치는 187눈 코호트에서만 재현된다 — 전체 275눈으로 대체 실행하지 말 것."
        )
    keys = set()
    with open(p, encoding="utf-8-sig") as f:
        for row in csv.reader(f):
            if row == _HEADER or not row:
                continue
            keys.add(tuple(row[:3]))
    return keys


def filter_rows(rows, path: str | Path | None = None):
    """analysis_master.csv 의 DictReader 행들을 코호트 D로 거른다.

    (patient_id, eye, vf_date) 3-튜플 기준. analysis_master 에서 (patient_id, eye)
    가 유일하므로 make_autovsgt_forest.py 의 (pid, eye)→vf_date 역참조 방식과
    결과가 동일하다(2026-08-15 실측: 양쪽 다 187행, 누락 0).
    """
    keys = load_cohort_keys(path)
    kept = [r for r in rows
            if (r["patient_id"], r["eye"], r["vf_date"]) in keys]
    missing = len(keys) - len(kept)
    if missing:
        raise ValueError(
            f"코호트 정의 {len(keys)}행 중 {missing}행이 analysis_master 에 없음 — "
            "입력 정본이 어긋났는지 확인할 것."
        )
    return kept
