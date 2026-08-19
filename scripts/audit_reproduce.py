# -*- coding: utf-8 -*-
"""재현 격차 감사 — '재현됨' 항목만 계산해 results/ registry 에 봉인.

부분재현·재현안됨 항목은 여기서 잠그지 않는다(격차 지도 문서에만 기록).
이 스크립트는 현재 데이터로부터 값을 '있는 그대로' 계산한다. KCC 숫자에 맞추지 않는다.

실행: PYTHONPATH=src python scripts/audit_reproduce.py
입력: kcc/gca_vf_matched.csv, zeiss_gca_ocr.csv  (PHI — 내용은 registry 에 안 남고 해시만)
"""
from __future__ import annotations
import csv
from pathlib import Path

import numpy as np
from sklearn.linear_model import LinearRegression

from hvf.config import project_root, set_seed
from hvf.registry import record_result

ROOT = project_root()
SF_CSV = ROOT / "kcc" / "gca_vf_matched.csv"
GCA_CSV = ROOT / "zeiss_gca_ocr.csv"


def _target_gcl(r: dict) -> float:
    return float(r["avg_gcl_od"] if r["eye"] == "OD" else r["avg_gcl_os"])


def audit_structure_function() -> None:
    rows = list(csv.DictReader(open(SF_CSV, encoding="utf-8-sig")))
    n_od = sum(1 for r in rows if r["eye"] == "OD")
    n_os = sum(1 for r in rows if r["eye"] == "OS")
    # 완전 케이스(ms·두께·나이·성별 파싱 가능)
    cc = []
    for r in rows:
        try:
            cc.append((r["eye"], float(r["ms"]), _target_gcl(r),
                       float(r["age"]), r["gender"].strip().upper()))
        except (ValueError, KeyError):
            pass

    record_result(
        "sf_sample_counts",
        {"n_matched_rows": len(rows), "n_od": n_od, "n_os": n_os,
         "n_complete_case": len(cc)},
        script="scripts/audit_reproduce.py", inputs=[SF_CSV],
        extra={"note": "n_matched_rows=파일행수(양안 매칭). 상관/회귀는 complete-case 로 계산.",
               "eye_thickness_rule": "target eye 의 avg_gcl (OD→avg_gcl_od, OS→avg_gcl_os)"},
    )

    pearson = {}
    r2_simple = {}
    for tag in ("OD", "OS", "ALL"):
        sub = [d for d in cc if tag == "ALL" or d[0] == tag]
        ms = np.array([d[1] for d in sub]); th = np.array([d[2] for d in sub])
        pearson[tag] = round(float(np.corrcoef(th, ms)[0, 1]), 3)
        r2_simple[tag] = round(float(
            LinearRegression().fit(th.reshape(-1, 1), ms).score(th.reshape(-1, 1), ms)), 3)

    record_result(
        "sf_pearson_r",
        {"n": len(cc), "by_eye": {"OD": len([d for d in cc if d[0]=='OD']),
                                  "OS": len([d for d in cc if d[0]=='OS'])},
         "r": pearson},
        script="scripts/audit_reproduce.py", inputs=[SF_CSV],
        extra={"x": "GCL+IPL 평균두께(target eye)", "y": "VF Mean Sensitivity",
               "method": "numpy.corrcoef, complete-case, no imputation"},
    )
    record_result(
        "sf_regression_simple_r2",
        {"r2": r2_simple},
        script="scripts/audit_reproduce.py", inputs=[SF_CSV],
        extra={"model": "MS ~ thickness (단순)", "impl": "sklearn LinearRegression.score"},
    )


def audit_gca_extraction() -> None:
    rows = list(csv.DictReader(open(GCA_CSV, encoding="utf-8-sig")))
    ok = sum(1 for r in rows if r.get("patient_id", "").strip())
    record_result(
        "gca_extraction",
        {"total": len(rows), "id_extracted": ok,
         "success_rate_pct": round(ok / len(rows) * 100, 1)},
        script="scripts/audit_reproduce.py", inputs=[GCA_CSV],
        extra={"success_def": "raw 헤더에서 patient_id 파싱 성공"},
    )


def main() -> None:
    set_seed()
    audit_structure_function()
    audit_gca_extraction()
    print("[audit] 재현됨 항목 registry 봉인 완료 → results/")


if __name__ == "__main__":
    main()
