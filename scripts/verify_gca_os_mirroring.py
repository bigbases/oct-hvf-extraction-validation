"""
GCA 6섹터 파이차트 추출 검증 + GT의 row-eye 의존적 os_s_* 컨벤션 문서화.

배경 (2026-07-16 하루 종일 검증, 결론 확정):
  - GCA 6섹터 파이차트는 스키마틱 오버레이로, T 라벨이 양안 모두 화면상
    동일 위치(30도 등)에 인쇄됨. B-scan/두께맵(사진)은 실제로 좌우 반전되어
    렌더링되지만(시신경유두 위치로 확인), 파이차트는 그와 무관하게 고정 배치.
    -> src/hvf/ocr_oct.py::parse_sectors()는 eye 분기 없이 동일 각도를 씀(정답).
  - GT(oct_tabular_90d.csv)는 os_s_* 컬럼을 **row의 eye가 무엇이냐에 따라
    다르게** 기록한다: 같은 리포트·같은 날짜인데 eye=OD인 행과 eye=OS인 행의
    os_s_sup_t 값이 서로 스왑되어 있다(대표 리포트 1건에서 직접 확인:
    OD행 os_s_sup_t=20, OS행 os_s_sup_t=45 — 실제 파이차트 값은 45).
    즉 GT는 "OD행에서는 os_s_*를 OD-normalized로, OS행에서는 on-screen
    그대로"로 기록되어 있다. 원인(의도적 정규화 vs 수기 검수 관행)은 불명이나
    현재 데이터의 재현 가능한 사실이다.

이 스크립트는 이 GT 특성을 반영한 올바른 채점 로직을 구현하고, 향후 동일한
혼란이 재발하지 않도록 **반드시 git에 커밋**한다.

입력:
  data/oct_values_raw.csv   (identical-angle 코드로 생성된 최신 추출)
  oct_tabular_90d.csv       (수동검수 GT)

출력:
  results/gca_os_mirroring_verification.json
"""
import csv
import json
import math
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]

SWAP_PAIRS = [('os_s_sup_t', 'os_s_sup_n'), ('os_s_inf_t', 'os_s_inf_n')]
ALL_SECTOR_COLS = {
    'OD': ['od_s_sup', 'od_s_sup_t', 'od_s_inf_t', 'od_s_inf', 'od_s_inf_n', 'od_s_sup_n'],
    'OS': ['os_s_sup', 'os_s_sup_t', 'os_s_inf_t', 'os_s_inf', 'os_s_inf_n', 'os_s_sup_n'],
}


def fnum(v):
    if v in ('', None, 'nan', 'None'):
        return None
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except (TypeError, ValueError):
        return None


def swapped_col(c):
    for a, b in SWAP_PAIRS:
        if c == a:
            return b
        if c == b:
            return a
    return c


def score_column(raw_idx, gt_idx, rows, col, swap):
    n = 0
    ex = 0
    src_col = swapped_col(col) if swap else col
    for k in rows:
        rv = fnum(raw_idx[k].get(src_col))
        gv = fnum(gt_idx[k].get(col))
        if rv is None or gv is None:
            continue
        n += 1
        if rv == gv:
            ex += 1
    return n, ex


def main():
    raw = list(csv.DictReader(open(_ROOT / 'data' / 'oct_values_raw.csv', encoding='utf-8-sig')))
    gt = list(csv.DictReader(open(_ROOT / 'oct_tabular_90d.csv', encoding='utf-8-sig')))
    raw_idx = {(r['patient_id'], r['eye']): r for r in raw}
    gt_idx = {(r['patient_id'], r['eye']): r for r in gt}
    common = sorted(set(raw_idx) & set(gt_idx))
    od_rows = [k for k in common if k[1] == 'OD']
    os_rows = [k for k in common if k[1] == 'OS']

    result = {'n_common': len(common), 'n_od_rows': len(od_rows), 'n_os_rows': len(os_rows)}

    # ---- 어느 raw 추출 컨벤션인지 자동판별 ----
    # (raw가 eye-dependent/mirrored 코드 산출이면 OD행=as-is/OS행=swap이 맞고,
    #  identical-angle 코드 산출이면 그 반대가 맞음. 매번 가정하지 않고 둘 다
    #  계산해서 더 잘 맞는 쪽을 채택 — 코드 변경에 안전하게.)
    swap_cols = {a for pair in SWAP_PAIRS for a in pair}

    def total_for(od_swap, os_swap):
        n = e = 0
        for c in swap_cols:
            n1, e1 = score_column(raw_idx, gt_idx, od_rows, c, swap=od_swap)
            n2, e2 = score_column(raw_idx, gt_idx, os_rows, c, swap=os_swap)
            n += n1 + n2
            e += e1 + e2
        return n, e

    n_a, e_a = total_for(od_swap=False, os_swap=True)   # eye-dependent/mirrored raw 가정
    n_b, e_b = total_for(od_swap=True, os_swap=False)   # identical-angle raw 가정
    pct_a = 100 * e_a / max(n_a, 1)
    pct_b = 100 * e_b / max(n_b, 1)
    od_swap, os_swap = (False, True) if pct_a >= pct_b else (True, False)
    result['convention_detected'] = (
        'raw=eye-dependent(mirrored) -> OD행 as-is/OS행 swap' if pct_a >= pct_b
        else 'raw=identical-angle -> OD행 swap/OS행 as-is'
    )
    result['convention_A_pct'] = round(pct_a, 1)
    result['convention_B_pct'] = round(pct_b, 1)

    per_col = {}
    tot_n = 0
    tot_ex = 0
    for eye_key, cols in ALL_SECTOR_COLS.items():
        for c in cols:
            if c in swap_cols:
                n1, e1 = score_column(raw_idx, gt_idx, od_rows, c, swap=od_swap)
                n2, e2 = score_column(raw_idx, gt_idx, os_rows, c, swap=os_swap)
                n, e = n1 + n2, e1 + e2
            else:
                n, e = score_column(raw_idx, gt_idx, common, c, swap=False)
            per_col[c] = {'n': n, 'exact': e, 'pct': round(100 * e / max(n, 1), 1)}
            tot_n += n
            tot_ex += e

    result['per_column_all_6_sectors'] = per_col
    result['overall_mean_pct'] = round(100 * tot_ex / max(tot_n, 1), 1)

    # ---- 참고용: row-eye 조건부 채점을 안 썼을 때(단순 직접비교, 공식스크립트 방식) ----
    naive_n = naive_ex = 0
    for c in swap_cols:
        n, e = score_column(raw_idx, gt_idx, common, c, swap=False)
        naive_n += n
        naive_ex += e
    result['naive_no_row_condition_pct_for_swap_cols'] = round(100 * naive_ex / max(naive_n, 1), 1)

    out_path = _ROOT / 'results' / 'gca_os_mirroring_verification.json'
    json.dump(result, open(out_path, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)

    print(f'n_common={len(common)} (OD행 {len(od_rows)}, OS행 {len(os_rows)})')
    print('\n6섹터 전체 (row-eye 조건부 채점):')
    for c, s in per_col.items():
        print(f'  {c:14s}: {s["exact"]}/{s["n"]} = {s["pct"]}%')
    print(f'\n전체 평균: {result["overall_mean_pct"]}%')
    print(f'(참고) row-eye 조건 무시한 단순비교: {result["naive_no_row_condition_pct_for_swap_cols"]}%')
    print(f'\nsaved: {out_path}')


if __name__ == '__main__':
    main()
