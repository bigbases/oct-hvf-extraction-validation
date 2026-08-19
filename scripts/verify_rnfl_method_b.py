"""
RNFL quadrant 정확도 — Method B (row-eye 조건부 flip 후 공식 GT 대조).

배경 (2026-07-16/17): oct_tabular_90d.csv의 rnfl_q_t/rnfl_q_n은 OS 행에서
OD-normalized 컨벤션으로 기록됨(coordinate_frame_convention.md 참고). 이를
모르고 raw를 GT와 직접 대조하면(scripts/phase3_extraction_accuracy.py의
단순 컬럼대조 방식) OS의 T/N이 구조적으로 틀리게 세어져 정확도가 실제보다
낮게 나온다. 이 스크립트는 OS 행에 한해 rnfl_q_t<->rnfl_q_n을 스왑한 뒤
대조하는 Method B를 구현한다.

입력:
  data/rnfl_detail_raw.csv  (현재 추출, on-screen 컨벤션)
  oct_tabular_90d.csv       (수동검수 GT)

출력:
  results/rnfl_method_b_accuracy.json
"""
import csv
import json
import math
from pathlib import Path
from collections import defaultdict

_ROOT = Path(__file__).resolve().parents[1]

RNFL_COLS = ['rnfl_q_s', 'rnfl_q_t', 'rnfl_q_i', 'rnfl_q_n']


def fnum(v):
    if v in ('', None, 'nan', 'None'):
        return None
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except (TypeError, ValueError):
        return None


def main():
    raw = list(csv.DictReader(open(_ROOT / 'data' / 'rnfl_detail_raw.csv', encoding='utf-8-sig')))
    gt = list(csv.DictReader(open(_ROOT / 'oct_tabular_90d.csv', encoding='utf-8-sig')))
    raw_idx = {(r['patient_id'], r['eye'], r.get('oct_date', '')): r for r in raw}
    gt_idx = {(r['patient_id'], r['eye'], r['oct_date']): r for r in gt}

    n_rows_matched = 0
    n_flagged = 0
    n_valid = n_gt_missing = n_raw_missing = n_exact = 0
    by_col = defaultdict(lambda: [0, 0])
    qi = [0, 0]

    for key, r in raw_idx.items():
        pid, eye, oct_date = key
        gt_row = gt_idx.get((pid, eye, oct_date))
        if not gt_row:
            continue
        n_rows_matched += 1
        if (r.get('cv_flags') or '').strip():
            n_flagged += 1

        vals = dict(r)
        if eye == 'OS':
            vals['rnfl_q_t'], vals['rnfl_q_n'] = vals.get('rnfl_q_n'), vals.get('rnfl_q_t')

        for c in RNFL_COLS:
            gv = fnum(gt_row.get(c))
            rv = fnum(vals.get(c))
            if gv is None:
                n_gt_missing += 1
                continue
            if rv is None:
                n_raw_missing += 1
                continue
            n_valid += 1
            by_col[c][1] += 1
            if c == 'rnfl_q_i':
                qi[1] += 1
            if gv == rv:
                n_exact += 1
                by_col[c][0] += 1
                if c == 'rnfl_q_i':
                    qi[0] += 1

    result = {
        'method': 'Method B: OS 행에서 rnfl_q_t<->rnfl_q_n 스왑 후 공식 GT(oct_tabular_90d.csv) 대조',
        'n_rows_matched_to_gt': n_rows_matched,
        'n_valid': n_valid, 'n_exact': n_exact,
        'exact_match_pct': round(100 * n_exact / max(n_valid, 1), 1),
        'n_gt_missing': n_gt_missing, 'n_raw_missing': n_raw_missing,
        'missing_rate_pct': round(100 * n_raw_missing / max(n_valid + n_raw_missing, 1), 1),
        'per_column': {c: {'exact': ex, 'n': n, 'pct': round(100 * ex / max(n, 1), 1)}
                       for c, (ex, n) in by_col.items()},
        'rnfl_q_i_alone': {'exact': qi[0], 'n': qi[1], 'pct': round(100 * qi[0] / max(qi[1], 1), 1)},
        'cv_flag_rate_pct': round(100 * n_flagged / max(n_rows_matched, 1), 1),
        'n_flagged': n_flagged,
    }

    out_path = _ROOT / 'results' / 'rnfl_method_b_accuracy.json'
    json.dump(result, open(out_path, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)

    print(f'rows matched to GT: {n_rows_matched}')
    print(f'exact-match: {n_exact}/{n_valid} = {result["exact_match_pct"]}%')
    print(f'missing_rate: {result["missing_rate_pct"]}%')
    for c, s in result['per_column'].items():
        print(f'  {c}: {s["exact"]}/{s["n"]} = {s["pct"]}%')
    print(f'rnfl_q_i alone: {result["rnfl_q_i_alone"]["pct"]}%')
    print(f'cv_flag rate: {result["cv_flag_rate_pct"]}%')
    print(f'saved: {out_path}')


if __name__ == '__main__':
    main()
