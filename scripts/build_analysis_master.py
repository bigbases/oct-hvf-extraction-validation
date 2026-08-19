"""
scripts/build_analysis_master.py
STEP 5 — data/analysis_master.csv (데이터 준비 최종 산출물)

조인:
  vf_results.csv (276행, anchor)
  x oct_canonical.csv (276행)  JOIN (patient_id, eye)  — 1:1
  LEFT JOIN demographics.csv (144명)  JOIN patient_id    — age/sex 결측 허용

age_at_oct = (oct_date - dob) / 365.25  소수점 1자리, oct_date 기준

covariate_set flag:
  'full'  — dob+sex 모두 있음 (age+sex 보정 가능, 183눈 예상)
  'age'   — dob만 있음       (age만 보정 가능)
  'none'  — dob 없음         (주 분석만)

출력: data/analysis_master.csv (276행)
봉인: results/step5_provenance.json
"""
import csv, json, hashlib, datetime, sys
from pathlib import Path
from datetime import datetime as dt

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))

from hvf.config import data_dir, results_dir
from hvf.registry import sha256_file

# ---------------------------------------------------------------------------
# 컬럼 정의
# ---------------------------------------------------------------------------
_META_COLS = ['patient_id', 'eye', 'vf_date', 'oct_date', 'gap_days']
_DEMO_COLS = ['age_at_oct', 'sex']
_VF_COLS   = (['vf_mean', 'n_absolute_scotoma'] +
               [f'p{i}' for i in range(1, 55)])  # p1..p54
_OCT_COLS  = [
    'avg_gcl_od', 'avg_gcl_os', 'min_gcl_od', 'min_gcl_os',
    'od_s_sup', 'od_s_sup_t', 'od_s_inf_t', 'od_s_inf', 'od_s_inf_n', 'od_s_sup_n',
    'os_s_sup', 'os_s_sup_t', 'os_s_inf_t', 'os_s_inf', 'os_s_inf_n', 'os_s_sup_n',
    'rnfl_q_s', 'rnfl_q_t', 'rnfl_q_i', 'rnfl_q_n',
    'rnfl_h01', 'rnfl_h02', 'rnfl_h03', 'rnfl_h04', 'rnfl_h05', 'rnfl_h06',
    'rnfl_h07', 'rnfl_h08', 'rnfl_h09', 'rnfl_h10', 'rnfl_h11', 'rnfl_h12',
    'ss_gca_od', 'ss_gca_os', 'ss_rnfl_od', 'ss_rnfl_os',
]
_FLAG_COLS = ['covariate_set', 'severity_tertile']

OUT_COLS = _META_COLS + _DEMO_COLS + _VF_COLS + _OCT_COLS + _FLAG_COLS


def _assign_severity_tertiles(rows_out):
    """
    severity_tertile 컬럼 계산 (2026-07-20 도입).

    정의: vf_mean(MS)의 33rd/67th 백분위수(numpy 기본 선형보간, 미반올림)를
    경계로 mild(>q67) / moderate(q33 초과 & q67 이하) / severe(q33 이하)로
    3분할. 논문 Methods §2.1의 "mild (>27.1 dB), moderate (21.6–27.1 dB),
    severe (<=21.6 dB)"는 이 정밀 경계의 반올림 표기일 뿐, 계산 기준은
    항상 이 함수의 미반올림 백분위수다(경계 근처 안이 반올림 경계로는
    다르게 분류될 수 있음 — 2026-07-20 감사에서 발견).
    """
    ms = np.array([float(r['vf_mean']) for r in rows_out if r.get('vf_mean') not in ('', None)])
    q33 = float(np.percentile(ms, 33))
    q67 = float(np.percentile(ms, 67))
    for r in rows_out:
        v = r.get('vf_mean')
        if v in ('', None):
            r['severity_tertile'] = ''
            continue
        v = float(v)
        if v > q67:
            r['severity_tertile'] = 'mild'
        elif v <= q33:
            r['severity_tertile'] = 'severe'
        else:
            r['severity_tertile'] = 'moderate'
    return q33, q67


def _read_csv(path):
    return list(csv.DictReader(open(path, encoding='utf-8-sig')))


def _age_at(oct_date_str, dob_str):
    """age_at_oct 계산. 실패 시 None."""
    if not dob_str or not oct_date_str:
        return None
    try:
        d_oct = dt.strptime(oct_date_str, '%Y%m%d').date()
        d_dob = dt.strptime(dob_str,      '%Y-%m-%d').date()
        return round((d_oct - d_dob).days / 365.25, 1)
    except Exception:
        return None


def main():
    data = data_dir()
    res  = results_dir()

    vf_rows   = _read_csv(_ROOT / 'data' / 'vf_results.csv')
    oct_rows  = _read_csv(_ROOT / 'data' / 'oct_canonical.csv')
    demo_rows = _read_csv(_ROOT / 'data' / 'demographics.csv')

    # 인덱스 빌드
    oct_map  = {(r['patient_id'], r['eye']): r for r in oct_rows}
    demo_map = {r['patient_id']: r          for r in demo_rows}

    # 소스 컬럼 검증
    assert set(_VF_COLS) - {'vf_mean', 'n_absolute_scotoma'} <= set(vf_rows[0].keys()), \
        'VF 컬럼 불일치'
    assert set(_OCT_COLS) <= set(oct_rows[0].keys()), 'OCT 컬럼 불일치'

    rows_out = []
    n_oct_miss = 0
    n_demo_miss = 0

    for vf_r in vf_rows:
        pid = vf_r['patient_id']
        eye = vf_r['eye']

        row = {}

        # --- meta ---
        oct_r = oct_map.get((pid, eye), {})
        if not oct_r:
            n_oct_miss += 1
        row['patient_id'] = pid
        row['eye']        = eye
        row['vf_date']    = vf_r.get('vf_date', '')
        row['oct_date']   = oct_r.get('oct_date', '')
        row['gap_days']   = vf_r.get('gap_days', '')

        # --- demographics ---
        demo_r = demo_map.get(pid, {})
        if not demo_r:
            n_demo_miss += 1
        dob = demo_r.get('dob', '')
        sex = demo_r.get('sex', '')
        row['age_at_oct'] = _age_at(row['oct_date'], dob)
        row['sex']        = sex if sex in ('Male', 'Female') else ''

        # --- covariate_set flag ---
        has_age = row['age_at_oct'] is not None
        has_sex = row['sex'] in ('Male', 'Female')
        if has_age and has_sex:
            row['covariate_set'] = 'full'
        elif has_age:
            row['covariate_set'] = 'age'
        else:
            row['covariate_set'] = 'none'

        # --- VF ---
        for col in _VF_COLS:
            row[col] = vf_r.get(col, '')

        # --- OCT ---
        for col in _OCT_COLS:
            row[col] = oct_r.get(col, '')

        rows_out.append(row)

    # ---- 275행 검증 ----
    # 2026-07-23: 특정 사례 OD 제외(gap_days 필드가 공백이라 어떤 필터도 걸러내지
    # 못했으나, 실제 |oct_date-vf_date|=205일로 ±90일 기준 위반 확인 — 코호트
    # 무결성 전수검사에서 발견). data/vf_results.csv에서 해당 행 제거 + 나머지
    # gap_days 공백 5건은 실제 계산값으로 채움(반영: _backup_vf_results_pre_gapfix_*.csv).
    assert len(rows_out) == 275, f'행수 오류: {len(rows_out)}'
    assert n_oct_miss == 0,      f'OCT 미매칭: {n_oct_miss}'
    print(f'행수: {len(rows_out)} / OCT 미매칭: {n_oct_miss} / demo 미매칭: {n_demo_miss}')

    # ---- severity tertile (정밀 백분위수 기준, 2026-07-20 도입) ----
    q33, q67 = _assign_severity_tertiles(rows_out)
    tertile_dist = {}
    for r in rows_out:
        k = r['severity_tertile']
        tertile_dist[k] = tertile_dist.get(k, 0) + 1
    print(f'\nseverity_tertile 경계(정밀): q33={q33:.3f}dB  q67={q67:.3f}dB')
    print(f'severity_tertile 분포: mild={tertile_dist.get("mild",0)}  '
          f'moderate={tertile_dist.get("moderate",0)}  severe={tertile_dist.get("severe",0)}')

    # ---- 출력 ----
    out_path = data / 'analysis_master.csv'
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)

    # ---- 결측률 요약 ----
    def _none_rate(col):
        n = sum(1 for r in rows_out if r.get(col) in (None, '', 'None'))
        return n, len(rows_out)

    print(f'\n컬럼별 결측률:')
    for col in ['age_at_oct', 'sex',
                'vf_mean', 'avg_gcl_od', 'avg_gcl_os',
                'min_gcl_od', 'min_gcl_os',
                'rnfl_q_s', 'rnfl_q_t', 'rnfl_q_i', 'rnfl_q_n',
                'ss_gca_od', 'ss_gca_os']:
        n, total = _none_rate(col)
        print(f'  {col}: {n}/{total} ({100*n/total:.1f}%)')

    # covariate_set 분포
    cov_dist = {}
    for r in rows_out:
        k = r['covariate_set']
        cov_dist[k] = cov_dist.get(k, 0) + 1
    print(f'\ncovariate_set 분포: {cov_dist}')
    print(f'  full (age+sex 보정): {cov_dist.get("full",0)}눈')
    print(f'  age  (age만 보정):   {cov_dist.get("age",0)}눈')
    print(f'  none (주 분석만):    {cov_dist.get("none",0)}눈')
    print(f'  full+age (age 가용): {cov_dist.get("full",0)+cov_dist.get("age",0)}눈')

    # ---- SHA-256 봉인 ----
    sha = sha256_file(str(out_path))

    # 소스 해시
    src_hashes = {
        'vf_results':    sha256_file(str(_ROOT / 'data' / 'vf_results.csv')),
        'oct_canonical': sha256_file(str(_ROOT / 'data' / 'oct_canonical.csv')),
        'demographics':  sha256_file(str(_ROOT / 'data' / 'demographics.csv')),
    }

    prov = {
        'generated':  datetime.date.today().isoformat(),
        'output':     'data/analysis_master.csv',
        'sha256':     sha,
        'n_rows':     len(rows_out),
        'n_cols':     len(OUT_COLS),
        'col_groups': {
            'meta':  _META_COLS,
            'demo':  _DEMO_COLS,
            'vf':    f'vf_mean, n_absolute_scotoma, p1..p54 ({len(_VF_COLS)} cols)',
            'oct':   f'GCL avg/min/6sec x2, RNFL q3+h12, SS x4 ({len(_OCT_COLS)} cols)',
            'flag':  _FLAG_COLS,
        },
        'covariate_set': cov_dist,
        'severity_tertile': {
            'method': 'precise 33rd/67th percentile of vf_mean (numpy default linear '
                      'interpolation, unrounded); mild = MS > q67, moderate = q33 < MS <= '
                      'q67, severe = MS <= q33',
            'q33_dB': round(q33, 3),
            'q67_dB': round(q67, 3),
            'distribution': tertile_dist,
        },
        'analysis_tiers': {
            'primary':         '276 (모든 눈, covariate 없음)',
            'age_adjusted':    f'{cov_dist.get("full",0)+cov_dist.get("age",0)} (covariate_set in full/age)',
            'age_sex_adjusted':f'{cov_dist.get("full",0)} (covariate_set == full)',
        },
        'join': {
            'key_vf_oct':  '(patient_id, eye) — 1:1, OCT 미매칭 0',
            'key_demo':    'patient_id — LEFT JOIN, 결측 허용',
            'age_formula': 'age_at_oct = (oct_date - dob).days / 365.25, 소수점 1자리',
            'sex_note':    'Male/Female만; Unknown/empty → empty (Korean OCR 구조적 한계)',
        },
        'source_sha256': src_hashes,
    }
    prov_path = res / 'step5_provenance.json'
    json.dump(prov, open(prov_path, 'w', encoding='utf-8'),
              indent=2, ensure_ascii=False)

    print(f'\n출력: {out_path} ({len(OUT_COLS)}컬럼)')
    print(f'봉인: {prov_path}')
    print(f'SHA256: {sha}')


if __name__ == '__main__':
    main()
