"""
scripts/phase5_sita_sensitivity.py
SITA 전략 혼재 민감도 분석 (방어용 비교, 코호트 축소 아님)

목적:
  시야검사 전략 혼재(Standard 251 / Fast 17 / Faster 2 / 미상 6)가 전역 구조-기능
  상관에 영향을 주는지 우리 데이터로 직접 검정한다. SITA-Standard 251안만으로 β를
  재계산하여 전체 276안 β와 비교한다. β 크기가 비슷하고 CI가 겹치면 "전략 혼재는
  결과에 영향 없음"이 확립된다. **코호트는 276안 유지.**

방법:
  phase5_4_gcl_rnfl_compare.lme_single 재사용 (Z-표준화, 환자 random intercept, REML).
  지표: gcl_avg (눈별 target), rnfl_q_i, rnfl_ch_mean.
  data/sita_per_eye.csv 를 (patient_id, eye, vf_date) 로 조인.

교차검증:
  full β 가 봉인값(187 코호트 단일마커 LME)과 일치하는지 확인.
    gcl_avg=0.5723, rnfl_q_i=0.5904, rnfl_ch_mean=0.5455
  (이전 275 코호트 봉인값 0.4805/0.5476/0.5216 은 187 재생성으로 대체됨.)

출력: results/sita_sensitivity.json (비 PHI)
"""
import csv, json, math, sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'scripts'))
sys.path.insert(0, str(_ROOT / 'src'))

from phase5_4_gcl_rnfl_compare import lme_single  # noqa: E402  (동일 통계 로직 재사용)
from hvf.registry import sha256_file              # noqa: E402
from hvf.cohort import filter_rows                # noqa: E402


def _f(v):
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except (TypeError, ValueError):
        return None


# 봉인값(187 코호트 단일마커 LME) — full β 교차검증용
_SEALED_FULL_BETA = {'gcl_avg': 0.5723, 'rnfl_q_i': 0.5904, 'rnfl_ch_mean': 0.5455}
_METRICS = ['gcl_avg', 'rnfl_q_i', 'rnfl_ch_mean']


def load_joined_df():
    """analysis_master(187눈 코호트) + sita_per_eye 조인 → DataFrame."""
    master = filter_rows(list(csv.DictReader(
        open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8'))))
    sita   = list(csv.DictReader(open(_ROOT / 'data' / 'sita_per_eye.csv', encoding='utf-8-sig')))
    sita_map = {(r['patient_id'], r['eye'], r['vf_date']): r['sita_strategy'] for r in sita}

    records = []
    for r in master:
        e = r['eye'].lower()
        ch_vals  = [_f(r.get(f'rnfl_h{h:02d}')) for h in range(1, 13)]
        ch_valid = [v for v in ch_vals if v is not None]
        records.append({
            'pid':          r['patient_id'],
            'vf_mean':      _f(r.get('vf_mean')),
            'gcl_avg':      _f(r.get(f'avg_gcl_{e}')),
            'rnfl_q_i':     _f(r.get('rnfl_q_i')),
            'rnfl_ch_mean': (sum(ch_valid) / len(ch_valid)) if ch_valid else None,
            'sita_strategy': sita_map.get((r['patient_id'], r['eye'], r['vf_date']), 'UNDETERMINED'),
        })
    return pd.DataFrame(records)


def _overlap(ci_a, ci_b):
    """두 신뢰구간이 겹치는지."""
    return ci_a[0] <= ci_b[1] and ci_b[0] <= ci_a[1]


def _within(beta, ci):
    """beta 가 ci 안에 드는지."""
    return ci[0] <= beta <= ci[1]


def main():
    df = load_joined_df()
    df_std = df[df['sita_strategy'] == 'SITA Standard'].copy()

    n_total = len(df)
    n_std   = len(df_std)
    n_nonstd = (df['sita_strategy'].isin(['SITA Fast', 'SITA Faster'])).sum()
    n_undet  = (df['sita_strategy'] == 'UNDETERMINED').sum()

    print('=' * 74)
    print('SITA 전략 혼재 민감도 분석')
    print('=' * 74)
    print(f'전체 {n_total}안  |  SITA-Standard {n_std}안  |  '
          f'Fast+Faster {n_nonstd}안  |  미상 {n_undet}안')
    print()

    results = {
        'generated': '2026-07-09',
        'input_sha256': {
            'analysis_master': sha256_file(str(_ROOT / 'data' / 'analysis_master.csv')),
        },
        'design': {
            'purpose': 'SITA strategy heterogeneity sensitivity (defensive comparison, cohort NOT reduced)',
            'n_total': int(n_total),
            'n_sita_standard': int(n_std),
            'n_fast_faster': int(n_nonstd),
            'n_undetermined': int(n_undet),
            'method': 'lme_single (Z-standardized REML mixedlm, random intercept per patient), reused from phase5_4',
        },
        'metrics': {},
        'crosscheck_full_beta_vs_sealed': {},
    }

    header = (f'{"지표":14s} {"n_full":>6s} {"β_full":>8s} {"CI_full":>18s}   '
              f'{"n_std":>5s} {"β_std":>8s} {"CI_std":>18s}   {"Δβ":>7s} {"판정":>10s}')
    print(header)
    print('-' * len(header))

    all_robust = True
    for m in _METRICS:
        full = lme_single(df,     m)
        std  = lme_single(df_std, m)

        bf, cf, nf = full['lme']['beta'], full['lme']['ci'], full['n']
        bs, cs, ns = std['lme']['beta'],  std['lme']['ci'],  std['n']
        rho_f = full['spearman']['rho']
        rho_s = std['spearman']['rho']
        d_beta = round(bs - bf, 4)

        ci_overlap   = _overlap(cf, cs)
        std_in_full  = _within(bs, cf)      # standard β 가 full CI 안에 드는가 (핵심)
        full_in_std  = _within(bf, cs)

        # 판정: CI 겹치고 standard β 가 full CI 안 → 전략 효과 없음
        if ci_overlap and std_in_full:
            verdict = 'ROBUST'
        elif ci_overlap:
            verdict = 'CI겹침만'      # β는 다소 다르나 구간 겹침 (표본효과 가능)
            all_robust = False
        else:
            verdict = 'DIVERGENT'    # 전략 효과 의심 → 정직 보고
            all_robust = False

        print(f'{m:14s} {nf:6d} {bf:+8.3f} [{cf[0]:+.3f},{cf[1]:+.3f}]   '
              f'{ns:5d} {bs:+8.3f} [{cs[0]:+.3f},{cs[1]:+.3f}]   '
              f'{d_beta:+7.3f} {verdict:>10s}')

        # 봉인값 교차검증
        sealed = _SEALED_FULL_BETA.get(m)
        cc_ok = (sealed is not None and abs(bf - sealed) < 0.005)
        results['crosscheck_full_beta_vs_sealed'][m] = {
            'full_beta': bf, 'sealed_beta': sealed,
            'match': bool(cc_ok),
        }
        if sealed is not None and not cc_ok:
            print(f'   ⚠ 교차검증 실패: {m} full β={bf} vs 봉인 {sealed}')

        results['metrics'][m] = {
            'full':          {'n': nf, 'beta': bf, 'se': full['lme']['se'],
                              'p': full['lme']['p'], 'ci': cf, 'spearman_rho': rho_f},
            'standard_only': {'n': ns, 'beta': bs, 'se': std['lme']['se'],
                              'p': std['lme']['p'], 'ci': cs, 'spearman_rho': rho_s},
            'delta_beta':          d_beta,
            'ci_overlap':          bool(ci_overlap),
            'standard_beta_in_full_ci': bool(std_in_full),
            'full_beta_in_std_ci':      bool(full_in_std),
            'verdict':             verdict,
        }

    overall = 'ROBUST_no_effect' if all_robust else 'REVIEW_NEEDED'
    results['overall_verdict'] = overall

    print()
    print('교차검증 (full β vs 봉인값):')
    for m, cc in results['crosscheck_full_beta_vs_sealed'].items():
        mark = 'OK' if cc['match'] else 'MISMATCH'
        print(f'  {m:14s}: full={cc["full_beta"]:.4f}  sealed={cc["sealed_beta"]}  [{mark}]')

    print()
    print(f'종합 판정: {overall}')
    print('  (ROBUST_no_effect = 3지표 모두 CI 겹침 + standard β가 full CI 안 = 전략 혼재 영향 없음)')

    out = _ROOT / 'results' / 'sita_sensitivity.json'
    json.dump(results, open(out, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print(f'\n저장: {out}')


if __name__ == '__main__':
    main()
