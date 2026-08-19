"""
scripts/phase5_3_covariate.py
Phase 5-3: 공변량 보정 3-tier LME 분석

대상: 187눈 분석 코호트(cohort_final_D.csv). load_df()가 필터를 적용한다.
아래 n은 지표별 결측 제거 전 기준이며, 실제 n은 지표마다 다르다.

Tier 1  187눈  공변량 없음           (covariate_set: all)
Tier 2         age 보정              (covariate_set: full|age)
Tier 3         age + sex 보정        (covariate_set: full)

대상 지표 (5-2 확정 대표값):
  GCL   : gcl_avg (전역 대표), gcl_min
  RNFL  : rnfl_q_i (국소 대표), rnfl_quad_mean, rnfl_ch_mean

핵심 질문: 공변량 보정 후에도 구조-기능 β가 유지되는가?
           (reviewer 방어: 구조 지표가 독립적 예측인자인가)

출력: results/phase5_3_covariate.json
"""
import csv, json, math, datetime, sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))
from hvf.registry import sha256_file
from hvf.cohort import filter_rows

# ─── OS VF flip (make_flip.py와 동일) ───────────────────────────
_VF_FLIP_ROWS = [
    list(range(1, 5)),  list(range(5, 11)),  list(range(11, 19)),
    list(range(19, 28)),list(range(28, 37)), list(range(37, 45)),
    list(range(45, 51)),list(range(51, 55)),
]

def _flip_vf(vals):
    out = list(vals)
    for grp in _VF_FLIP_ROWS:
        idxs = [i-1 for i in grp]
        sub  = [vals[i] for i in idxs]
        for idx, val in zip(idxs, reversed(sub)):
            out[idx] = val
    return out

def _f(v):
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except (TypeError, ValueError):
        return None


# ─── 데이터 로드 ─────────────────────────────────────────────────
def load_df(cohort_filter=True):
    """analysis_master.csv → 분석 DataFrame.

    cohort_filter=True(기본)면 187눈 분석 코호트(cohort_final_D.csv)로 거른다.
    tier1 은 공변량 필터가 없을 뿐 코호트 필터는 그대로 적용된다 — 이걸 빠뜨리면
    275눈 전체가 tier1 로 들어가 봉인 JSON 과 어긋난다(2026-08-15 감사).
    """
    rows = list(csv.DictReader(
        open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    if cohort_filter:
        rows = filter_rows(rows)

    records = []
    for r in rows:
        e = r['eye'].lower()

        # VF mean (OS flip, 맹점 제외 52점 평균 — master 기준 vf_mean 사용)
        vf_mean = _f(r.get('vf_mean'))

        # GCL (target eye)
        gcl_avg = _f(r.get(f'avg_gcl_{e}'))
        gcl_min = _f(r.get(f'min_gcl_{e}'))

        # RNFL 사분면
        qs  = _f(r.get('rnfl_q_s'))
        qi  = _f(r.get('rnfl_q_i'))
        qt  = _f(r.get('rnfl_q_t'))
        qn  = _f(r.get('rnfl_q_n'))
        # quad_mean: 4사분면 평균 (non-null만)
        quad_vals = [v for v in [qs, qi, qt, qn] if v is not None]
        quad_mean = sum(quad_vals)/len(quad_vals) if quad_vals else None

        # clock-hour mean (h01-h12)
        ch_vals = [_f(r.get(f'rnfl_h{h:02d}')) for h in range(1, 13)]
        ch_valid = [v for v in ch_vals if v is not None]
        ch_mean = sum(ch_valid)/len(ch_valid) if ch_valid else None

        # 공변량
        age  = _f(r.get('age_at_oct'))
        sex_raw = r.get('sex', '').strip()
        sex  = 1 if sex_raw == 'Male' else (0 if sex_raw == 'Female' else None)

        records.append({
            'pid':           r['patient_id'],
            'eye':           r['eye'],
            'covariate_set': r.get('covariate_set', 'none'),
            'vf_mean':       vf_mean,
            'gcl_avg':       gcl_avg,
            'gcl_min':       gcl_min,
            'rnfl_q_i':      qi,
            'rnfl_quad_mean':quad_mean,
            'rnfl_ch_mean':  ch_mean,
            'age':           age,
            'sex':           sex,
        })
    return pd.DataFrame(records)


# ─── LME 함수 ────────────────────────────────────────────────────
def run_lme(df_sub, struct_col, vf_col, covariates):
    """
    LME: z(vf) ~ z(struct) [+ z(age)] [+ sex]
    random intercept per patient (pid).
    Returns dict with beta, SE, p, CI95, n, n_patients, covariate_effects.
    """
    cols = [struct_col, vf_col, 'pid'] + covariates
    sub = df_sub[cols].dropna()
    n = len(sub)
    n_pat = sub['pid'].nunique()

    if n < 10:
        return {'n': n, 'note': 'too few'}

    # Z-standardize continuous vars
    sub = sub.copy()
    sub['zx'] = (sub[struct_col] - sub[struct_col].mean()) / sub[struct_col].std()
    sub['zy'] = (sub[vf_col]    - sub[vf_col].mean())    / sub[vf_col].std()

    formula_parts = ['zy ~ zx']
    if 'age' in covariates:
        sub['z_age'] = (sub['age'] - sub['age'].mean()) / sub['age'].std()
        formula_parts.append('z_age')
    if 'sex' in covariates:
        sub['sex_c'] = sub['sex']   # 0/1 binary
        formula_parts.append('sex_c')

    formula = ' + '.join(formula_parts)

    try:
        m = smf.mixedlm(formula, sub, groups=sub['pid']).fit(reml=True, disp=False)
    except Exception as ex:
        return {'n': n, 'error': str(ex)}

    def _coef(name):
        if name not in m.params.index:
            return None
        return {
            'beta': round(float(m.params[name]), 4),
            'se':   round(float(m.bse[name]),    4),
            'p':    round(float(m.pvalues[name]), 5),
            'ci':   [round(float(m.conf_int().loc[name, 0]), 4),
                     round(float(m.conf_int().loc[name, 1]), 4)],
        }

    result = {
        'n': n, 'n_patients': n_pat,
        'struct_beta': _coef('zx'),
    }

    if 'age' in covariates:
        result['age_coef'] = _coef('z_age')
    if 'sex' in covariates:
        result['sex_coef'] = _coef('sex_c')

    return result


# ─── 3-tier 분석 ──────────────────────────────────────────────────
METRICS = {
    'gcl_avg':        'gcl_avg',
    'gcl_min':        'gcl_min',
    'rnfl_q_i':       'rnfl_q_i',
    'rnfl_quad_mean': 'rnfl_quad_mean',
    'rnfl_ch_mean':   'rnfl_ch_mean',
}

TIERS = [
    ('tier1_full',  None,           []),
    ('tier2_age',   ['full','age'],  ['age']),
    ('tier3_age_sex',['full'],       ['age','sex']),
]

def main():
    df = load_df()
    print(f'로드: {len(df)}행')
    print(f'covariate_set 분포: {df["covariate_set"].value_counts().to_dict()}')

    results = {
        'generated': datetime.date.today().isoformat(),
        'input_sha256': sha256_file(str(_ROOT / 'data' / 'analysis_master.csv')),
        'design': {
            'tiers': {
                'tier1_full':   '전체 187, 공변량 없음',
                'tier2_age':    'covariate_set in (full,age), age 보정',
                'tier3_age_sex':'covariate_set==full, age+sex 보정',
            },
            'outcome': 'vf_mean (Z-std)',
            'lme': 'random intercept per patient (REML)',
            'covariates_z': 'age Z-std; sex binary 0/1 (Female=0, Male=1)',
        },
        'by_metric': {},
        'summary_table': [],
    }

    # ── 3-tier × 5-metric 계산 ──
    tbl = []   # for summary table
    for metric, col in METRICS.items():
        row_res = {}
        for tier_name, cov_filter, covs in TIERS:
            if cov_filter is None:
                sub = df.copy()
            else:
                sub = df[df['covariate_set'].isin(cov_filter)].copy()

            res = run_lme(sub, col, 'vf_mean', covs)
            row_res[tier_name] = res

            # console 출력
            if 'struct_beta' in res and res['struct_beta']:
                sb = res['struct_beta']
                cov_str = ''
                if 'age_coef' in res and res['age_coef']:
                    ac = res['age_coef']
                    cov_str += f'  age β={ac["beta"]:+.3f}(p={ac["p"]:.3f})'
                if 'sex_coef' in res and res['sex_coef']:
                    sc = res['sex_coef']
                    cov_str += f'  sex β={sc["beta"]:+.3f}(p={sc["p"]:.3f})'
                print(f'  {metric:18s} | {tier_name:18s}: β={sb["beta"]:+.3f} SE={sb["se"]:.3f} p={sb["p"]:.4f} n={res["n"]}{cov_str}')
            else:
                print(f'  {metric:18s} | {tier_name:18s}: ERR {res.get("error","?")} n={res.get("n","?")}')

        results['by_metric'][metric] = row_res

        # summary table row (tier1 → tier2 → tier3 beta 변화)
        def _b(t): return row_res[t].get('struct_beta', {}).get('beta') if row_res[t].get('struct_beta') else None
        def _p(t): return row_res[t].get('struct_beta', {}).get('p')    if row_res[t].get('struct_beta') else None
        def _n(t): return row_res[t].get('n')

        b1, b2, b3 = _b('tier1_full'), _b('tier2_age'), _b('tier3_age_sex')
        p1, p2, p3 = _p('tier1_full'), _p('tier2_age'), _p('tier3_age_sex')
        tbl.append({
            'metric': metric,
            'tier1_beta': round(b1, 3) if b1 is not None else None,
            'tier1_p':    round(p1, 5) if p1 is not None else None,
            'tier1_n':    _n('tier1_full'),
            'tier2_beta': round(b2, 3) if b2 is not None else None,
            'tier2_p':    round(p2, 5) if p2 is not None else None,
            'tier2_n':    _n('tier2_age'),
            'tier3_beta': round(b3, 3) if b3 is not None else None,
            'tier3_p':    round(p3, 5) if p3 is not None else None,
            'tier3_n':    _n('tier3_age_sex'),
            'beta_stable': (
                b1 is not None and b3 is not None and
                abs(b3 - b1) / abs(b1) < 0.20  # 20% 미만 변화 = 안정
            ),
        })

    results['summary_table'] = tbl

    # ── 공변량 효과 요약 ──
    print('\n── 공변량 효과 요약 ──────────────────────────────────')
    for metric, col in METRICS.items():
        t2 = results['by_metric'][metric]['tier2_age']
        t3 = results['by_metric'][metric]['tier3_age_sex']
        ac2 = t2.get('age_coef')
        ac3 = t3.get('age_coef')
        sc3 = t3.get('sex_coef')
        if ac2:
            print(f'  {metric:18s} age(T2): β={ac2["beta"]:+.4f} p={ac2["p"]:.4f}')
        if sc3:
            print(f'  {metric:18s} sex(T3): β={sc3["beta"]:+.4f} p={sc3["p"]:.4f}  (Male=1 vs Female=0)')

    # ── 결과 요약 테이블 출력 ──
    print('\n── 3-tier 요약 테이블 ───────────────────────────────────')
    print(f'  {"지표":18s}  T1(276) β/p         T2(~260) β/p        T3(~181) β/p        stable?')
    print('  ' + '-'*90)
    for row in tbl:
        def fmt(b, p):
            if b is None or p is None: return '  N/A              '
            p_disp = '<.0001' if p < 0.0001 else f'{p:.4f}'
            sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
            return f'β={b:+.3f} p={p_disp} {sig}'
        s = 'STABLE' if row['beta_stable'] else 'SHIFT'
        print(f'  {row["metric"]:18s}  {fmt(row["tier1_beta"],row["tier1_p"])}  '
              f'{fmt(row["tier2_beta"],row["tier2_p"])}  '
              f'{fmt(row["tier3_beta"],row["tier3_p"])}  {s}  '
              f'n={row["tier1_n"]}/{row["tier2_n"]}/{row["tier3_n"]}')

    # 저장
    out_path = _ROOT / 'results' / 'phase5_3_covariate.json'
    json.dump(results, open(out_path, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print(f'\n출력: {out_path}')


if __name__ == '__main__':
    main()
