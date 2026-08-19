"""
scripts/phase5_4_gcl_rnfl_compare.py
Phase 5-4: GCL vs RNFL 직접 비교

대상: 187눈 분석 코호트(cohort_final_D.csv). load_df()가 필터를 적용한다.

1. 전역 head-to-head: gcl_avg vs rnfl_q_i / rnfl_quad_mean / rnfl_ch_mean
   - LME β, 95% CI, Spearman ρ (같은 187눈, random intercept per patient)
   - Hotelling T² 유사 검정: paired difference test (각 눈에서 두 struct β 비교)

2. 중증도 층화: mild(상위 1/3) vs moderate+severe(하위 2/3)
   - gcl_avg vs rnfl_q_i 각 층에서
   - "언제 어느 게 나은가" = 조기 진단 함의

3. 상보성: GCL + RNFL 동시 투입 LME
   - zy ~ z_gcl + z_rnfl + (1|pid)
   - 두 β 모두 유의하면 → 상보적 (독립 기여)
   - 한쪽 β 소멸하면 → 중복

출력: results/phase5_4_gcl_rnfl_compare.json
"""
import csv, json, math, datetime, sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.stats as sps
import statsmodels.formula.api as smf

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))
from hvf.registry import sha256_file
from hvf.cohort import filter_rows

# ─── VF flip (make_flip.py 동일) ─────────────────────────────────
_VF_FLIP_ROWS = [
    list(range(1, 5)), list(range(5, 11)), list(range(11, 19)),
    list(range(19, 28)), list(range(28, 37)), list(range(37, 45)),
    list(range(45, 51)), list(range(51, 55)),
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
    논문 Table 6/7 은 이 코호트 기준이며, 필터 없이(275눈) 돌리면 계수가 달라진다
    — 2026-08-15 감사에서 이 필터 누락으로 봉인 JSON 재현이 실패했었다.
    """
    rows = list(csv.DictReader(
        open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    if cohort_filter:
        rows = filter_rows(rows)
    records = []
    for r in rows:
        e = r['eye'].lower()
        vf_mean = _f(r.get('vf_mean'))
        gcl_avg = _f(r.get(f'avg_gcl_{e}'))
        gcl_min = _f(r.get(f'min_gcl_{e}'))
        qs  = _f(r.get('rnfl_q_s'))
        qi  = _f(r.get('rnfl_q_i'))
        qt  = _f(r.get('rnfl_q_t'))
        qn  = _f(r.get('rnfl_q_n'))
        qv  = [v for v in [qs, qi, qt, qn] if v is not None]
        quad_mean = sum(qv)/len(qv) if qv else None
        ch_vals   = [_f(r.get(f'rnfl_h{h:02d}')) for h in range(1, 13)]
        ch_valid  = [v for v in ch_vals if v is not None]
        ch_mean   = sum(ch_valid)/len(ch_valid) if ch_valid else None
        records.append({
            'pid':      r['patient_id'],
            'eye':      r['eye'],
            'vf_mean':  vf_mean,
            'gcl_avg':  gcl_avg,
            'gcl_min':  gcl_min,
            'rnfl_q_i': qi,
            'rnfl_quad_mean': quad_mean,
            'rnfl_ch_mean':   ch_mean,
        })
    return pd.DataFrame(records)


# ─── LME 단일 예측변수 ────────────────────────────────────────────
def lme_single(df_sub, x_col, y_col='vf_mean'):
    sub = df_sub[[x_col, y_col, 'pid']].dropna()
    n, n_pat = len(sub), sub['pid'].nunique()
    if n < 10:
        return {'n': n, 'note': 'too few'}
    sub = sub.copy()
    sub['zx'] = (sub[x_col]  - sub[x_col].mean())  / sub[x_col].std()
    sub['zy'] = (sub[y_col]  - sub[y_col].mean())  / sub[y_col].std()
    # Spearman
    sr, sp = sps.spearmanr(sub[x_col].values, sub[y_col].values)
    pr, pp = sps.pearsonr(sub[x_col].values, sub[y_col].values)
    def _fci(r, n):
        if abs(r) >= 1: return [None, None]
        z  = 0.5 * math.log((1+r)/(1-r))
        se = 1 / math.sqrt(n - 3)
        return [round(math.tanh(z - 1.96*se), 4), round(math.tanh(z + 1.96*se), 4)]
    try:
        m = smf.mixedlm('zy ~ zx', sub, groups=sub['pid']).fit(reml=True, disp=False)
        beta = float(m.params['zx'])
        se_b = float(m.bse['zx'])
        ci   = [float(m.conf_int().loc['zx', 0]), float(m.conf_int().loc['zx', 1])]
        pval = float(m.pvalues['zx'])
    except Exception as ex:
        return {'n': n, 'error': str(ex)}
    return {
        'n': n, 'n_patients': n_pat,
        'spearman': {'rho': round(float(sr), 4), 'p': round(float(sp), 5), 'ci': _fci(float(sr), n)},
        'pearson':  {'r':   round(float(pr), 4), 'p': round(float(pp), 5)},
        'lme':      {'beta': round(beta, 4), 'se': round(se_b, 4),
                     'p': round(pval, 5), 'ci': [round(c, 4) for c in ci]},
    }


# ─── LME 이중 예측변수 (상보성) ──────────────────────────────────
def lme_joint(df_sub, x1_col, x2_col, y_col='vf_mean'):
    sub = df_sub[[x1_col, x2_col, y_col, 'pid']].dropna()
    n, n_pat = len(sub), sub['pid'].nunique()
    if n < 15:
        return {'n': n, 'note': 'too few'}
    sub = sub.copy()
    for col in [x1_col, x2_col, y_col]:
        sub[f'z_{col}'] = (sub[col] - sub[col].mean()) / sub[col].std()
    try:
        formula = f'z_{y_col} ~ z_{x1_col} + z_{x2_col}'
        m = smf.mixedlm(formula, sub, groups=sub['pid']).fit(reml=True, disp=False)
        def _coef(name):
            if name not in m.params.index: return None
            ci = [float(m.conf_int().loc[name, 0]), float(m.conf_int().loc[name, 1])]
            return {'beta': round(float(m.params[name]), 4),
                    'se':   round(float(m.bse[name]),    4),
                    'p':    round(float(m.pvalues[name]),5),
                    'ci':   [round(c, 4) for c in ci]}
        return {
            'n': n, 'n_patients': n_pat,
            x1_col: _coef(f'z_{x1_col}'),
            x2_col: _coef(f'z_{x2_col}'),
        }
    except Exception as ex:
        return {'n': n, 'error': str(ex)}


# ─── 상보성 진단: VIF + AIC (ML 기반) ────────────────────────────
def complementarity_stats(df_sub, x1_col, x2_col, y_col='vf_mean'):
    """
    GCL+RNFL 결합 모델의 공선성(VIF)과 모델 적합도(AIC) 진단.
    - VIF: 2예측변수이므로 VIF = 1/(1-r²), r=두 예측변수 상관 (양쪽 동일)
    - AIC: 단독(GCL만, RNFL만) vs 결합 모델. 고정효과가 다르므로 ML(reml=False).
           세 모델 모두 동일 dropna 표본에서 적합해야 비교 가능.
    - ΔAIC = AIC_단독 - AIC_결합  (양수 = 결합이 더 우수)
    """
    sub = df_sub[[x1_col, x2_col, y_col, 'pid']].dropna()
    n, n_pat = len(sub), sub['pid'].nunique()
    if n < 15:
        return {'n': n, 'note': 'too few'}
    sub = sub.copy()
    for col in [x1_col, x2_col, y_col]:
        sub[f'z_{col}'] = (sub[col] - sub[col].mean()) / sub[col].std()

    # VIF: 두 예측변수 간 상관 기반
    r_xx = float(np.corrcoef(sub[f'z_{x1_col}'], sub[f'z_{x2_col}'])[0, 1])
    vif = 1.0 / (1.0 - r_xx**2) if abs(r_xx) < 1 else float('inf')

    try:
        # ML 적합 (reml=False) — AIC 비교용, 동일 표본
        m_gcl = smf.mixedlm(f'z_{y_col} ~ z_{x1_col}', sub,
                            groups=sub['pid']).fit(reml=False, disp=False)
        m_rnfl = smf.mixedlm(f'z_{y_col} ~ z_{x2_col}', sub,
                             groups=sub['pid']).fit(reml=False, disp=False)
        m_joint = smf.mixedlm(f'z_{y_col} ~ z_{x1_col} + z_{x2_col}', sub,
                              groups=sub['pid']).fit(reml=False, disp=False)
        aic_gcl   = float(m_gcl.aic)
        aic_rnfl  = float(m_rnfl.aic)
        aic_joint = float(m_joint.aic)
    except Exception as ex:
        return {'n': n, 'vif': round(vif, 4), 'aic_error': str(ex)}

    return {
        'n': n, 'n_patients': n_pat,
        'predictor_corr_r': round(r_xx, 4),
        'vif': round(vif, 4),
        'aic_gcl_only':  round(aic_gcl, 2),
        'aic_rnfl_only': round(aic_rnfl, 2),
        'aic_joint':     round(aic_joint, 2),
        'delta_aic_vs_gcl_only':  round(aic_gcl  - aic_joint, 2),
        'delta_aic_vs_rnfl_only': round(aic_rnfl - aic_joint, 2),
        'aic_note': 'ML(reml=False), same dropna sample; positive ΔAIC = joint model better',
    }


# ─── 중증도 3분위 기준 ────────────────────────────────────────────
def get_severity_groups(df):
    vm = df['vf_mean'].dropna().sort_values().values
    q33 = vm[int(len(vm)*0.333)]
    q67 = vm[int(len(vm)*0.667)]
    mild   = df[df['vf_mean'] > q67].copy()
    modsev = df[df['vf_mean'] <= q67].copy()
    return mild, modsev, q33, q67


def main():
    df = load_df()
    mild, modsev, q33, q67 = get_severity_groups(df)
    print(f'전체 n={len(df)}  mild n={len(mild)}(vf>{q67:.1f}dB)  mod+sev n={len(modsev)}')

    results = {
        'generated': datetime.date.today().isoformat(),
        'input_sha256': sha256_file(str(_ROOT / 'data' / 'analysis_master.csv')),
        'severity_cutoffs': {'mild_threshold_dB': round(float(q67), 1),
                             'lower_tertile_dB': round(float(q33), 1)},
        'head_to_head': {},
        'severity_stratified': {},
        'complementarity': {},
    }

    # ══════════════════════════════════════════════════════════════
    # 1. 전역 head-to-head
    # ══════════════════════════════════════════════════════════════
    print('\n' + '='*70)
    print(f'1. 전역 Head-to-Head (n={len(df)})')
    print('='*70)
    struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']
    print(f'\n  {"지표":18s}  {"ρ":6s}  {"LME β":7s}  {"95% CI β":18s}  {"p":8s}  n')
    print('  '+'-'*72)
    h2h_res = {}
    for col in struct_cols:
        res = lme_single(df, col)
        h2h_res[col] = res
        if 'lme' in res:
            rho = res['spearman']['rho']
            b   = res['lme']['beta']
            ci  = res['lme']['ci']
            p   = res['lme']['p']
            n   = res['n']
            sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
            sep = ' ← GCL' if col.startswith('gcl') else ' ← RNFL'
            print(f'  {col:18s}  {rho:6.3f}  {b:+7.3f}  [{ci[0]:+.3f}, {ci[1]:+.3f}]  {p:.5f} {sig:3s}{sep}')
    results['head_to_head'] = h2h_res

    # CI 겹침 확인 (gcl_avg vs rnfl_q_i)
    ci_gcl  = h2h_res['gcl_avg']['lme']['ci']
    ci_rnfl = h2h_res['rnfl_q_i']['lme']['ci']
    overlap = ci_gcl[0] < ci_rnfl[1] and ci_rnfl[0] < ci_gcl[1]
    b_gcl  = h2h_res['gcl_avg']['lme']['beta']
    b_rnfl = h2h_res['rnfl_q_i']['lme']['beta']
    print(f'\n  gcl_avg β={b_gcl:+.3f} [{ci_gcl[0]:+.3f},{ci_gcl[1]:+.3f}]')
    print(f'  rnfl_q_i β={b_rnfl:+.3f} [{ci_rnfl[0]:+.3f},{ci_rnfl[1]:+.3f}]')
    print(f'  CI {"겹침" if overlap else "비겹침"} → {"유의한 차이 없음" if overlap else "유의한 차이 있음"}')
    results['head_to_head']['gcl_vs_rnfl_ci_overlap'] = bool(overlap)

    # ══════════════════════════════════════════════════════════════
    # 2. 중증도 층화
    # ══════════════════════════════════════════════════════════════
    print('\n' + '='*70)
    print(f'2. 중증도 층화  mild(vf>{q67:.1f}dB, n={len(mild)}) vs mod+sev(n={len(modsev)})')
    print('='*70)

    for group_name, group_df in [('mild', mild), ('mod+severe', modsev)]:
        print(f'\n  ── {group_name} (n={len(group_df)}) ──')
        print(f'  {"지표":18s}  {"ρ":6s}  {"LME β":7s}  {"95% CI β":18s}  p')
        print('  '+'-'*60)
        strat_res = {}
        for col in ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']:
            res = lme_single(group_df, col)
            strat_res[col] = res
            if 'lme' in res:
                rho = res['spearman']['rho']
                b   = res['lme']['beta']
                ci  = res['lme']['ci']
                p   = res['lme']['p']
                sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else '(ns)'))
                print(f'  {col:18s}  {rho:6.3f}  {b:+7.3f}  [{ci[0]:+.3f}, {ci[1]:+.3f}]  {sig}')
        results['severity_stratified'][group_name] = strat_res

    # ══════════════════════════════════════════════════════════════
    # 3. 상보성: GCL + RNFL 동시 투입
    # ══════════════════════════════════════════════════════════════
    print('\n' + '='*70)
    print('3. 상보성: GCL + RNFL 동시 투입 LME')
    print('='*70)

    combos = [
        ('gcl_avg',    'rnfl_q_i',       'primary combo'),
        ('gcl_avg',    'rnfl_quad_mean',  'avg vs quad_mean'),
        ('gcl_min',    'rnfl_q_i',       'min vs q_i'),
        ('gcl_avg',    'rnfl_ch_mean',   'avg vs ch_mean'),
    ]
    print(f'\n  {"조합":30s}  {"GCL β":7s} {"p":7s}  {"RNFL β":7s} {"p":7s}  n  판정')
    print('  '+'-'*80)
    comp_res = {}
    for x1, x2, label in combos:
        res = lme_joint(df, x1, x2)
        # VIF/AIC 진단 추가
        diag = complementarity_stats(df, x1, x2)
        for k in ('vif', 'predictor_corr_r', 'aic_gcl_only', 'aic_rnfl_only',
                  'aic_joint', 'delta_aic_vs_gcl_only', 'delta_aic_vs_rnfl_only',
                  'aic_note'):
            if k in diag:
                res[k] = diag[k]
        comp_res[f'{x1}__{x2}'] = res
        if x1 in res and x2 in res and res[x1] and res[x2]:
            b1, p1 = res[x1]['beta'], res[x1]['p']
            b2, p2 = res[x2]['beta'], res[x2]['p']
            n = res['n']
            s1 = '***' if p1<0.001 else ('**' if p1<0.01 else ('*' if p1<0.05 else 'ns'))
            s2 = '***' if p2<0.001 else ('**' if p2<0.01 else ('*' if p2<0.05 else 'ns'))
            if s1 != 'ns' and s2 != 'ns':
                verdict = 'COMPLEMENTARY'
            elif s1 == 'ns' and s2 != 'ns':
                verdict = 'RNFL dominates'
            elif s1 != 'ns' and s2 == 'ns':
                verdict = 'GCL dominates'
            else:
                verdict = 'neither sig'
            vif = res.get('vif', '?')
            da_g = res.get('delta_aic_vs_gcl_only', '?')
            da_r = res.get('delta_aic_vs_rnfl_only', '?')
            print(f'  {x1:12s}+{x2:15s}  {b1:+7.3f} {s1:3s}  {b2:+7.3f} {s2:3s}  {n}  {verdict}')
            print(f'  {"":30s}  VIF={vif}  ΔAIC(vs GCL단독)={da_g}  ΔAIC(vs RNFL단독)={da_r}')
    results['complementarity'] = comp_res

    # ══════════════════════════════════════════════════════════════
    # 단독 β 비교 vs 동시 투입 β (변화 확인)
    # ══════════════════════════════════════════════════════════════
    print('\n── 단독 vs 동시 투입: β 감소량 (gcl_avg + rnfl_q_i) ──────────────')
    solo_gcl  = h2h_res['gcl_avg']['lme']['beta']
    solo_rnfl = h2h_res['rnfl_q_i']['lme']['beta']
    key = 'gcl_avg__rnfl_q_i'
    joint = results['complementarity'].get(key, {})
    if 'gcl_avg' in joint and 'rnfl_q_i' in joint and joint['gcl_avg'] and joint['rnfl_q_i']:
        j_gcl  = joint['gcl_avg']['beta']
        j_rnfl = joint['rnfl_q_i']['beta']
        d_gcl  = solo_gcl  - j_gcl
        d_rnfl = solo_rnfl - j_rnfl
        print(f'  gcl_avg:  단독 β={solo_gcl:+.3f} → 동시 β={j_gcl:+.3f}  감소 {d_gcl:+.3f}({d_gcl/solo_gcl*100:+.0f}%)')
        print(f'  rnfl_q_i: 단독 β={solo_rnfl:+.3f} → 동시 β={j_rnfl:+.3f}  감소 {d_rnfl:+.3f}({d_rnfl/solo_rnfl*100:+.0f}%)')
        print(f'  → 각 β 감소량이 작으면 중복 적음 = 상보적')
        results['complementarity']['beta_attenuation'] = {
            'gcl_avg_solo':  solo_gcl,  'gcl_avg_joint':  j_gcl,
            'rnfl_q_i_solo': solo_rnfl, 'rnfl_q_i_joint': j_rnfl,
            'gcl_attenuation_pct':  round(d_gcl/solo_gcl*100, 1),
            'rnfl_attenuation_pct': round(d_rnfl/solo_rnfl*100, 1),
        }

    # 저장
    out_path = _ROOT / 'results' / 'phase5_4_gcl_rnfl_compare.json'
    json.dump(results, open(out_path, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print(f'\n출력: {out_path}')


if __name__ == '__main__':
    main()
