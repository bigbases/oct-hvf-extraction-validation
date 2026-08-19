"""
scripts/phase5_1_correlations.py
Phase 5-1: 기본 구조-기능 상관 (GCL+IPL, RNFL vs vf_mean)

분석 설계:
  - GCL avg/min (eye-specific: OD행→avg_gcl_od, OS행→avg_gcl_os)
  - RNFL: 사분면 평균 (q_s+q_t+q_n, q_i 제외) + 개별 사분면
  - 층: OD-only / OS-only / 전체 276눈 naive / 전체 LME (환자 random intercept)
  - Pearson r + Spearman ρ, 95% CI (Fisher z), normality test
  - 재현성: seed 불필요 (순수 통계), 입력 hash 기록

출력: results/phase5_1_correlations.json
"""
import csv, json, math, datetime, sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import scipy.stats as sps

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))

from hvf.registry import sha256_file
from hvf.cohort import filter_rows


# ---------------------------------------------------------------------------
# 유틸
# ---------------------------------------------------------------------------
def _f(v):
    try:
        fv = float(v)
        return None if math.isnan(fv) else fv
    except (TypeError, ValueError):
        return None


def _fisher_ci(r, n, alpha=0.05):
    """Fisher z 변환 기반 95% CI. n >= 4 필요."""
    if n < 4 or abs(r) >= 1.0:
        return None, None
    z   = 0.5 * math.log((1 + r) / (1 - r))
    se  = 1.0 / math.sqrt(n - 3)
    zc  = 1.96
    lo  = math.tanh(z - zc * se)
    hi  = math.tanh(z + zc * se)
    return round(lo, 4), round(hi, 4)


def _corr_block(x, y, label, groups=None):
    """
    x, y: list (None 포함) — 공통 non-None 쌍만 사용
    groups: list (patient_id) — LME용, None이면 naive만
    반환: dict
    """
    pairs = [(xi, yi, gi if groups else None)
             for xi, yi, gi in zip(x, y, groups if groups else x)
             if xi is not None and yi is not None]
    if not pairs:
        return {'label': label, 'n': 0}

    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    gs = [p[2] for p in pairs]
    n  = len(xs)

    # Normality (Shapiro-Wilk, n<=5000)
    sw_x = sps.shapiro(xs) if n >= 3 else (None, None)
    sw_y = sps.shapiro(ys) if n >= 3 else (None, None)

    # Pearson
    pr, pp = sps.pearsonr(xs, ys)
    p_lo, p_hi = _fisher_ci(pr, n)

    # Spearman
    sr, sp = sps.spearmanr(xs, ys)
    s_lo, s_hi = _fisher_ci(sr, n)

    result = {
        'label':       label,
        'n':           n,
        'normality':   {
            'shapiro_x_p': round(float(sw_x[1]), 4) if sw_x[1] is not None else None,
            'shapiro_y_p': round(float(sw_y[1]), 4) if sw_y[1] is not None else None,
            'note': ('Shapiro p>0.05 → normal' if (sw_x[1] and sw_y[1] and
                     sw_x[1] > 0.05 and sw_y[1] > 0.05) else 'non-normal → Spearman primary'),
        },
        'pearson':  {'r':  round(float(pr),  4), 'p': round(float(pp),  4),
                     'ci': [p_lo, p_hi]},
        'spearman': {'rho': round(float(sr), 4), 'p': round(float(sp),  4),
                     'ci': [s_lo, s_hi]},
    }

    # LME (전체 눈 분석용 — groups 있을 때만)
    if groups is not None and n >= 10:
        try:
            import pandas as pd
            import statsmodels.formula.api as smf

            df = pd.DataFrame({'x': xs, 'y': ys, 'pid': gs})
            # Z-표준화 → standardized beta = 상관계수 등가
            df['zx'] = (df['x'] - df['x'].mean()) / df['x'].std()
            df['zy'] = (df['y'] - df['y'].mean()) / df['y'].std()

            m   = smf.mixedlm('zy ~ zx', df, groups=df['pid']).fit(reml=True)
            beta = float(m.params['zx'])
            se   = float(m.bse['zx'])
            pval = float(m.pvalues['zx'])
            lo   = float(m.conf_int().loc['zx', 0])
            hi   = float(m.conf_int().loc['zx', 1])

            result['lme'] = {
                'beta':  round(beta, 4),
                'se':    round(se,   4),
                'p':     round(pval, 4),
                'ci':    [round(lo, 4), round(hi, 4)],
                'note':  'random intercept per patient (REML)',
            }
        except Exception as e:
            result['lme'] = {'error': str(e)}

    return result


# ---------------------------------------------------------------------------
# 메인
# ---------------------------------------------------------------------------
def main():
    res_dir = _ROOT / 'results'

    # 187눈 분석 코호트로 제한 — 필터 없이 돌리면 275눈이 들어와 논문 값과 어긋난다
    # (2026-08-15 재현성 감사).
    rows = filter_rows(list(csv.DictReader(
        open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8'))))

    # eye-specific GCL 컬럼 생성
    for r in rows:
        if r['eye'] == 'OD':
            r['_gcl_avg'] = _f(r['avg_gcl_od'])
            r['_gcl_min'] = _f(r['min_gcl_od'])
        else:
            r['_gcl_avg'] = _f(r['avg_gcl_os'])
            r['_gcl_min'] = _f(r['min_gcl_os'])

        # RNFL 사분면 평균 (q_s, q_t, q_i, q_n 모두 포함) — 유효값 ≥2 필요
        rnfl_vals = [_f(r[f'rnfl_q_{q}']) for q in ('s', 't', 'i', 'n')]
        rnfl_ok   = [v for v in rnfl_vals if v is not None]
        r['_rnfl_mean'] = sum(rnfl_ok) / len(rnfl_ok) if len(rnfl_ok) >= 2 else None

        # clock-hour 평균 (h01-h12, 보조)
        ch_vals = [_f(r[f'rnfl_h{str(i).zfill(2)}']) for i in range(1, 13)]
        ch_ok   = [v for v in ch_vals if v is not None]
        r['_rnfl_ch_mean'] = sum(ch_ok) / len(ch_ok) if ch_ok else None

        r['_vf'] = _f(r['vf_mean'])
        r['_pid'] = r['patient_id']

    od  = [r for r in rows if r['eye'] == 'OD']
    os_ = [r for r in rows if r['eye'] == 'OS']
    all_rows = rows

    results = {'generated': datetime.date.today().isoformat(),
               'kcc_reference': 0.734,
               'design': {
                   'n_total': len(rows),
                   'n_od': len(od), 'n_os': len(os_),
                   'bilateral_patients': sum(1 for pid, cnt in
                       defaultdict(int, {}).items()),  # placeholder
                   'gcl_note': 'OD행→avg_gcl_od, OS행→avg_gcl_os (eye-specific)',
                   'rnfl_note': 'rnfl_quad_mean = mean(q_s, q_t, q_i, q_n) 4개 사분면, 수동검수본 q_i 복원',
                   'ci_method': 'Fisher z 변환 95% CI',
                   'lme_method': 'random intercept per patient, REML, Z-standardized beta',
               },
               'correlations': {},
               }

    # bilateral count
    from collections import Counter
    pid_cnt = Counter(r['patient_id'] for r in rows)
    results['design']['bilateral_patients'] = sum(1 for v in pid_cnt.values() if v == 2)
    results['design']['unilateral_patients'] = sum(1 for v in pid_cnt.values() if v == 1)

    def _run(pred_key, label_base, subset_rows, subset_label, pid_list=None):
        x = [r[pred_key] for r in subset_rows]
        y = [r['_vf']    for r in subset_rows]
        return _corr_block(x, y,
                           label=f'{label_base} vs vf_mean ({subset_label})',
                           groups=pid_list)

    # ---- GCL avg ----
    results['correlations']['gcl_avg'] = {
        'OD':  _run('_gcl_avg', 'avg_gcl', od,      'OD-only'),
        'OS':  _run('_gcl_avg', 'avg_gcl', os_,     'OS-only'),
        'all_naive': _run('_gcl_avg', 'avg_gcl', all_rows, 'all 276 naive'),
        'all_lme':   _run('_gcl_avg', 'avg_gcl', all_rows, 'all 276 LME',
                          pid_list=[r['_pid'] for r in all_rows]),
    }

    # ---- GCL min ----
    results['correlations']['gcl_min'] = {
        'OD':  _run('_gcl_min', 'min_gcl', od,      'OD-only'),
        'OS':  _run('_gcl_min', 'min_gcl', os_,     'OS-only'),
        'all_naive': _run('_gcl_min', 'min_gcl', all_rows, 'all 276 naive'),
        'all_lme':   _run('_gcl_min', 'min_gcl', all_rows, 'all 276 LME',
                          pid_list=[r['_pid'] for r in all_rows]),
    }

    # ---- RNFL 사분면 평균 ----
    results['correlations']['rnfl_quad_mean'] = {
        'OD':  _run('_rnfl_mean', 'rnfl_quad_mean', od,      'OD-only'),
        'OS':  _run('_rnfl_mean', 'rnfl_quad_mean', os_,     'OS-only'),
        'all_naive': _run('_rnfl_mean', 'rnfl_quad_mean', all_rows, 'all 276 naive'),
        'all_lme':   _run('_rnfl_mean', 'rnfl_quad_mean', all_rows, 'all 276 LME',
                          pid_list=[r['_pid'] for r in all_rows]),
    }

    # ---- RNFL 개별 사분면 ----
    for q in ('s', 't', 'i', 'n'):
        key = f'rnfl_q_{q}'
        qkey = f'rnfl_{q}'
        for r in all_rows:
            r[f'_{qkey}'] = _f(r[key])
        results['correlations'][qkey] = {
            'all_naive': _run(f'_{qkey}', f'rnfl_q_{q}', all_rows, 'all 276 naive'),
        }

    # ---- clock-hour 평균 ----
    results['correlations']['rnfl_ch_mean'] = {
        'all_naive': _run('_rnfl_ch_mean', 'rnfl_ch_mean', all_rows, 'all 276 naive'),
        'all_lme':   _run('_rnfl_ch_mean', 'rnfl_ch_mean', all_rows, 'all 276 LME',
                          pid_list=[r['_pid'] for r in all_rows]),
    }

    # ---- 저장 ----
    out_path = res_dir / 'phase5_1_correlations.json'
    json.dump(results, open(out_path, 'w', encoding='utf-8'),
              indent=2, ensure_ascii=False)

    # ---- 요약 출력 ----
    print('=' * 65)
    print('Phase 5-1: 기본 구조-기능 상관 결과')
    print('KCC 참조값: r = 0.734')
    print('=' * 65)

    def _fmt(blk, key='spearman'):
        if not blk or blk.get('n', 0) == 0:
            return 'n=0'
        n   = blk['n']
        d   = blk.get(key, {})
        r   = d.get('rho', d.get('r', '?'))
        p   = d.get('p', '?')
        ci  = d.get('ci', [None, None])
        p_s = f'p={p:.3f}' if isinstance(p, float) else f'p={p}'
        ci_s = (f'[{ci[0]:.3f},{ci[1]:.3f}]'
                if ci[0] is not None else '[?]')
        lme = blk.get('lme', {})
        lme_s = (f'  LME β={lme.get("beta","?"):.3f}[{lme.get("ci",["?","?"])[0]:.3f},{lme.get("ci",["?","?"])[1]:.3f}]'
                 if isinstance(lme.get('beta'), float) else '')
        return f'n={n}  ρ={r:.3f} {ci_s} {p_s}{lme_s}'

    for var, blks in results['correlations'].items():
        print(f'\n[ {var} ]')
        for sublabel, blk in blks.items():
            norm_note = ''
            if 'normality' in blk:
                if 'non-normal' in blk['normality'].get('note', ''):
                    norm_note = ' *non-normal'
            key = 'spearman'
            line = _fmt(blk, key)
            pearson_r = blk.get('pearson', {}).get('r', None)
            p_s = f'  Pearson r={pearson_r:.3f}' if isinstance(pearson_r, float) else ''
            print(f'  {sublabel:15s}: {line}{p_s}{norm_note}')

    print()
    print(f'결과 저장: {out_path}')
    # Input hash
    input_sha = sha256_file(str(_ROOT / 'data' / 'analysis_master.csv'))
    print(f'입력 SHA256: {input_sha}')

    # JSON에 input hash 추가
    results['input_sha256'] = {'analysis_master': input_sha}
    json.dump(results, open(out_path, 'w', encoding='utf-8'),
              indent=2, ensure_ascii=False)


if __name__ == '__main__':
    main()
