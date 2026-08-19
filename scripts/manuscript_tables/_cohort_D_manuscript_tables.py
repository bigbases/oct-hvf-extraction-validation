# -*- coding: utf-8 -*-
"""
코호트 D(n=187) 확정 기준 원고용 최종 수치 일괄 생성.
tab:sf/tab:comp/tab:cohort/276-vs-187 병치/민감도분석(GT vs 자동) 전부 187 기준.
275/340 base·analysis_master.csv·phase5_* 등 기존 파이프라인 파일은 읽기만 함.
출력: _cohort_D_manuscript_tables.txt (한 파일에 전부)
"""
import csv, sys, io
from pathlib import Path
from collections import Counter

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
from phase5_4_gcl_rnfl_compare import load_df, lme_single, lme_joint, complementarity_stats  # noqa: E402
import pandas as pd
import numpy as np


def to_f(v):
    try: return float(v)
    except (TypeError, ValueError): return None


def build_eye_metrics():
    rel = {(r['patient_id'], r['eye'], r['vf_date']): r
           for r in csv.DictReader(open(_ROOT / 'vf_reliability.csv', encoding='utf-8-sig'))}
    am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    out = {}
    for r in am:
        pid, eye, vfdate = r['patient_id'], r['eye'], r['vf_date']
        ss_gca  = to_f(r['ss_gca_od'])  if eye == 'OD' else to_f(r['ss_gca_os'])
        ss_rnfl = to_f(r['ss_rnfl_od']) if eye == 'OD' else to_f(r['ss_rnfl_os'])
        bad_ss = (ss_gca is not None and ss_gca < 6) or (ss_rnfl is not None and ss_rnfl < 6)
        key = (pid, eye, vfdate)
        rv = rel.get(key)
        if rv is None:
            out[(pid, eye)] = dict(fp=None, fn=None, bad_ss=bad_ss)
            continue
        fp = to_f(rv['false_pos']); fn = to_f(rv['false_neg'])
        out[(pid, eye)] = dict(fp=fp, fn=fn, bad_ss=bad_ss)
    return out


def main():
    out = io.StringIO()
    def P(*a):
        print(*a)
        print(*a, file=out)

    m = build_eye_metrics()
    bad_D = {}
    for k, d in m.items():
        b_fp = d['fp'] is not None and d['fp'] > 15
        b_fn = d['fn'] is not None and d['fn'] > 20
        bad_D[k] = b_fp or b_fn or d['bad_ss']

    am_rows = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    n_A = len(am_rows)
    n_D = sum(1 for k in bad_D if not bad_D[k])

    df = load_df()
    df['bad_D'] = df.apply(lambda r: bad_D.get((r['pid'], r['eye']), False), axis=1)
    df_A, df_D = df, df[~df['bad_D']].copy()

    P('=' * 90)
    P(f'코호트 D 확정: n={n_D}  (A=276원본/{n_A}실사용 base 대비)')
    P('=' * 90)

    # ── tab:sf ──
    P(f'\n[tab:sf] 5개 지표 Spearman rho / LME beta / 95%CI -- D(n={n_D}) 기준')
    struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']
    tab_sf = {}
    P(f'  {"지표":16s}  {"n":5s}  {"Spearman rho":14s}  {"LME beta":10s}  {"95% CI":20s}  p')
    for col in struct_cols:
        r = lme_single(df_D, col)
        tab_sf[col] = r
        rho = r['spearman']['rho']; b = r['lme']['beta']; ci = r['lme']['ci']; p = r['lme']['p']
        P(f'  {col:16s}  {r["n"]:5d}  {rho:14.4f}  {b:+10.4f}  [{ci[0]:+.4f},{ci[1]:+.4f}]  {p:.5f}')

    # ── tab:comp ──
    P(f'\n[tab:comp] joint model (gcl_avg + rnfl_q_i) -- D(n={n_D})')
    j_D = lme_joint(df_D, 'gcl_avg', 'rnfl_q_i')
    diag_D = complementarity_stats(df_D, 'gcl_avg', 'rnfl_q_i')
    g = j_D['gcl_avg']; rr = j_D['rnfl_q_i']
    P(f'  n={j_D["n"]}  gcl_avg beta={g["beta"]:+.4f}(p={g["p"]:.5f})  '
      f'rnfl_q_i beta={rr["beta"]:+.4f}(p={rr["p"]:.5f})')
    P(f'  VIF={diag_D.get("vif")}  '
      f'deltaAIC(vs GCL-only)={diag_D.get("delta_aic_vs_gcl_only")}  '
      f'deltaAIC(vs RNFL-only)={diag_D.get("delta_aic_vs_rnfl_only")}')
    P(f'\n  참고: 사용자 언급 "176 기준 deltaAIC=21.09" 대비 '
      f'D(187) deltaAIC(vsGCL)={diag_D.get("delta_aic_vs_gcl_only")} '
      f'-> {"개선" if diag_D.get("delta_aic_vs_gcl_only",0) > 21.09 else "약화"}')

    # ── tab:cohort (D=187) ──
    P(f'\n[tab:cohort] 코호트 D(n={n_D}) 특성')
    d_pids = {k[0] for k in bad_D if not bad_D[k]}
    d_rows = [r for r in am_rows if (r['patient_id'], r['eye']) in {k for k in bad_D if not bad_D[k]}]
    ages = [to_f(r['age_at_oct']) for r in d_rows if to_f(r['age_at_oct']) is not None]
    sexes = Counter(r['sex'] for r in d_rows if r['sex'])
    ms = [to_f(r['vf_mean']) for r in d_rows if to_f(r['vf_mean']) is not None]
    sev = Counter(r['severity_tertile'] for r in d_rows if r['severity_tertile'])
    n_patients_D = len({r['patient_id'] for r in d_rows})

    sita = {(r['patient_id'], r['eye'], r['vf_date']): r['sita_strategy']
            for r in csv.DictReader(open(_ROOT / 'data' / 'sita_per_eye.csv', encoding='utf-8-sig'))}
    strat_counter = Counter()
    for r in d_rows:
        s = sita.get((r['patient_id'], r['eye'], r['vf_date']), 'unknown')
        strat_counter[s] += 1

    P(f'  n(안)={len(d_rows)}  n(환자)={n_patients_D}')
    P(f'  나이: mean={np.mean(ages):.1f} sd={np.std(ages,ddof=1):.1f} (n={len(ages)}, 결측 {len(d_rows)-len(ages)})')
    P(f'  성별: {dict(sexes)}')
    P(f'  MS(vf_mean): mean={np.mean(ms):.2f} sd={np.std(ms,ddof=1):.2f} (n={len(ms)})')
    P(f'  severity 분포: {dict(sev)}')
    P(f'  VF strategy 분포: {dict(strat_counter)}')

    # ── 276 vs 187 (138 참조) ──
    P('\n[276 vs 187 병치] (138 참조용)')
    bad_B = {}
    rel = {(r['patient_id'], r['eye'], r['vf_date']): r
           for r in csv.DictReader(open(_ROOT / 'vf_reliability.csv', encoding='utf-8-sig'))}
    for r in am_rows:
        pid, eye, vfdate = r['patient_id'], r['eye'], r['vf_date']
        ss_gca  = to_f(r['ss_gca_od'])  if eye == 'OD' else to_f(r['ss_gca_os'])
        ss_rnfl = to_f(r['ss_rnfl_od']) if eye == 'OD' else to_f(r['ss_rnfl_os'])
        bad_ss = (ss_gca is not None and ss_gca < 6) or (ss_rnfl is not None and ss_rnfl < 6)
        key = (pid, eye, vfdate)
        rv = rel.get(key)
        if rv is None:
            bad_B[(pid, eye)] = bad_ss
            continue
        num, den = to_f(rv['fix_loss_num']), to_f(rv['fix_loss_denom'])
        fl = 100*num/den if (num is not None and den) else None
        fp = to_f(rv['false_pos']); fn = to_f(rv['false_neg']); low = rv['low_reliability']
        b_fl = fl is not None and fl > 33
        b_fp = fp is not None and fp > 15
        b_fn = fn is not None and fn > 20
        b_low = low == 'Y'
        bad_B[(pid, eye)] = b_fl or b_fp or b_fn or b_low or bad_ss
    n_B = sum(1 for v in bad_B.values() if not v)
    df['bad_B'] = df.apply(lambda r: bad_B.get((r['pid'], r['eye']), False), axis=1)
    df_B = df[~df['bad_B']].copy()

    P(f'  n: A={n_A}  B(참조)={n_B}  D={n_D}')
    for col in struct_cols:
        rA = lme_single(df_A, col); rB = lme_single(df_B, col); rD = lme_single(df_D, col)
        P(f'  {col:16s}  A_beta={rA["lme"]["beta"]:+.4f}  B_beta={rB["lme"]["beta"]:+.4f}  D_beta={rD["lme"]["beta"]:+.4f}')

    out.seek(0)
    with open(_ROOT / '_cohort_D_manuscript_tables.txt', 'w', encoding='utf-8') as f:
        f.write(out.getvalue())


if __name__ == '__main__':
    main()
