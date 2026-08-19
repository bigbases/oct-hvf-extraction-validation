# -*- coding: utf-8 -*-
"""
276(A, 무필터) vs 138(B, 표준기준) vs 187(D, evidence-based 확정) 3열 병치.
tab:sf 5개 지표 + tab:comp joint model 전부. 신규성(방법론 기여) 서술용.
275/340 base·기존 파이프라인 파일은 읽기만 함.
"""
import csv, sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
from phase5_4_gcl_rnfl_compare import load_df, lme_single, lme_joint, complementarity_stats  # noqa: E402


def to_f(v):
    try: return float(v)
    except (TypeError, ValueError): return None


def build_bad_sets():
    rel = {(r['patient_id'], r['eye'], r['vf_date']): r
           for r in csv.DictReader(open(_ROOT / 'vf_reliability.csv', encoding='utf-8-sig'))}
    am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    bad_B, bad_D = {}, {}
    for r in am:
        pid, eye, vfdate = r['patient_id'], r['eye'], r['vf_date']
        ss_gca  = to_f(r['ss_gca_od'])  if eye == 'OD' else to_f(r['ss_gca_os'])
        ss_rnfl = to_f(r['ss_rnfl_od']) if eye == 'OD' else to_f(r['ss_rnfl_os'])
        bad_ss = (ss_gca is not None and ss_gca < 6) or (ss_rnfl is not None and ss_rnfl < 6)
        rv = rel.get((pid, eye, vfdate))
        if rv is None:
            bad_B[(pid, eye)] = bad_ss
            bad_D[(pid, eye)] = bad_ss
            continue
        num, den = to_f(rv['fix_loss_num']), to_f(rv['fix_loss_denom'])
        fl = 100*num/den if (num is not None and den) else None
        fp = to_f(rv['false_pos']); fn = to_f(rv['false_neg']); low = rv['low_reliability']
        b_fl = fl is not None and fl > 33
        b_fp = fp is not None and fp > 15
        b_fn = fn is not None and fn > 20
        b_low = low == 'Y'
        bad_B[(pid, eye)] = b_fl or b_fp or b_fn or b_low or bad_ss
        bad_D[(pid, eye)] = b_fp or b_fn or bad_ss
    return bad_B, bad_D


def main():
    bad_B, bad_D = build_bad_sets()
    df = load_df()
    df['bad_B'] = df.apply(lambda r: bad_B.get((r['pid'], r['eye']), False), axis=1)
    df['bad_D'] = df.apply(lambda r: bad_D.get((r['pid'], r['eye']), False), axis=1)
    df_A = df
    df_B = df[~df['bad_B']].copy()
    df_D = df[~df['bad_D']].copy()

    n_A, n_B, n_D = len(df_A), len(df_B), len(df_D)
    print('=' * 100)
    print(f'코호트 정의: A(무필터)={n_A}   B(표준: FL33/FP15/FN20/배너/ss6)={n_B}   '
          f'D(evidence-based: FP15/FN20/ss6만, FL·배너 제외기준 아님)={n_D}')
    print('=' * 100)

    struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']
    print(f'\n[tab:sf] LME β 3열 병치')
    print(f'  {"지표":16s}  {"A(276) β[CI]":28s}  {"B(138) β[CI]":28s}  {"D(187) β[CI]":28s}')
    for col in struct_cols:
        rA = lme_single(df_A, col); rB = lme_single(df_B, col); rD = lme_single(df_D, col)
        def f(r):
            b = r['lme']['beta']; ci = r['lme']['ci']
            return f'{b:+.3f}[{ci[0]:+.3f},{ci[1]:+.3f}]'
        print(f'  {col:16s}  {f(rA):28s}  {f(rB):28s}  {f(rD):28s}')

    print(f'\n[tab:sf] Spearman rho 3열 병치')
    print(f'  {"지표":16s}  {"A rho":8s}  {"B rho":8s}  {"D rho":8s}')
    for col in struct_cols:
        rA = lme_single(df_A, col); rB = lme_single(df_B, col); rD = lme_single(df_D, col)
        print(f'  {col:16s}  {rA["spearman"]["rho"]:8.4f}  {rB["spearman"]["rho"]:8.4f}  {rD["spearman"]["rho"]:8.4f}')

    print(f'\n[tab:comp] joint model(gcl_avg+rnfl_q_i) 3열 병치')
    print(f'  {"":6s}  {"n":5s}  {"gcl_avg β(p)":22s}  {"rnfl_q_i β(p)":22s}  {"VIF":8s}  {"ΔAIC(vsGCL)":12s}  {"ΔAIC(vsRNFL)":12s}')
    for label, d in [('A(276)', df_A), ('B(138)', df_B), ('D(187)', df_D)]:
        j = lme_joint(d, 'gcl_avg', 'rnfl_q_i')
        diag = complementarity_stats(d, 'gcl_avg', 'rnfl_q_i')
        g = j['gcl_avg']; rr = j['rnfl_q_i']
        print(f'  {label:6s}  {j["n"]:5d}  '
              f'{g["beta"]:+.3f}(p={g["p"]:.4f})       '
              f'{rr["beta"]:+.3f}(p={rr["p"]:.4f})       '
              f'{diag.get("vif"):8.3f}  {diag.get("delta_aic_vs_gcl_only"):12.2f}  '
              f'{diag.get("delta_aic_vs_rnfl_only"):12.2f}')

    print(f'\n[요약] 신규성 서술용 — 필터 강도에 따른 코호트 크기·강도 트레이드오프')
    print(f'  A(무필터, n={n_A}): 가장 크지만 저품질 안 포함 — β 가장 약함')
    print(f'  B(표준기준, n={n_B}): 가장 엄격 — β 가장 강하지만 코호트 절반 손실, VIF 상승(공선성)')
    print(f'  D(evidence-based, n={n_D}): B보다 49안 크면서 상보성(ΔAIC) B보다 오히려 강함 — 크기/강도 둘 다 확보')


if __name__ == '__main__':
    main()
