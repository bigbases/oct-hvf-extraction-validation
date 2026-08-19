# -*- coding: utf-8 -*-
"""
STEP 0: clean-138(VF 신뢰도+signal strength 필터 통과)에서 tab:sf/tab:comp
전체(LME beta 5개, joint beta, VIF, ΔAIC)를 재계산해 276안 값과 병치.

phase5_4_gcl_rnfl_compare.py의 함수(load_df/lme_single/lme_joint/
complementarity_stats)를 그대로 재사용 — 로직 중복 금지, data/analysis_master.csv
는 건드리지 않음. 138 필터는 이 스크립트 내 별도 변수로만 존재(하드코딩 아님,
탐색용 산출물).
"""
import csv, sys, json
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))

from phase5_4_gcl_rnfl_compare import (  # noqa: E402
    load_df, lme_single, lme_joint, complementarity_stats, get_severity_groups
)


def to_f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def build_bad_set():
    """(patient_id, eye) -> True(제외 대상) / False(clean). vf_reliability 미매칭 6안은
    ss 기준만으로 판정(그 부분은 STEP0 출력에 별도 명시)."""
    rel = {(r['patient_id'], r['eye'], r['vf_date']): r
           for r in csv.DictReader(open(_ROOT / 'vf_reliability.csv', encoding='utf-8-sig'))}
    am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    bad = {}
    unknown = []
    for r in am:
        eye = r['eye']
        ss_gca  = to_f(r['ss_gca_od'])  if eye == 'OD' else to_f(r['ss_gca_os'])
        ss_rnfl = to_f(r['ss_rnfl_od']) if eye == 'OD' else to_f(r['ss_rnfl_os'])
        bad_ss = (ss_gca is not None and ss_gca < 6) or (ss_rnfl is not None and ss_rnfl < 6)

        key = (r['patient_id'], eye, r['vf_date'])
        rv = rel.get(key)
        if rv is None:
            bad[(r['patient_id'], eye)] = bad_ss
            unknown.append((r['patient_id'], eye))
            continue
        num, den = to_f(rv['fix_loss_num']), to_f(rv['fix_loss_denom'])
        fl = 100 * num / den if (num is not None and den) else None
        fp = to_f(rv['false_pos'])
        fn = to_f(rv['false_neg'])
        low = rv['low_reliability']
        bad_vf = (fl is not None and fl > 33) or (fp is not None and fp > 15) \
            or (fn is not None and fn > 20) or (low == 'Y')
        bad[(r['patient_id'], eye)] = bad_vf or bad_ss
    return bad, unknown


def fmt_lme(res):
    if 'lme' not in res:
        return 'n/a'
    b = res['lme']['beta']; ci = res['lme']['ci']; p = res['lme']['p']
    return f'β={b:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}] p={p:.4f}'


def main():
    bad, unknown = build_bad_set()
    df = load_df()
    df['bad'] = df.apply(lambda r: bad.get((r['pid'], r['eye']), False), axis=1)
    df_clean = df[~df['bad']].copy()

    print('=' * 78)
    print(f'전체(276) n={len(df)}   clean(138 잠정) n={len(df_clean)}   '
          f'(VF신뢰도 미매칭 {len(unknown)}안은 ss 기준만 적용)')
    print('=' * 78)

    struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']

    print(f'\n[ 1. Head-to-head LME β : 276 vs clean ]')
    print(f'  {"지표":16s}  {"276 β [CI]":32s}  {"clean β [CI]":32s}  Δβ')
    print('  ' + '-' * 100)
    h2h_all = {}
    h2h_cln = {}
    for col in struct_cols:
        r_all = lme_single(df, col)
        r_cln = lme_single(df_clean, col)
        h2h_all[col] = r_all
        h2h_cln[col] = r_cln
        b_all = r_all.get('lme', {}).get('beta')
        b_cln = r_cln.get('lme', {}).get('beta')
        d = (b_cln - b_all) if (b_all is not None and b_cln is not None) else None
        print(f'  {col:16s}  {fmt_lme(r_all):32s}  {fmt_lme(r_cln):32s}  '
              f'{f"{d:+.3f}" if d is not None else "?"}')

    print(f'\n[ 2. 상보성(gcl_avg + rnfl_q_i) : 276 vs clean ]')
    joint_all = lme_joint(df, 'gcl_avg', 'rnfl_q_i')
    joint_cln = lme_joint(df_clean, 'gcl_avg', 'rnfl_q_i')
    diag_all = complementarity_stats(df, 'gcl_avg', 'rnfl_q_i')
    diag_cln = complementarity_stats(df_clean, 'gcl_avg', 'rnfl_q_i')

    def show_joint(label, joint, diag):
        print(f'\n  -- {label} --')
        if 'gcl_avg' in joint and joint['gcl_avg']:
            g = joint['gcl_avg']; r = joint['rnfl_q_i']
            print(f'  n={joint["n"]}  gcl_avg β={g["beta"]:+.3f} p={g["p"]:.5f}   '
                  f'rnfl_q_i β={r["beta"]:+.3f} p={r["p"]:.5f}')
        if 'vif' in diag:
            print(f'  VIF={diag["vif"]}  predictor_corr_r={diag.get("predictor_corr_r")}')
            print(f'  AIC gcl-only={diag.get("aic_gcl_only")}  rnfl-only={diag.get("aic_rnfl_only")}  '
                  f'joint={diag.get("aic_joint")}')
            print(f'  ΔAIC(vs gcl-only)={diag.get("delta_aic_vs_gcl_only")}   '
                  f'ΔAIC(vs rnfl-only)={diag.get("delta_aic_vs_rnfl_only")}')

    show_joint('276 (baseline)', joint_all, diag_all)
    show_joint('clean-138 (잠정)', joint_cln, diag_cln)

    # 판정 기준: 두 β 모두 유의(p<0.05), VIF 낮음(<5 통상 기준), ΔAIC 둘 다 양수
    def verdict(joint, diag):
        if 'gcl_avg' not in joint or not joint.get('gcl_avg'):
            return 'UNKNOWN(모델 실패)'
        g_p = joint['gcl_avg']['p']; r_p = joint['rnfl_q_i']['p']
        vif = diag.get('vif', 99)
        da_g = diag.get('delta_aic_vs_gcl_only', -1)
        da_r = diag.get('delta_aic_vs_rnfl_only', -1)
        sig_both = g_p < 0.05 and r_p < 0.05
        low_vif = vif < 5
        aic_better = da_g > 0 and da_r > 0
        if sig_both and low_vif and aic_better:
            return 'COMPLEMENTARITY HOLDS'
        return f'COMPLEMENTARITY QUESTIONABLE (sig_both={sig_both}, low_vif={low_vif}({vif}), aic_better={aic_better})'

    v_all = verdict(joint_all, diag_all)
    v_cln = verdict(joint_cln, diag_cln)
    print(f'\n[ 판정 ]')
    print(f'  276:        {v_all}')
    print(f'  clean-138:  {v_cln}')

    out = {
        'n_all': len(df), 'n_clean': len(df_clean),
        'n_vf_reliability_unmatched': len(unknown),
        'head_to_head_all': {k: v for k, v in h2h_all.items()},
        'head_to_head_clean': {k: v for k, v in h2h_cln.items()},
        'complementarity_all': {'joint': joint_all, 'diag': diag_all},
        'complementarity_clean': {'joint': joint_cln, 'diag': diag_cln},
        'verdict_all': v_all, 'verdict_clean': v_cln,
    }
    with open(_ROOT / '_vf_reliability_step0_result.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)
    print(f'\n저장: _vf_reliability_step0_result.json')


if __name__ == '__main__':
    main()
