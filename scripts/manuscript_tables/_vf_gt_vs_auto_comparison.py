# -*- coding: utf-8 -*-
"""
275 코호트에서 GT(수동검수본) vs 자동추출본으로 구조-기능 분석(tab:sf/comp)을
각각 돌려 추출오차가 결론을 바꾸는지 정량화.

y(vf_mean)는 두 런 모두 GT 고정 — analysis_master.csv 값 사용(VF는 이미 98.9%
정확도라 y까지 바꾸면 어느 쪽 오차 때문인지 분리 안 됨. x(OCT 5개 지표)만 GT/자동
스왑). 좌표계 컨벤션(GCA sup/inf swap, RNFL Method B) 은 phase3에서 이미 검증된
방향 그대로 적용 — 안 하면 라벨링 관례 차이를 추출오차로 오인함.

275 base 유지. data/analysis_master.csv 등 기존 파이프라인 파일은 읽기만 함.
"""
import csv, sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
from phase5_4_gcl_rnfl_compare import lme_single, lme_joint, complementarity_stats  # noqa: E402
import pandas as pd


def to_f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# GCA swap 방향 (phase3_extraction_accuracy.py 자동판별 결과 그대로 재사용)
GCA_SWAP_MAP = {'sup_t': 'sup_n', 'sup_n': 'sup_t', 'inf_t': 'inf_n', 'inf_n': 'inf_t'}
OD_GCA_SWAP, OS_GCA_SWAP = True, False   # phase3 auto-detected: A vs B 중 B(54.6%) 채택

# RNFL Method B: OS행 T/N 스왑
RNFL_TN_SWAP = {'rnfl_q_t': 'rnfl_q_n', 'rnfl_q_n': 'rnfl_q_t'}


def auto_gca_col(base_suffix, eye):
    """base_suffix: 'sup','sup_t','inf_t','inf','inf_n','sup_n' -> raw 컬럼명(스왑 적용)."""
    swap = OD_GCA_SWAP if eye == 'OD' else OS_GCA_SWAP
    key = base_suffix
    if swap and key in GCA_SWAP_MAP:
        key = GCA_SWAP_MAP[key]
    prefix = 'od_s_' if eye == 'OD' else 'os_s_'
    return prefix + key


def auto_rnfl_q_col(q, eye):
    col = f'rnfl_q_{q}'
    if eye == 'OS' and col in RNFL_TN_SWAP:
        return RNFL_TN_SWAP[col]
    return col


def build_bad_D():
    """코호트 D 확정 기준: FP>15 OR FN>20 OR ss<6 (FL·배너 제외기준 아님)."""
    rel = {(r['patient_id'], r['eye'], r['vf_date']): r
           for r in csv.DictReader(open(_ROOT / 'vf_reliability.csv', encoding='utf-8-sig'))}
    am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    bad = {}
    for r in am:
        pid, eye, vfdate = r['patient_id'], r['eye'], r['vf_date']
        ss_gca  = to_f(r['ss_gca_od'])  if eye == 'OD' else to_f(r['ss_gca_os'])
        ss_rnfl = to_f(r['ss_rnfl_od']) if eye == 'OD' else to_f(r['ss_rnfl_os'])
        bad_ss = (ss_gca is not None and ss_gca < 6) or (ss_rnfl is not None and ss_rnfl < 6)
        rv = rel.get((pid, eye, vfdate))
        if rv is None:
            bad[(pid, eye)] = bad_ss
            continue
        fp = to_f(rv['false_pos']); fn = to_f(rv['false_neg'])
        b_fp = fp is not None and fp > 15
        b_fn = fn is not None and fn > 20
        bad[(pid, eye)] = b_fp or b_fn or bad_ss
    return bad


def main(cohort='275'):
    am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    octv = {(r['patient_id'], r['eye'], r['oct_date']): r
            for r in csv.DictReader(open(_ROOT / 'data' / 'oct_values_raw.csv', encoding='utf-8-sig'))}
    rnfl = {(r['patient_id'], r['eye'], r['oct_date']): r
            for r in csv.DictReader(open(_ROOT / 'data' / 'rnfl_detail_raw.csv', encoding='utf-8-sig'))}

    if cohort == 'D':
        bad_D = build_bad_D()
        am = [r for r in am if not bad_D.get((r['patient_id'], r['eye']), False)]
        print(f'코호트 D 필터 적용: {len(am)}안')

    gt_records, auto_records = [], []
    n_miss = 0
    for r in am:
        pid, eye, octd = r['patient_id'], r['eye'], r['oct_date']
        vf_mean = to_f(r.get('vf_mean'))
        e = eye.lower()

        # ---- GT (기존 analysis_master.csv 값) ----
        gcl_avg_gt = to_f(r.get(f'avg_gcl_{e}'))
        gcl_min_gt = to_f(r.get(f'min_gcl_{e}'))
        qs = to_f(r.get('rnfl_q_s')); qi = to_f(r.get('rnfl_q_i'))
        qt = to_f(r.get('rnfl_q_t')); qn = to_f(r.get('rnfl_q_n'))
        qv = [v for v in [qs, qi, qt, qn] if v is not None]
        quad_mean_gt = sum(qv) / len(qv) if qv else None
        ch_vals = [to_f(r.get(f'rnfl_h{h:02d}')) for h in range(1, 13)]
        ch_ok = [v for v in ch_vals if v is not None]
        ch_mean_gt = sum(ch_ok) / len(ch_ok) if ch_ok else None

        gt_records.append({'pid': pid, 'eye': eye, 'vf_mean': vf_mean,
                            'gcl_avg': gcl_avg_gt, 'gcl_min': gcl_min_gt,
                            'rnfl_q_i': qi, 'rnfl_quad_mean': quad_mean_gt,
                            'rnfl_ch_mean': ch_mean_gt})

        # ---- 자동 추출본 ----
        key = (pid, eye, octd)
        ov = octv.get(key)
        rv = rnfl.get(key)
        if ov is None or rv is None:
            n_miss += 1
            continue

        gcl_avg_auto = to_f(ov.get(f'avg_gcl_{e}'))
        gcl_min_auto = to_f(ov.get(f'min_gcl_{e}'))

        qi_auto = to_f(rv.get(auto_rnfl_q_col('i', eye)))
        qs_auto = to_f(rv.get(auto_rnfl_q_col('s', eye)))
        qt_auto = to_f(rv.get(auto_rnfl_q_col('t', eye)))
        qn_auto = to_f(rv.get(auto_rnfl_q_col('n', eye)))
        qv_auto = [v for v in [qs_auto, qi_auto, qt_auto, qn_auto] if v is not None]
        quad_mean_auto = sum(qv_auto) / len(qv_auto) if qv_auto else None

        ch_vals_auto = [to_f(rv.get(f'rnfl_h{h:02d}')) for h in range(1, 13)]
        ch_ok_auto = [v for v in ch_vals_auto if v is not None]
        ch_mean_auto = sum(ch_ok_auto) / len(ch_ok_auto) if ch_ok_auto else None

        auto_records.append({'pid': pid, 'eye': eye, 'vf_mean': vf_mean,
                              'gcl_avg': gcl_avg_auto, 'gcl_min': gcl_min_auto,
                              'rnfl_q_i': qi_auto, 'rnfl_quad_mean': quad_mean_auto,
                              'rnfl_ch_mean': ch_mean_auto})

    print(f'GT n={len(gt_records)}  자동추출 매칭 n={len(auto_records)}  (미매칭 {n_miss})')

    df_gt = pd.DataFrame(gt_records)
    df_auto = pd.DataFrame(auto_records)

    struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']

    print('\n' + '=' * 100)
    print('Head-to-head LME β : GT vs 자동추출 (275 base, y=vf_mean은 GT 고정)')
    print('=' * 100)
    print(f'  {"지표":16s}  {"GT β [CI]":30s}  {"자동 β [CI]":30s}  {"Δβ":8s}  {"CI겹침":6s}')
    print('  ' + '-' * 100)

    summary = []
    for col in struct_cols:
        r_gt = lme_single(df_gt, col)
        r_auto = lme_single(df_auto, col)

        def f(r):
            if 'lme' not in r:
                return 'n/a', None, None
            b = r['lme']['beta']; ci = r['lme']['ci']
            return f'{b:+.3f}[{ci[0]:+.3f},{ci[1]:+.3f}]', b, ci

        s_gt, b_gt, ci_gt = f(r_gt)
        s_auto, b_auto, ci_auto = f(r_auto)
        if b_gt is not None and b_auto is not None:
            d = b_auto - b_gt
            overlap = ci_gt[0] < ci_auto[1] and ci_auto[0] < ci_gt[1]
        else:
            d, overlap = None, None
        summary.append((col, b_gt, b_auto, d, overlap))
        print(f'  {col:16s}  {s_gt:30s}  {s_auto:30s}  '
              f'{f"{d:+.3f}" if d is not None else "?":8s}  {"YES" if overlap else "NO":6s}')

    print('\n' + '=' * 90)
    print('상보성(gcl_avg + rnfl_q_i) : GT vs 자동추출')
    print('=' * 90)
    for label, df in [('GT', df_gt), ('자동추출', df_auto)]:
        j = lme_joint(df, 'gcl_avg', 'rnfl_q_i')
        diag = complementarity_stats(df, 'gcl_avg', 'rnfl_q_i')
        if 'gcl_avg' in j and j.get('gcl_avg'):
            g = j['gcl_avg']; r_ = j['rnfl_q_i']
            print(f'  {label:8s}  n={j["n"]}  gcl_avg β={g["beta"]:+.3f}(p={g["p"]:.5f})  '
                  f'rnfl_q_i β={r_["beta"]:+.3f}(p={r_["p"]:.5f})  '
                  f'VIF={diag.get("vif")}  ΔAIC(vsGCL)={diag.get("delta_aic_vs_gcl_only")}  '
                  f'ΔAIC(vsRNFL)={diag.get("delta_aic_vs_rnfl_only")}')

    print('\n' + '=' * 70)
    print('핵심 판정')
    print('=' * 70)
    all_small = all(d is not None and abs(d) < 0.03 for _, _, _, d, _ in summary)
    all_overlap = all(ov for _, _, _, _, ov in summary if ov is not None)
    for col, b_gt, b_auto, d, overlap in summary:
        if d is None:
            print(f'  {col:16s}  비교불가(자동추출본 n=0 또는 실패) — GT β={b_gt}')
            continue
        flag = ''
        if abs(d) >= 0.03 or not overlap:
            flag = '  <<< 벌어짐(주의)'
        print(f'  {col:16s}  Δβ={d:+.3f}  CI겹침={overlap}{flag}')
    print(f'\n전체 판정: 모든|Δβ|<0.03 이고 CI 전부 겹침 = {all_small and all_overlap}')
    if all_small and all_overlap:
        print('→ "추출 오차가 구조-기능 추정을 편향시키지 않는다" 성립')
    else:
        print('→ 일부 지표에서 벌어짐 — 아래에서 원인(추출정확도) 연결')


if __name__ == '__main__':
    cohort = sys.argv[1] if len(sys.argv) > 1 else '275'
    main(cohort)
