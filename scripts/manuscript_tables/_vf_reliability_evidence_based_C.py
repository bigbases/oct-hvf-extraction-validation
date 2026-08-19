# -*- coding: utf-8 -*-
"""
코호트 A(무필터 276) / B(표준기준 138) / C(evidence-based: FP>15% 또는 배너Y만,
FL·FN 제외기준에서 뺌; ss<6는 세 코호트 공통) 3열 병치.
Yohannan 2017 근거: FL은 established glaucoma에서 신뢰도 지표로 약함,
FN은 중증에서 병 자체를 반영. FP만 전 단계에서 확실한 지표.

기존 276 파이프라인(data/analysis_master.csv, scripts/phase5_*.py)은 읽기만
하고 수정하지 않음. 세 코호트 중 무엇도 "최종"으로 하드코딩하지 않음 —
이 스크립트 실행 시점의 탐색 산출물일 뿐.
"""
import csv, sys, json
from pathlib import Path
from collections import Counter

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
from phase5_4_gcl_rnfl_compare import (  # noqa: E402
    load_df, lme_single, lme_joint, complementarity_stats
)


def to_f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def build_eye_metrics():
    """(patient_id, eye) -> dict(fl_pct, fp, fn, low, ss_gca, ss_rnfl, bad_ss,
    vf_reliability_known). analysis_master 자신의 vf_date로만 매칭(교차-방문 오염 방지)."""
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
            out[(pid, eye)] = dict(fl=None, fp=None, fn=None, low=None,
                                    bad_ss=bad_ss, known=False)
            continue
        num, den = to_f(rv['fix_loss_num']), to_f(rv['fix_loss_denom'])
        fl = 100 * num / den if (num is not None and den) else None
        fp = to_f(rv['false_pos'])
        fn = to_f(rv['false_neg'])
        low = rv['low_reliability']
        out[(pid, eye)] = dict(fl=fl, fp=fp, fn=fn, low=low, bad_ss=bad_ss, known=True)
    return out


def main():
    m = build_eye_metrics()

    bad_B, bad_C = {}, {}
    reason_B = {}  # eye -> set of reasons it failed B
    for k, d in m.items():
        b_fl  = d['fl'] is not None and d['fl'] > 33
        b_fp  = d['fp'] is not None and d['fp'] > 15
        b_fn  = d['fn'] is not None and d['fn'] > 20
        b_low = d['low'] == 'Y'
        bad_B[k] = b_fl or b_fp or b_fn or b_low or d['bad_ss']
        bad_C[k] = b_fp or b_low or d['bad_ss']
        reasons = set()
        if b_fl:  reasons.add('FL')
        if b_fp:  reasons.add('FP')
        if b_fn:  reasons.add('FN')
        if b_low: reasons.add('banner')
        if d['bad_ss']: reasons.add('ss')
        reason_B[k] = reasons

    n_A = len(m)
    n_B = sum(1 for v in bad_B.values() if not v)
    n_C = sum(1 for v in bad_C.values() if not v)
    n_pat_A = len({k[0] for k in m})
    n_pat_B = len({k[0] for k, v in bad_B.items() if not v})
    n_pat_C = len({k[0] for k, v in bad_C.items() if not v})

    print('=' * 78)
    print(f'코호트 크기:  A(무필터)={n_A}안/{n_pat_A}명   '
          f'B(표준)={n_B}안/{n_pat_B}명   C(evidence-based)={n_C}안/{n_pat_C}명')
    print('=' * 78)

    # B에서 빠졌다가 C에서 살아난 안
    restored = [k for k in m if bad_B[k] and not bad_C[k]]
    # C에서도 여전히 빠진 안 (B에서도 빠졌던 것 중)
    still_out = [k for k in m if bad_B[k] and bad_C[k]]
    print(f'\nB→C 복원: {len(restored)}안 (B에서 제외됐다가 C에서 살아남)')
    print(f'B에서도 C에서도 제외 유지: {len(still_out)}안')

    # 복원된 안들의 원래 제외사유 분해
    reason_counter = Counter()
    for k in restored:
        reasons = reason_B[k] - {'ss'}  # ss는 C에도 공통 적용되므로 복원 안 됐을 것 — 확인용 제외
        reason_counter[tuple(sorted(reasons))] += 1
    print('\nB→C 복원 안의 원래(B기준) 제외사유 조합:')
    for combo, cnt in reason_counter.most_common():
        print(f'  {"+".join(combo):20s}  {cnt}안')

    # 배너 Y ∩ FP<=15% (FL 완화-배너 충돌 점검)
    banner_y = [k for k, d in m.items() if d['low'] == 'Y']
    banner_y_fp_ok = [k for k in banner_y if m[k]['fp'] is not None and m[k]['fp'] <= 15]
    print(f'\n배너 Y 전체: {len(banner_y)}안   그중 FP<=15%(FL/FN 문제로 배너 뜬 것으로 추정): '
          f'{len(banner_y_fp_ok)}안')
    print('  → 이 안들은 evidence-based 기준에서도 "배너 Y"로 여전히 제외됨(배너는 유지 원칙이므로 충돌 없음)')

    # ══════════════════════════════════════════════════════════════
    # STEP 2+3: tab:sf/comp + Pearson r, A/B/C 3열
    # ══════════════════════════════════════════════════════════════
    df = load_df()
    df['bad_B'] = df.apply(lambda r: bad_B.get((r['pid'], r['eye']), False), axis=1)
    df['bad_C'] = df.apply(lambda r: bad_C.get((r['pid'], r['eye']), False), axis=1)
    df_A = df
    df_B = df[~df['bad_B']].copy()
    df_C = df[~df['bad_C']].copy()

    print('\n' + '=' * 90)
    print('STEP 2: Head-to-head LME β  (A / B / C)')
    print('=' * 90)
    struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']
    print(f'  {"지표":16s}  {"A(276) β[CI]":26s}  {"B(138) β[CI]":26s}  {"C β[CI]":26s}')
    print('  ' + '-' * 100)
    h2h = {}
    for col in struct_cols:
        rA = lme_single(df_A, col); rB = lme_single(df_B, col); rC = lme_single(df_C, col)
        h2h[col] = {'A': rA, 'B': rB, 'C': rC}
        def f(r):
            if 'lme' not in r: return 'n/a'
            b=r['lme']['beta']; ci=r['lme']['ci']
            return f'{b:+.3f}[{ci[0]:+.3f},{ci[1]:+.3f}]'
        print(f'  {col:16s}  {f(rA):26s}  {f(rB):26s}  {f(rC):26s}')

    print('\n상보성(gcl_avg + rnfl_q_i): A/B/C')
    joint = {}
    for lbl, d in [('A', df_A), ('B', df_B), ('C', df_C)]:
        j = lme_joint(d, 'gcl_avg', 'rnfl_q_i')
        diag = complementarity_stats(d, 'gcl_avg', 'rnfl_q_i')
        joint[lbl] = {'joint': j, 'diag': diag}
        if 'gcl_avg' in j and j.get('gcl_avg'):
            g=j['gcl_avg']; rr=j['rnfl_q_i']
            print(f'  {lbl}: n={j["n"]}  gcl_avg β={g["beta"]:+.3f}(p={g["p"]:.5f})  '
                  f'rnfl_q_i β={rr["beta"]:+.3f}(p={rr["p"]:.5f})  '
                  f'VIF={diag.get("vif")}  ΔAIC(vsGCL)={diag.get("delta_aic_vs_gcl_only")}  '
                  f'ΔAIC(vsRNFL)={diag.get("delta_aic_vs_rnfl_only")}')

    print('\n' + '=' * 78)
    print('STEP 3: Pearson r  (A / B / C)')
    print('=' * 78)
    for col in ['gcl_avg', 'gcl_min', 'rnfl_quad_mean']:
        rA = lme_single(df_A, col)['pearson']['r']
        rB = lme_single(df_B, col)['pearson']['r']
        rC = lme_single(df_C, col)['pearson']['r']
        print(f'  {col:16s}  A={rA:.4f}   B={rB:.4f}   C={rC:.4f}')

    out = {
        'n': {'A': n_A, 'B': n_B, 'C': n_C},
        'n_patients': {'A': n_pat_A, 'B': n_pat_B, 'C': n_pat_C},
        'restored_B_to_C': len(restored),
        'restored_reason_breakdown': {'+'.join(k) if k else 'none': v for k, v in reason_counter.items()},
        'banner_y_total': len(banner_y),
        'banner_y_fp_ok': len(banner_y_fp_ok),
        'head_to_head': h2h,
        'complementarity': joint,
    }
    with open(_ROOT / '_vf_reliability_cohortC_result.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)
    print('\n저장: _vf_reliability_cohortC_result.json')


if __name__ == '__main__':
    main()
