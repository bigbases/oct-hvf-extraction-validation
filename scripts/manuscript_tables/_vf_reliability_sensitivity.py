# -*- coding: utf-8 -*-
"""
VF 신뢰도(FL/FP/FN/배너) + signal strength<6 제외 시 tab:sf 핵심 상관(gcl_avg,
gcl_min, rnfl_quad_mean vs vf_mean)의 beta가 바뀌는지 확인하는 민감도 분석.
phase5_1_correlations.py의 _corr_block을 그대로 재사용(로직 중복 방지).
data/analysis_master.csv는 건드리지 않음 — 필터링은 이 스크립트 내에서만.
"""
import csv, sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
sys.path.insert(0, str(_ROOT / 'src'))

from phase5_1_correlations import _corr_block, _f  # noqa: E402


def to_f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    rel = {(r['patient_id'], r['eye'], r['vf_date']): r
           for r in csv.DictReader(open('vf_reliability.csv', encoding='utf-8-sig'))}
    rows = list(csv.DictReader(open('data/analysis_master.csv', encoding='utf-8')))

    for r in rows:
        if r['eye'] == 'OD':
            r['_gcl_avg'] = _f(r['avg_gcl_od'])
            r['_gcl_min'] = _f(r['min_gcl_od'])
            r['_ss_gca']  = to_f(r['ss_gca_od'])
            r['_ss_rnfl'] = to_f(r['ss_rnfl_od'])
        else:
            r['_gcl_avg'] = _f(r['avg_gcl_os'])
            r['_gcl_min'] = _f(r['min_gcl_os'])
            r['_ss_gca']  = to_f(r['ss_gca_os'])
            r['_ss_rnfl'] = to_f(r['ss_rnfl_os'])

        rnfl_vals = [_f(r[f'rnfl_q_{q}']) for q in ('s', 't', 'i', 'n')]
        rnfl_ok   = [v for v in rnfl_vals if v is not None]
        r['_rnfl_mean'] = sum(rnfl_ok) / len(rnfl_ok) if len(rnfl_ok) >= 2 else None

        r['_vf']  = _f(r['vf_mean'])
        r['_pid'] = r['patient_id']

        # signal strength 필터는 vf_reliability 매칭 여부와 무관하게 276안 전체 적용 가능
        bad_ss = (r['_ss_gca'] is not None and r['_ss_gca'] < 6) or \
                 (r['_ss_rnfl'] is not None and r['_ss_rnfl'] < 6)

        key = (r['patient_id'], r['eye'], r['vf_date'])
        rv = rel.get(key)
        if rv is None:
            r['_vf_reliability_known'] = False
            r['_bad'] = bad_ss   # VF 신뢰도 정보 없는 6안은 그 부분만 "판정 불가"로 두고 ss 필터는 그대로 적용
        else:
            r['_vf_reliability_known'] = True
            num, den = to_f(rv['fix_loss_num']), to_f(rv['fix_loss_denom'])
            fl = 100 * num / den if (num is not None and den) else None
            fp = to_f(rv['false_pos'])
            fn = to_f(rv['false_neg'])
            low = rv['low_reliability']
            bad_vf = (fl is not None and fl > 33) or (fp is not None and fp > 15) \
                or (fn is not None and fn > 20) or (low == 'Y')
            r['_bad'] = bad_vf or bad_ss

    all_rows = rows
    clean_rows = [r for r in rows if not r['_bad']]
    known = sum(1 for r in rows if r['_vf_reliability_known'])

    print('=' * 70)
    print(f'전체(baseline) n={len(all_rows)}   |   VF신뢰도 정보 확보 n={known}/276')
    print(f'필터 통과(clean) n={len(clean_rows)}   |   제외 n={len(all_rows)-len(clean_rows)}')
    print('=' * 70)

    def cmp(pred_key, label):
        x_all = [r[pred_key] for r in all_rows]
        y_all = [r['_vf'] for r in all_rows]
        x_cln = [r[pred_key] for r in clean_rows]
        y_cln = [r['_vf'] for r in clean_rows]
        b_all = _corr_block(x_all, y_all, label + ' (all 276)')
        b_cln = _corr_block(x_cln, y_cln, label + ' (clean)')
        ra = b_all['spearman']['rho']; rc = b_cln['spearman']['rho']
        pa = b_all['pearson']['r'];    pc = b_cln['pearson']['r']
        print(f'\n[{label}]')
        print(f'  all276  n={b_all["n"]:3d}  Spearman rho={ra:.4f}  Pearson r={pa:.4f}')
        print(f'  clean   n={b_cln["n"]:3d}  Spearman rho={rc:.4f}  Pearson r={pc:.4f}')
        print(f'  Δrho={rc-ra:+.4f}   Δr={pc-pa:+.4f}')

    cmp('_gcl_avg', 'avg_gcl vs vf_mean')
    cmp('_gcl_min', 'min_gcl vs vf_mean')
    cmp('_rnfl_mean', 'rnfl_quad_mean vs vf_mean')


if __name__ == '__main__':
    main()
