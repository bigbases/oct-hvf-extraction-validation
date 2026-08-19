# -*- coding: utf-8 -*-
"""
코호트 D 최종 확정 + 파일 저장(휘발 방지).
입력: data/analysis_master.csv(275, gap<=90일 이미 만족) + vf_reliability.csv(340)
제외: FP>15% OR FN>20% OR ss<6(GCA 또는 RNFL). FL·배너는 제외기준 아님.
출력: cohort_final_D.csv (patient_id,eye,vf_date 키만 — analysis_master.csv 전체 컬럼은
      필요시 이 키로 다시 join해서 사용. 이 파일 자체가 "확정된 187안 목록".)
"""
import csv

def to_f(v):
    try: return float(v)
    except (TypeError, ValueError): return None

rel = {(r['patient_id'], r['eye'], r['vf_date']): r
       for r in csv.DictReader(open('vf_reliability.csv', encoding='utf-8-sig'))}
am = list(csv.DictReader(open('data/analysis_master.csv', encoding='utf-8')))
print(f'입력: analysis_master.csv {len(am)}행, vf_reliability.csv {len(rel)}고유키')

kept = []
excluded = []
n_unmatched = 0
for r in am:
    pid, eye, vfdate = r['patient_id'], r['eye'], r['vf_date']
    ss_gca  = to_f(r['ss_gca_od'])  if eye == 'OD' else to_f(r['ss_gca_os'])
    ss_rnfl = to_f(r['ss_rnfl_od']) if eye == 'OD' else to_f(r['ss_rnfl_os'])
    bad_ss = (ss_gca is not None and ss_gca < 6) or (ss_rnfl is not None and ss_rnfl < 6)

    rv = rel.get((pid, eye, vfdate))
    if rv is None:
        n_unmatched += 1
        bad = bad_ss
        reason = 'ss<6' if bad_ss else ''
    else:
        fp = to_f(rv['false_pos']); fn = to_f(rv['false_neg'])
        b_fp = fp is not None and fp > 15
        b_fn = fn is not None and fn > 20
        bad = b_fp or b_fn or bad_ss
        reasons = []
        if b_fp: reasons.append('FP>15')
        if b_fn: reasons.append('FN>20')
        if bad_ss: reasons.append('ss<6')
        reason = '+'.join(reasons)

    if bad:
        excluded.append((pid, eye, vfdate, reason))
    else:
        kept.append((pid, eye, vfdate))

print(f'vf_reliability 미매칭(ss만 적용): {n_unmatched}안')
print(f'제외: {len(excluded)}안')
print(f'유지(코호트 D): {len(kept)}안')

with open('cohort_final_D.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['patient_id', 'eye', 'vf_date'])
    for row in kept:
        w.writerow(row)

with open('cohort_final_D_excluded.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['patient_id', 'eye', 'vf_date', 'exclude_reason'])
    for row in excluded:
        w.writerow(row)

print(f'\n저장: cohort_final_D.csv ({len(kept)}행), cohort_final_D_excluded.csv ({len(excluded)}행)')
print(f'검증: {len(kept)} + {len(excluded)} = {len(kept)+len(excluded)} (원본 {len(am)}과 일치해야 함: {len(kept)+len(excluded)==len(am)})')
