# -*- coding: utf-8 -*-
"""
원고용 최종 표 일괄 생성 — 코호트 D의 유일 출처는 cohort_final_D.csv(187행).
이 파일 외 다른 곳에서 코호트를 재정의하지 않음. vfacc만 별도로 324 pool 사용.
"""
import csv, sys, json
from pathlib import Path
from collections import Counter
import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
from phase5_4_gcl_rnfl_compare import load_df, lme_single, lme_joint, complementarity_stats  # noqa: E402


def to_f(v):
    try: return float(v)
    except (TypeError, ValueError): return None


# ── 코호트 D 유일 출처 로드 ──
D_KEYS = set()
with open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig') as f:
    for row in csv.reader(f):
        if row == ['patient_id', 'eye', 'vf_date']:
            continue
        D_KEYS.add(tuple(row))
print(f'[출처 확인] cohort_final_D.csv 로드: {len(D_KEYS)}안')
assert len(D_KEYS) == 187, f'cohort_final_D.csv 행수 이상: {len(D_KEYS)}'

am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
am_idx = {(r['patient_id'], r['eye'], r['vf_date']): r for r in am}
d_rows = [am_idx[k] for k in D_KEYS if k in am_idx]
print(f'[정합성] analysis_master.csv와 매칭: {len(d_rows)}/187')
assert len(d_rows) == 187, 'cohort_final_D.csv가 analysis_master.csv와 완전히 매칭되지 않음'

print('\n' + '=' * 90)
print('1. tab:cohort  (n=187, 출처: cohort_final_D.csv x data/analysis_master.csv)')
print('=' * 90)

n_patients = len({r['patient_id'] for r in d_rows})
print(f'환자 수: {n_patients}   안 수: {len(d_rows)}')

patients = {}
for r in d_rows:
    patients.setdefault(r['patient_id'], []).append(r)
sex_counter = Counter()
for pid, rows in patients.items():
    sexes = {r['sex'] for r in rows if r['sex']}
    if len(sexes) == 1:
        sex_counter[list(sexes)[0]] += 1
    elif len(sexes) == 0:
        sex_counter['NR'] += 1
    else:
        sex_counter['conflict'] += 1
print(f'성별(환자단위): Male={sex_counter.get("Male",0)}  Female={sex_counter.get("Female",0)}  NR={sex_counter.get("NR",0)}')

ages = [to_f(r['age_at_oct']) for r in d_rows if to_f(r['age_at_oct']) is not None]
print(f'나이(안단위): mean={np.mean(ages):.1f} SD={np.std(ages,ddof=1):.1f} '
      f'range=[{min(ages):.1f},{max(ages):.1f}] N={len(ages)} (결측 {len(d_rows)-len(ages)})')

ms = [to_f(r['vf_mean']) for r in d_rows if to_f(r['vf_mean']) is not None]
print(f'MS(vf_mean): mean={np.mean(ms):.2f} SD={np.std(ms,ddof=1):.2f} '
      f'range=[{min(ms):.1f},{max(ms):.1f}] N={len(ms)}')

sev = Counter(r['severity_tertile'] for r in d_rows if r['severity_tertile'])
print(f'severity: mild={sev.get("mild",0)}  moderate={sev.get("moderate",0)}  severe={sev.get("severe",0)}')

sita = {(r['patient_id'], r['eye'], r['vf_date']): r['sita_strategy']
        for r in csv.DictReader(open(_ROOT / 'data' / 'sita_per_eye.csv', encoding='utf-8-sig'))}
strat = Counter()
for r in d_rows:
    s = sita.get((r['patient_id'], r['eye'], r['vf_date']), 'NR')
    strat[s] += 1
print(f'VF strategy: {dict(strat)}')
print('instrument SW version: 계산 불가 — 어느 파일에도 체계적으로 기록되지 않음')

# ── df 구성(phase5_4 스타일) — cohort_final_D.csv로만 필터 ──
df_full = load_df()
am_key_by_pideye = {(r['patient_id'], r['eye']): (r['patient_id'], r['eye'], r['vf_date']) for r in am}
df_full['in_D'] = df_full.apply(lambda r: am_key_by_pideye.get((r['pid'], r['eye'])) in D_KEYS, axis=1)
df_D = df_full[df_full['in_D']].copy()
assert len(df_D) == 187

print('\n' + '=' * 90)
print('2. tab:sf  (n=187, 출처: cohort_final_D.csv)')
print('=' * 90)
struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']
sf_result = {}
for col in struct_cols:
    r = lme_single(df_D, col)
    sf_result[col] = r
    rho = r['spearman']['rho']; b = r['lme']['beta']; ci = r['lme']['ci']; p = r['lme']['p']
    print(f'  {col:16s}  n={r["n"]:3d}  rho={rho:.4f}  beta={b:+.4f}  CI=[{ci[0]:+.4f},{ci[1]:+.4f}]  p={p:.5f}')

print('\n' + '=' * 90)
print('3. tab:comp  (n=187, 출처: cohort_final_D.csv)')
print('=' * 90)
j = lme_joint(df_D, 'gcl_avg', 'rnfl_q_i')
diag = complementarity_stats(df_D, 'gcl_avg', 'rnfl_q_i')
g = j['gcl_avg']; rr = j['rnfl_q_i']
print(f'  n={j["n"]}  gcl_avg joint_beta={g["beta"]:+.4f}(p={g["p"]:.5f})  '
      f'rnfl_q_i joint_beta={rr["beta"]:+.4f}(p={rr["p"]:.5f})')
print(f'  VIF={diag.get("vif")}  deltaAIC(vsGCL)={diag.get("delta_aic_vs_gcl_only")}  '
      f'deltaAIC(vsRNFL)={diag.get("delta_aic_vs_rnfl_only")}')
solo_gcl = sf_result['gcl_avg']['lme']['beta']; solo_rnfl = sf_result['rnfl_q_i']['lme']['beta']
att_gcl = (solo_gcl - g['beta']) / solo_gcl * 100
att_rnfl = (solo_rnfl - rr['beta']) / solo_rnfl * 100
print(f'  solo->joint 감쇠율: gcl_avg {solo_gcl:.3f}->{g["beta"]:.3f} ({att_gcl:.1f}%)   '
      f'rnfl_q_i {solo_rnfl:.3f}->{rr["beta"]:.3f} ({att_rnfl:.1f}%)')

print('\n' + '=' * 90)
print('4. tab:autovsgt  (n=187, 출처: cohort_final_D.csv, y=vf_mean GT고정)')
print('=' * 90)
octv = {(r['patient_id'], r['eye'], r['oct_date']): r
        for r in csv.DictReader(open(_ROOT / 'data' / 'oct_values_raw.csv', encoding='utf-8-sig'))}
rnfl_raw = {(r['patient_id'], r['eye'], r['oct_date']): r
            for r in csv.DictReader(open(_ROOT / 'data' / 'rnfl_detail_raw.csv', encoding='utf-8-sig'))}
GCA_SWAP_MAP = {'sup_t': 'sup_n', 'sup_n': 'sup_t', 'inf_t': 'inf_n', 'inf_n': 'inf_t'}
RNFL_TN_SWAP = {'rnfl_q_t': 'rnfl_q_n', 'rnfl_q_n': 'rnfl_q_t'}

def auto_rnfl_q_col(q, eye):
    col = f'rnfl_q_{q}'
    return RNFL_TN_SWAP[col] if (eye == 'OS' and col in RNFL_TN_SWAP) else col

auto_records = []
n_auto_miss = 0
for k in D_KEYS:
    r = am_idx[k]
    pid, eye, octd = r['patient_id'], r['eye'], r['oct_date']
    vf_mean = to_f(r.get('vf_mean'))
    ov = octv.get((pid, eye, octd)); rv = rnfl_raw.get((pid, eye, octd))
    if ov is None or rv is None:
        n_auto_miss += 1
        continue
    e = eye.lower()
    gcl_avg_auto = to_f(ov.get(f'avg_gcl_{e}'))
    gcl_min_auto = to_f(ov.get(f'min_gcl_{e}'))
    qi = to_f(rv.get(auto_rnfl_q_col('i', eye))); qs = to_f(rv.get(auto_rnfl_q_col('s', eye)))
    qt = to_f(rv.get(auto_rnfl_q_col('t', eye))); qn = to_f(rv.get(auto_rnfl_q_col('n', eye)))
    qv = [v for v in [qs, qi, qt, qn] if v is not None]
    quad_mean_auto = sum(qv) / len(qv) if qv else None
    ch_vals = [to_f(rv.get(f'rnfl_h{h:02d}')) for h in range(1, 13)]
    ch_ok = [v for v in ch_vals if v is not None]
    ch_mean_auto = sum(ch_ok) / len(ch_ok) if ch_ok else None
    auto_records.append({'pid': pid, 'eye': eye, 'vf_mean': vf_mean,
                          'gcl_avg': gcl_avg_auto, 'gcl_min': gcl_min_auto,
                          'rnfl_q_i': qi, 'rnfl_quad_mean': quad_mean_auto,
                          'rnfl_ch_mean': ch_mean_auto})
import pandas as pd
df_auto = pd.DataFrame(auto_records)
print(f'  187안 중 자동추출 매칭: {len(df_auto)}  (미매칭 {n_auto_miss})')
for col in struct_cols:
    r_gt = sf_result[col]
    r_auto = lme_single(df_auto, col)
    b_gt = r_gt['lme']['beta']; ci_gt = r_gt['lme']['ci']
    if 'lme' not in r_auto:
        print(f'  {col:16s}  GT={b_gt:+.4f}  자동=비교불가(n부족)')
        continue
    b_auto = r_auto['lme']['beta']; ci_auto = r_auto['lme']['ci']
    d = b_auto - b_gt
    overlap = ci_gt[0] < ci_auto[1] and ci_auto[0] < ci_gt[1]
    print(f'  {col:16s}  GT_beta={b_gt:+.4f}[{ci_gt[0]:+.4f},{ci_gt[1]:+.4f}]  '
          f'auto_beta={b_auto:+.4f}[{ci_auto[0]:+.4f},{ci_auto[1]:+.4f}]  '
          f'delta={d:+.4f}  CI겹침={overlap}')

print('\n' + '=' * 90)
print('5. tab:vfacc  (n=324, VF validation pool — 187과 무관, 출처: results/phase3_extraction_accuracy.json)')
print('=' * 90)
r3 = json.load(open(_ROOT / 'results' / 'phase3_extraction_accuracy.json', encoding='utf-8'))
va = r3['vf_accuracy']['accuracy_normal_eyes']
print(f'  exact={va["exact_match_pct"]}%  within1dB={va["within1dB_pct"]}%  within2dB={va["within2dB_pct"]}%')
print(f'  MAE={va["mae"]}dB  RMSE={va["rmse"]}dB  missing_rate={va["missing_rate_raw_pct"]}%')
print(f'  유효포인트={va["n_valid"]}/{va["n_pairs"]}  안수=324(offset 0)')

print('\n' + '=' * 90)
print('6. tab:octacc  (GT 있는 280안 pool — 187/324 둘 다 아님, 출처: 동일 JSON)')
print('=' * 90)
ra = r3['rnfl_accuracy']['per_column']
ex_sum = n_sum = 0
for k in ['rnfl_q_s', 'rnfl_q_t', 'rnfl_q_i', 'rnfl_q_n']:
    ex_sum += ra[k]['exact_match_pct'] * ra[k]['n_valid']; n_sum += ra[k]['n_valid']
print(f'  RNFL 4분면(S/T/I/N) 가중평균 exact: {ex_sum/n_sum:.1f}% (n={n_sum})')
print(f'  RNFL inferior(rnfl_q_i) MAE: {ra["rnfl_q_i"]["mae"]}um (n={ra["rnfl_q_i"]["n_valid"]})')
ch_ex = ch_n = ch_mae = 0
for h in range(1, 13):
    k = f'rnfl_h{h:02d}'
    ch_ex += ra[k]['exact_match_pct'] * ra[k]['n_valid']
    ch_mae += ra[k]['mae'] * ra[k]['n_valid']
    ch_n += ra[k]['n_valid']
print(f'  clock-hour 가중평균 exact: {ch_ex/ch_n:.1f}%  MAE: {ch_mae/ch_n:.3f}um (n={ch_n})')

print('\n' + '=' * 90)
print('7. 본문값')
print('=' * 90)
print(f'  187안 MS: mean={np.mean(ms):.2f} +/- SD={np.std(ms,ddof=1):.2f} dB')
