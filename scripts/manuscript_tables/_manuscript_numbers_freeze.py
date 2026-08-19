# -*- coding: utf-8 -*-
"""
원고 수치 전량 재계산 + freeze. 코호트 D 유일 출처 = cohort_final_D.csv(187).
vfacc=324 pool, octacc=280 GT pool. 결과를 _manuscript_numbers_freeze.txt 로 저장.
"""
import csv, sys, json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
from phase5_4_gcl_rnfl_compare import load_df, lme_single, lme_joint, complementarity_stats  # noqa: E402

OUT = []
def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s); OUT.append(s)

def to_f(v):
    try: return float(v)
    except (TypeError, ValueError): return None

# ── 코호트 D 유일 출처 ──
D_KEYS = set()
with open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig') as f:
    for row in csv.reader(f):
        if row == ['patient_id', 'eye', 'vf_date']: continue
        D_KEYS.add(tuple(row))
assert len(D_KEYS) == 187, len(D_KEYS)

am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
am_idx = {(r['patient_id'], r['eye'], r['vf_date']): r for r in am}
am_key_by_pideye = {(r['patient_id'], r['eye']): (r['patient_id'], r['eye'], r['vf_date']) for r in am}
d_rows = [am_idx[k] for k in D_KEYS]

rel = {(r['patient_id'], r['eye'], r['vf_date']): r
       for r in csv.DictReader(open(_ROOT / 'vf_reliability.csv', encoding='utf-8-sig'))}

# ============ 코호트 정의 함수들 (275 base = analysis_master) ============
def bad_flags(r):
    pid, eye, vfdate = r['patient_id'], r['eye'], r['vf_date']
    ss_gca  = to_f(r['ss_gca_od'])  if eye == 'OD' else to_f(r['ss_gca_os'])
    ss_rnfl = to_f(r['ss_rnfl_od']) if eye == 'OD' else to_f(r['ss_rnfl_os'])
    bad_ss = (ss_gca is not None and ss_gca < 6) or (ss_rnfl is not None and ss_rnfl < 6)
    rv = rel.get((pid, eye, vfdate))
    fl = fp = fn = None; low = None
    if rv is not None:
        num, den = to_f(rv['fix_loss_num']), to_f(rv['fix_loss_denom'])
        fl = 100*num/den if (num is not None and den) else None
        fp = to_f(rv['false_pos']); fn = to_f(rv['false_neg']); low = rv['low_reliability']
    return dict(fl=fl, fp=fp, fn=fn, low=low, bad_ss=bad_ss,
                sev=r.get('severity_tertile', ''))

def in_cohort(r, mode):
    b = bad_flags(r)
    fp_bad = b['fp'] is not None and b['fp'] > 15
    fn_bad = b['fn'] is not None and b['fn'] > 20
    fl_bad = b['fl'] is not None and b['fl'] > 33
    low_bad = b['low'] == 'Y'
    if mode == 'A':       # 무필터
        return True
    if mode == 'B':       # 표준(instrument-flag 계열): FL33/FP15/FN20/배너/ss
        return not (fl_bad or fp_bad or fn_bad or low_bad or b['bad_ss'])
    if mode == 'D':       # evidence-based: FP15/FN20/ss
        return not (fp_bad or fn_bad or b['bad_ss'])
    if mode == 'R':       # relaxed FN (severity-specific, Yohannan): mild/mod FN<25, severe FN<50
        sev = b['sev']
        if sev == 'severe':
            fn_bad_r = b['fn'] is not None and b['fn'] > 50
        else:
            fn_bad_r = b['fn'] is not None and b['fn'] > 25
        return not (fp_bad or fn_bad_r or b['bad_ss'])
    raise ValueError(mode)

df_full = load_df()
df_full['keyv'] = df_full.apply(lambda r: am_key_by_pideye.get((r['pid'], r['eye'])), axis=1)

def df_for(mode):
    keep = {am_key_by_pideye[(r['patient_id'], r['eye'])] for r in am if in_cohort(r, mode)}
    return df_full[df_full['keyv'].isin(keep)].copy(), len({(r['patient_id'],r['eye']) for r in am if in_cohort(r, mode)})

# ============ tab:cohortcurve ============
P('='*95)
P('tab:cohortcurve — 4개 코호트, ΔAIC(vsGCL, N에 비례) vs rnfl_q_i joint β(N-비교가능)')
P('='*95)
P(f'  {"코호트":32s}  {"N":5s}  {"ΔAIC(vsGCL)":12s}  {"ΔAIC(vsRNFL)":12s}  {"rnfl_q_i joint β":16s}  {"gcl_avg joint β":16s}')
curve = {}
for mode, label in [('A','무필터'), ('B','instrument flag(FL33/FP/FN/배너/ss)'),
                    ('D','evidence-based(FP15/FN20/ss)'), ('R','relaxed FN(severity-spec)')]:
    d, n = df_for(mode)
    j = lme_joint(d, 'gcl_avg', 'rnfl_q_i')
    diag = complementarity_stats(d, 'gcl_avg', 'rnfl_q_i')
    curve[mode] = dict(n=n, daic_gcl=diag.get('delta_aic_vs_gcl_only'),
                       daic_rnfl=diag.get('delta_aic_vs_rnfl_only'),
                       rnfl_b=j['rnfl_q_i']['beta'], gcl_b=j['gcl_avg']['beta'])
    P(f'  {label:32s}  {n:5d}  {diag.get("delta_aic_vs_gcl_only"):12.2f}  '
      f'{diag.get("delta_aic_vs_rnfl_only"):12.2f}  {j["rnfl_q_i"]["beta"]:16.4f}  {j["gcl_avg"]["beta"]:16.4f}')
P('  ※ ΔAIC는 표본크기 N에 비례해 커짐 → 서로 다른 N 코호트 간 "강도" 비교 부적절.')
P('  ※ rnfl_q_i joint β(표준화 계수)는 N-비교가능 → 이걸로 보면 D가 최대여야 서사 성립.')

# ============ df_D 고정 ============
df_D = df_for('D')[0]
assert len(df_D) == 187

# ============ tab:sf (187) ============
P('\n'+'='*95); P('tab:sf (n=187, cohort_final_D.csv)'); P('='*95)
struct_cols = ['gcl_avg', 'gcl_min', 'rnfl_q_i', 'rnfl_quad_mean', 'rnfl_ch_mean']
sf = {}
for col in struct_cols:
    r = lme_single(df_D, col); sf[col] = r
    P(f'  {col:16s}  n={r["n"]:3d}  rho={r["spearman"]["rho"]:.4f}  beta={r["lme"]["beta"]:+.4f}  '
      f'CI=[{r["lme"]["ci"][0]:+.4f},{r["lme"]["ci"][1]:+.4f}]  p={r["lme"]["p"]:.5f}')
# 순위
ranked = sorted(struct_cols, key=lambda c: -sf[c]['lme']['beta'])
P('  β 순위: ' + ' > '.join(ranked))

# ============ tab:comp (187) ============
P('\n'+'='*95); P('tab:comp (n=187)'); P('='*95)
j = lme_joint(df_D, 'gcl_avg', 'rnfl_q_i')
diag = complementarity_stats(df_D, 'gcl_avg', 'rnfl_q_i')
g = j['gcl_avg']; rr = j['rnfl_q_i']
solo_g = sf['gcl_avg']['lme']['beta']; solo_r = sf['rnfl_q_i']['lme']['beta']
P(f'  gcl_avg joint β={g["beta"]:.4f}(p={g["p"]:.5f})  rnfl_q_i joint β={rr["beta"]:.4f}(p={rr["p"]:.5f})')
P(f'  VIF={diag.get("vif")}  predictor_corr_r={diag.get("predictor_corr_r")}  '
  f'ΔAIC(vsGCL)={diag.get("delta_aic_vs_gcl_only")}  ΔAIC(vsRNFL)={diag.get("delta_aic_vs_rnfl_only")}')
P(f'  감쇠율: gcl_avg {solo_g:.4f}->{g["beta"]:.4f} ({(solo_g-g["beta"])/solo_g*100:.1f}%)  '
  f'rnfl_q_i {solo_r:.4f}->{rr["beta"]:.4f} ({(solo_r-rr["beta"])/solo_r*100:.1f}%)')

# ============ tab:autovsgt (187) ============
P('\n'+'='*95); P('tab:autovsgt (n=187, y=vf_mean GT고정)'); P('='*95)
octv = {(r['patient_id'], r['eye'], r['oct_date']): r
        for r in csv.DictReader(open(_ROOT / 'data' / 'oct_values_raw.csv', encoding='utf-8-sig'))}
rnfl_raw = {(r['patient_id'], r['eye'], r['oct_date']): r
            for r in csv.DictReader(open(_ROOT / 'data' / 'rnfl_detail_raw.csv', encoding='utf-8-sig'))}
RNFL_TN_SWAP = {'rnfl_q_t': 'rnfl_q_n', 'rnfl_q_n': 'rnfl_q_t'}
def auto_rnfl_q_col(q, eye):
    col = f'rnfl_q_{q}'
    return RNFL_TN_SWAP[col] if (eye == 'OS' and col in RNFL_TN_SWAP) else col
arec = []
for k in D_KEYS:
    r = am_idx[k]; pid, eye, octd = r['patient_id'], r['eye'], r['oct_date']
    ov = octv.get((pid,eye,octd)); rv = rnfl_raw.get((pid,eye,octd))
    if ov is None or rv is None: continue
    e = eye.lower()
    qi=to_f(rv.get(auto_rnfl_q_col('i',eye))); qs=to_f(rv.get(auto_rnfl_q_col('s',eye)))
    qt=to_f(rv.get(auto_rnfl_q_col('t',eye))); qn=to_f(rv.get(auto_rnfl_q_col('n',eye)))
    qv=[v for v in [qs,qi,qt,qn] if v is not None]
    ch=[to_f(rv.get(f'rnfl_h{h:02d}')) for h in range(1,13)]; chok=[v for v in ch if v is not None]
    arec.append({'pid':pid,'eye':eye,'vf_mean':to_f(r.get('vf_mean')),
                 'gcl_avg':to_f(ov.get(f'avg_gcl_{e}')),'gcl_min':to_f(ov.get(f'min_gcl_{e}')),
                 'rnfl_q_i':qi,'rnfl_quad_mean':sum(qv)/len(qv) if qv else None,
                 'rnfl_ch_mean':sum(chok)/len(chok) if chok else None})
df_auto = pd.DataFrame(arec)
for col in struct_cols:
    b_gt=sf[col]['lme']['beta']; ci_gt=sf[col]['lme']['ci']
    ra=lme_single(df_auto,col); b_a=ra['lme']['beta']; ci_a=ra['lme']['ci']
    ov_=ci_gt[0]<ci_a[1] and ci_a[0]<ci_gt[1]
    P(f'  {col:16s}  GT_β={b_gt:+.4f}  auto_β={b_a:+.4f}  Δβ={b_a-b_gt:+.4f}  CI겹침={"yes" if ov_ else "no"}')

# ============ tab:cohort (187) ============
P('\n'+'='*95); P('tab:cohort (n=187)'); P('='*95)
npat=len({r['patient_id'] for r in d_rows})
# bilateral/unilateral
pc=Counter(r['patient_id'] for r in d_rows)
bilat=sum(1 for v in pc.values() if v==2); uni=sum(1 for v in pc.values() if v==1)
P(f'  eyes/patients: 187/{npat} ({bilat} bilateral, {uni} unilateral)')
patients={}
for r in d_rows: patients.setdefault(r['patient_id'],[]).append(r)
sc=Counter()
for pid,rows in patients.items():
    sx={r['sex'] for r in rows if r['sex']}
    sc[list(sx)[0] if len(sx)==1 else ('NR' if not sx else 'conflict')]+=1
P(f'  sex(patient): Male={sc.get("Male",0)} Female={sc.get("Female",0)} NR={sc.get("NR",0)}')
ages=[to_f(r['age_at_oct']) for r in d_rows if to_f(r['age_at_oct']) is not None]
P(f'  age: {np.mean(ages):.1f}±{np.std(ages,ddof=1):.1f} range[{min(ages):.1f},{max(ages):.1f}] n={len(ages)}')
ms=[to_f(r['vf_mean']) for r in d_rows if to_f(r['vf_mean']) is not None]
P(f'  MS: {np.mean(ms):.2f}±{np.std(ms,ddof=1):.2f} range[{min(ms):.1f},{max(ms):.1f}] n={len(ms)}')
sev=Counter(r['severity_tertile'] for r in d_rows if r['severity_tertile'])
P(f'  severity: mild={sev.get("mild",0)} moderate={sev.get("moderate",0)} severe={sev.get("severe",0)}')
sita={(r['patient_id'],r['eye'],r['vf_date']):r['sita_strategy'] for r in csv.DictReader(open(_ROOT/'data'/'sita_per_eye.csv',encoding='utf-8-sig'))}
st=Counter(sita.get((r['patient_id'],r['eye'],r['vf_date']),'NR') for r in d_rows)
P(f'  VF strategy: {dict(st)}')

# ============ tab:vfacc (324) & tab:octacc (280) ============
r3=json.load(open(_ROOT/'results'/'phase3_extraction_accuracy.json',encoding='utf-8'))
va=r3['vf_accuracy']['accuracy_normal_eyes']
P('\n'+'='*95); P('tab:vfacc (n=324 pool)'); P('='*95)
P(f'  exact={va["exact_match_pct"]} ±1dB={va["within1dB_pct"]} ±2dB={va["within2dB_pct"]} '
  f'MAE={va["mae"]} RMSE={va["rmse"]} missing={va["missing_rate_raw_pct"]} nvalid={va["n_valid"]}/{va["n_pairs"]}')

P('\n'+'='*95); P('tab:octacc (n=280 GT pool)'); P('='*95)
ga=r3['gcl_accuracy']['per_column']; ra=r3['rnfl_accuracy']['per_column']
for k in ['avg_gcl_od','avg_gcl_os','min_gcl_od','min_gcl_os']:
    s=ga[k]; P(f'  {k}: exact={s["exact_match_pct"]} MAE={s["mae"]} miss={s["missing_rate_raw_pct"]} n={s["n_valid"]}')
# GCIPL 상/하 그룹 pooled
sup=['od_s_sup','od_s_sup_t','od_s_sup_n','os_s_sup','os_s_sup_t','os_s_sup_n']
inf=['od_s_inf','od_s_inf_t','od_s_inf_n','os_s_inf','os_s_inf_t','os_s_inf_n']
def pooled(cols):
    ex=n=0
    exs=[]
    for c in cols:
        s=ga[c]; ex+=s['exact_match_pct']*s['n_valid']; n+=s['n_valid']; exs.append(s['exact_match_pct'])
    return ex/n, min(exs), max(exs)
sp=pooled(sup); ip=pooled(inf)
P(f'  GCIPL superior group: pooled_exact={sp[0]:.1f} range[{sp[1]:.1f},{sp[2]:.1f}]')
P(f'  GCIPL inferior group: pooled_exact={ip[0]:.1f} range[{ip[1]:.1f},{ip[2]:.1f}]')
qex=qn=0
for k in ['rnfl_q_s','rnfl_q_t','rnfl_q_i','rnfl_q_n']:
    qex+=ra[k]['exact_match_pct']*ra[k]['n_valid']; qn+=ra[k]['n_valid']
P(f'  RNFL 4quad pooled exact={qex/qn:.1f} (n={qn})')
P(f'  RNFL each: S={ra["rnfl_q_s"]["exact_match_pct"]} T={ra["rnfl_q_t"]["exact_match_pct"]} '
  f'I={ra["rnfl_q_i"]["exact_match_pct"]} N={ra["rnfl_q_n"]["exact_match_pct"]}')
P(f'  RNFL inferior(q_i): exact={ra["rnfl_q_i"]["exact_match_pct"]} MAE={ra["rnfl_q_i"]["mae"]} n={ra["rnfl_q_i"]["n_valid"]}')
chex=chn=chmae=0
for h in range(1,13):
    k=f'rnfl_h{h:02d}'; chex+=ra[k]['exact_match_pct']*ra[k]['n_valid']; chn+=ra[k]['n_valid']; chmae+=ra[k]['mae']*ra[k]['n_valid']
chvals=[ra[f'rnfl_h{h:02d}']['exact_match_pct'] for h in range(1,13)]
P(f'  clock-hour: pooled exact={chex/chn:.1f} MAE={chmae/chn:.3f} range[{min(chvals):.1f},{max(chvals):.1f}] n={chn}')

with open(_ROOT/'_manuscript_numbers_freeze.txt','w',encoding='utf-8') as f:
    f.write('\n'.join(OUT))
P('\n[저장] _manuscript_numbers_freeze.txt')

# tab:cohortcurve 구조화 로그 (재현용, 비 PHI 집계)
_labels = {'A': 'No reliability filter', 'B': 'Instrument flag (FL33/FP15/FN20/banner/ss)',
           'D': 'Evidence-based (FP>15, FN>20, ss)', 'R': 'Relaxed FN (severity-specific)'}
_cc = {'generated': __import__('datetime').date.today().isoformat(),
       'note': 'Standardized joint-model beta of rnfl_q_i (adjusted for gcl_avg); comparable across cohort sizes. '
               'DeltaAIC scales with N and is not used for cross-cohort strength comparison.',
       'rows': {m: {'label': _labels[m], 'n_eyes': curve[m]['n'],
                    'rnfl_q_i_joint_beta': round(curve[m]['rnfl_b'], 4),
                    'gcl_avg_joint_beta': round(curve[m]['gcl_b'], 4),
                    'delta_aic_vs_gcl': curve[m]['daic_gcl'], 'delta_aic_vs_rnfl': curve[m]['daic_rnfl']}
                for m in ['A', 'B', 'D', 'R']}}
json.dump(_cc, open(_ROOT/'results'/'cohortcurve.json', 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
P('[저장] results/cohortcurve.json')
