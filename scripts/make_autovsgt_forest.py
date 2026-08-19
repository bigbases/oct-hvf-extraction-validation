# -*- coding: utf-8 -*-
"""
Fig 4: 자동추출 vs 수동검수(GT) 구조-기능 계수 비교 (187 코호트), 2패널.
좌: 두 방법의 β + marginal 95% CI (겹침 시각화).
우: Δβ(auto-ref)의 CI를 "환자-클러스터 paired bootstrap"으로 직접 추정 —
    같은 187안을 두 방법으로 잰 paired 데이터라 두 추정치가 강상관이므로 차이의 CI는
    marginal CI보다 훨씬 좁다. 0을 포함/근접하면 equivalence(추출오차가 결론 안 바꿈).
187 코호트(cohort_final_D.csv) 기준으로 직접 재계산(재현 가능).
출력: paper/manuscript/figures/autovsgt_forest.pdf (+ .png)
"""
import csv, sys, random, json
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.formula.api as smf

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'scripts'))
from phase5_4_gcl_rnfl_compare import load_df, lme_single  # noqa: E402

B = 400           # bootstrap 반복 (LME, per fit ~56ms)
random.seed(42); np.random.seed(42)


def to_f(v):
    try: return float(v)
    except (TypeError, ValueError): return None


# ── 187 코호트 ──
D_KEYS = set()
for row in csv.reader(open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig')):
    if row == ['patient_id', 'eye', 'vf_date']: continue
    D_KEYS.add(tuple(row))
am = list(csv.DictReader(open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
am_idx = {(r['patient_id'], r['eye'], r['vf_date']): r for r in am}
am_key = {(r['patient_id'], r['eye']): (r['patient_id'], r['eye'], r['vf_date']) for r in am}

df = load_df()
df['keyv'] = df.apply(lambda r: am_key.get((r['pid'], r['eye'])), axis=1)
df_gt = df[df['keyv'].isin(D_KEYS)].copy().reset_index(drop=True)

octv = {(r['patient_id'], r['eye'], r['oct_date']): r
        for r in csv.DictReader(open(_ROOT / 'data' / 'oct_values_raw.csv', encoding='utf-8-sig'))}
rnfl_raw = {(r['patient_id'], r['eye'], r['oct_date']): r
            for r in csv.DictReader(open(_ROOT / 'data' / 'rnfl_detail_raw.csv', encoding='utf-8-sig'))}
RNFL_TN_SWAP = {'rnfl_q_t': 'rnfl_q_n', 'rnfl_q_n': 'rnfl_q_t'}
def auto_q(q, eye):
    c = f'rnfl_q_{q}'
    return RNFL_TN_SWAP[c] if (eye == 'OS' and c in RNFL_TN_SWAP) else c

arec = []
for k in D_KEYS:
    r = am_idx[k]; pid, eye, octd = r['patient_id'], r['eye'], r['oct_date']
    ov = octv.get((pid, eye, octd)); rv = rnfl_raw.get((pid, eye, octd))
    if ov is None or rv is None: continue
    e = eye.lower()
    qi = to_f(rv.get(auto_q('i', eye))); qs = to_f(rv.get(auto_q('s', eye)))
    qt = to_f(rv.get(auto_q('t', eye))); qn = to_f(rv.get(auto_q('n', eye)))
    qv = [v for v in [qs, qi, qt, qn] if v is not None]
    ch = [to_f(rv.get(f'rnfl_h{h:02d}')) for h in range(1, 13)]; chok = [v for v in ch if v is not None]
    arec.append({'pid': pid, 'eye': eye, 'vf_mean': to_f(r.get('vf_mean')),
                 'gcl_avg': to_f(ov.get(f'avg_gcl_{e}')), 'gcl_min': to_f(ov.get(f'min_gcl_{e}')),
                 'rnfl_q_i': qi, 'rnfl_quad_mean': sum(qv)/len(qv) if qv else None,
                 'rnfl_ch_mean': sum(chok)/len(chok) if chok else None})
df_auto = pd.DataFrame(arec)

METRICS = [('rnfl_q_i', 'RNFL inferior quadrant'),
           ('gcl_avg', 'GCIPL average'),
           ('rnfl_quad_mean', 'RNFL mean quadrant'),
           ('rnfl_ch_mean', 'RNFL mean clock-hour'),
           ('gcl_min', 'GCIPL minimum')]


def lme_beta(s, col):
    """표준화 LME 고정효과 β (point/bootstrap 공통 — lme_single과 동일 정의)."""
    s = s[[col, 'vf_mean', 'pid']].dropna().copy()
    if len(s) < 10 or s[col].std() == 0:
        return np.nan
    s['zx'] = (s[col] - s[col].mean()) / s[col].std()
    s['zy'] = (s['vf_mean'] - s['vf_mean'].mean()) / s['vf_mean'].std()
    try:
        m = smf.mixedlm('zy ~ zx', s, groups=s['pid']).fit(reml=True, disp=False)
        return float(m.params['zx'])
    except Exception:
        return np.nan


# ── point estimates (LME, 원고 표와 동일) ──
rows = []
for col, label in METRICS:
    g = lme_single(df_gt, col); a = lme_single(df_auto, col)
    rows.append(dict(col=col, label=label,
                     g_b=g['lme']['beta'], g_lo=g['lme']['ci'][0], g_hi=g['lme']['ci'][1],
                     a_b=a['lme']['beta'], a_lo=a['lme']['ci'][0], a_hi=a['lme']['ci'][1]))

# ── 환자-클러스터 paired bootstrap → Δβ CI (LME, 점추정과 동일 방법) ──
# 주의: 복원추출로 같은 환자가 2번 뽑히면 mixedlm(groups=pid)이 둘을 한 그룹으로
# 병합해 cluster bootstrap이 깨진다. 뽑을 때마다 고유 그룹 id(pid__draw)를 부여해
# 중복 환자를 독립 클러스터로 취급.
# 결과 캐시: cosmetic 재실행 시 4분 부트스트랩을 다시 안 돌리게 JSON에 저장/로드.
# 데이터/코드가 바뀌면 캐시 파일을 지우고 다시 실행할 것.
CACHE = _ROOT / 'results' / '_fig4_boot_cache.json'
if CACHE.exists():
    cached = json.load(open(CACHE, encoding='utf-8'))
    for r in rows:
        c = cached[r['col']]
        r['d_point'] = r['a_b'] - r['g_b']; r['d_lo'] = c['lo']; r['d_hi'] = c['hi']; r['d_nboot'] = c['n']
    print(f'[cache] {CACHE.name} 로드 (부트스트랩 건너뜀)')
else:
    patients = sorted(df_gt['pid'].unique())
    gt_by_pat = defaultdict(list); auto_by_pat = defaultdict(list)
    for _, r in df_gt.iterrows():   gt_by_pat[r['pid']].append(r.to_dict())
    for _, r in df_auto.iterrows(): auto_by_pat[r['pid']].append(r.to_dict())
    boot = {col: [] for col, _ in METRICS}
    for b in range(B):
        samp = [random.choice(patients) for _ in patients]
        grows, arows = [], []
        for di, p in enumerate(samp):
            gid = f'{p}__{di}'
            for x in gt_by_pat[p]:   xx = dict(x); xx['pid'] = gid; grows.append(xx)
            for x in auto_by_pat[p]: xx = dict(x); xx['pid'] = gid; arows.append(xx)
        gsub = pd.DataFrame(grows); asub = pd.DataFrame(arows)
        for col, _ in METRICS:
            bg = lme_beta(gsub, col); ba = lme_beta(asub, col)
            if not (np.isnan(bg) or np.isnan(ba)):
                boot[col].append(ba - bg)
    cache_out = {}
    for r in rows:
        arr = np.array(boot[r['col']])
        r['d_point'] = r['a_b'] - r['g_b']
        r['d_lo'] = float(np.percentile(arr, 2.5)); r['d_hi'] = float(np.percentile(arr, 97.5))
        r['d_nboot'] = len(arr)
        cache_out[r['col']] = {'lo': r['d_lo'], 'hi': r['d_hi'], 'n': r['d_nboot']}
    json.dump(cache_out, open(CACHE, 'w', encoding='utf-8'), indent=2)
    print(f'[cache] {CACHE.name} 저장')

# ── 점추정 로그 출력 저장 (tab:autovsgt 재현용; 비 PHI 집계) ──
COEF_OUT = _ROOT / 'results' / 'autovsgt_coefficients.json'
coef = {'generated': __import__('datetime').date.today().isoformat(),
        'cohort': '187 analytic (cohort_final_D.csv)',
        'n_gt': int(len(df_gt)), 'n_auto': int(len(df_auto)),
        'note': 'reference = manual-review (analysis_master.csv); automated = oct_values_raw.csv + rnfl_detail_raw.csv. '
                'Delta-beta CI = patient-clustered paired bootstrap (B={}).'.format(B),
        'metrics': {r['col']: {
            'reference_beta': round(r['g_b'], 4), 'reference_ci': [round(r['g_lo'], 4), round(r['g_hi'], 4)],
            'automated_beta': round(r['a_b'], 4), 'automated_ci': [round(r['a_lo'], 4), round(r['a_hi'], 4)],
            'delta_beta': round(r['d_point'], 4), 'delta_boot_ci': [round(r['d_lo'], 4), round(r['d_hi'], 4)],
            'delta_boot_ci_includes_zero': bool(r['d_lo'] <= 0 <= r['d_hi']), 'n_boot': r['d_nboot'],
        } for r in rows}}
json.dump(coef, open(COEF_OUT, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
print(f'[coef] {COEF_OUT.name} 저장 (tab:autovsgt 점추정)')

# ── plot ──
fig = plt.figure(figsize=(10.5, 5.0))
gs = fig.add_gridspec(1, 2, width_ratios=[2.6, 1.25], wspace=0.26)
C_GT = '#1f4e79'; C_AUTO = '#c55a11'
n = len(rows); ypos = list(range(n, 0, -1)); off = 0.17
TOP = n + 1.25          # 상단 헤드룸(범례/라벨용 빈 밴드)

# 좌: β forest
axL = fig.add_subplot(gs[0, 0])
for y, r in zip(ypos, rows):
    axL.plot([r['g_lo'], r['g_hi']], [y + off, y + off], color=C_GT, lw=2.4, solid_capstyle='round')
    axL.plot(r['g_b'], y + off, 'o', color=C_GT, ms=7, zorder=3)
    axL.plot([r['a_lo'], r['a_hi']], [y - off, y - off], color=C_AUTO, lw=2.4, solid_capstyle='round')
    axL.plot(r['a_b'], y - off, 's', color=C_AUTO, ms=6.5, zorder=3)
axL.set_yticks(ypos); axL.set_yticklabels([r['label'] for r in rows], fontsize=10)
axL.set_ylim(0.4, TOP); axL.set_xlim(0.36, 0.80)
axL.set_xlabel('Standardized $\\beta$ (95% CI, marginal)', fontsize=10)
axL.set_title('(a) Coefficient by extraction source', fontsize=10.5, pad=6)
for sp in ('top', 'right'): axL.spines[sp].set_visible(False)
axL.grid(axis='x', ls=':', alpha=0.4)
# 범례: 데이터 위 헤드룸 밴드에 가로 배치 (막대와 겹치지 않음)
legend = [Line2D([0], [0], color=C_GT, marker='o', lw=2.4, ms=7, label='Manual reference'),
          Line2D([0], [0], color=C_AUTO, marker='s', lw=2.4, ms=6.5, label='Automated')]
axL.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.5, 1.0),
           ncol=2, fontsize=9, frameon=False, handletextpad=0.4, columnspacing=1.6)

# 우: Δβ (paired bootstrap CI) forest, 0선
axR = fig.add_subplot(gs[0, 1])
axR.axvline(0, color='0.55', lw=1)
for y, r in zip(ypos, rows):
    axR.plot([r['d_lo'], r['d_hi']], [y, y], color='#4a4a4a', lw=2.4, solid_capstyle='round')
    axR.plot(r['d_point'], y, 'D', color='#4a4a4a', ms=6, zorder=3)
    # Δβ 수치는 막대 '위'에 (오른쪽 여백 충돌 회피)
    axR.text(r['d_point'], y + 0.34, f"{r['d_point']:+.3f}", va='center', ha='center',
             fontsize=8, color='black')
axR.set_yticks(ypos); axR.set_yticklabels([])
axR.set_ylim(0.4, TOP); axR.set_xlim(-0.17, 0.21)
axR.set_xticks([-0.1, 0, 0.1])
axR.set_xlabel('$\\Delta\\beta$ (auto $-$ ref)\npaired-bootstrap 95% CI', fontsize=9.5)
axR.set_title('(b) Paired difference', fontsize=10.5, pad=6)
for sp in ('top', 'right'): axR.spines[sp].set_visible(False)
axR.grid(axis='x', ls=':', alpha=0.4)

# 그림 내부 제목 없음 — LaTeX 캡션이 제목 역할(중복 제거). (a)/(b) 패널 라벨은 캡션이
# 직접 참조하는 식별자라 유지.
out_dir = _ROOT / 'paper' / 'manuscript' / 'figures'
for ext in ('pdf', 'png'):
    fig.savefig(out_dir / f'autovsgt_forest.{ext}', dpi=300, bbox_inches='tight')
print('Fig 4 저장:', out_dir / 'autovsgt_forest.pdf')
print(f'\n[freeze] B={B} 환자-클러스터 bootstrap')
for r in rows:
    contains0 = r['d_lo'] <= 0 <= r['d_hi']
    print(f"  {r['label']:24s} β_ref={r['g_b']:.3f} β_auto={r['a_b']:.3f}  "
          f"Δβ={r['d_point']:+.3f} [{r['d_lo']:+.3f},{r['d_hi']:+.3f}]  0포함={contains0}")
