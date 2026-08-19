# -*- coding: utf-8 -*-
"""
Fig 2 (sf_heatmap.pdf) 재설계: RNFL 4x4 + GCIPL 6x6 를 둘 다 full cross-correlation
matrix 로 그린다(apples-to-apples). 데이터는 두 지표 모두 topographic(대각) 특이성이
없음을 보여줌 — 대각을 최댓값으로 강조하지 않고, 해부학적 매칭 셀만 외곽선으로 표시.
색 스케일은 데이터 실범위에 맞춰 좁힘(0.10-0.55)으로써 RNFL 행효과(상/하 강, 이/비측 약)
와 GCIPL 균일 패턴이 시각적으로 드러나게 함.
187 코호트(cohort_final_D.csv) 기준으로 이 스크립트가 직접 재계산.
출력: paper/manuscript/figures/sf_heatmap.pdf (+ .png)
"""
import csv, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'scripts'))
import phase5_2_sector_vf as p52

D_KEYS = set()
for row in csv.reader(open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig')):
    if row == ['patient_id', 'eye', 'vf_date']: continue
    D_KEYS.add(tuple(row))

rows = p52.load_data()
rows = [r for r in rows if (r['patient_id'], r['eye'], r['vf_date']) in D_KEYS]
assert len(rows) == 187

def rho(xv, yv):
    return p52.corr_block(xv, yv, 'x').get('spearman', {}).get('rho', float('nan'))

# ── RNFL 4x4 ──
quads = list(p52.RNFL_QUAD_VF.keys())  # q_s, q_i, q_t, q_n
RNFL_ROWLAB = {'q_s': 'S', 'q_i': 'I', 'q_t': 'T', 'q_n': 'N'}
RNFL_COLLAB = {'q_s': 'inf VF', 'q_i': 'sup VF', 'q_t': 'central', 'q_n': 'periph-T'}
rnfl = [[rho([r[f'_rnfl_q{rq[-1]}'] for r in rows],
             [p52.vf_group_mean(r['_vf'], p52.RNFL_QUAD_VF[cq][0]) for r in rows])
         for cq in quads] for rq in quads]

# ── GCIPL 6x6 ──
secs = list(p52.GCA_VF_DEF.keys())  # sup_t, sup, sup_n, inf_n, inf, inf_t
GCL_LAB = {'sup_t': 'ST', 'sup': 'S', 'sup_n': 'SN', 'inf_n': 'IN', 'inf': 'I', 'inf_t': 'IT'}
gcl = [[rho([r[f'_gcl_{rs}'] for r in rows],
            [p52.vf_group_mean(r['_vf'], p52.GCA_VF_PTS[cs]) for r in rows])
        for cs in secs] for rs in secs]

# ── plot: raw 절대크기 (viridis 공유 스케일). residual 패널은 순열검정에서
# 새 정보가 확인되지 않아 제외(대각평균=raw matched-minus-unmatched 의 상수배,
# GCIPL 블록대칭은 double-centering 강제, RNFL 교차는 라벨공간 4개로 검정력 부족).
vmin, vmax = 0.10, 0.56
cmap = plt.get_cmap('viridis')
norm = Normalize(vmin=vmin, vmax=vmax)

def _tcolor(v):
    # 글자색은 '반올림된 표시값'의 셀 밝기로 결정 → 같은 표시값=같은 글자색(경계 불일치 제거).
    rr, gg, bb, _ = cmap(norm(round(v, 2)))
    return 'black' if (0.299*rr + 0.587*gg + 0.114*bb) > 0.6 else 'white'

def draw(ax, mat, rowlabs, collab, title):
    ax.imshow(mat, cmap=cmap, norm=norm, aspect='equal')
    n = len(mat)
    for i in range(n):
        for j in range(n):
            v = mat[i][j]
            ax.text(j, i, f'{v:.2f}', ha='center', va='center', color=_tcolor(v), fontsize=9)
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False,
                               edgecolor='0.15', lw=1.6))  # 대각(해부학 매칭) 위치 표시
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(collab, fontsize=9); ax.set_yticklabels(rowlabs, fontsize=9)
    ax.set_xlabel('VF region', fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.tick_params(length=0)

fig = plt.figure(figsize=(10.5, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[4.2, 6, 0.28], wspace=0.32)

axr = fig.add_subplot(gs[0, 0])
draw(axr, rnfl, [RNFL_ROWLAB[q] for q in quads], [RNFL_COLLAB[q] for q in quads],
     'RNFL quadrant $\\times$ VF region (4$\\times$4)')
axr.set_ylabel('RNFL quadrant (disc)', fontsize=9)

axg = fig.add_subplot(gs[0, 1])
# x축(열)은 GCIPL 섹터명이 아니라 그 섹터에 Garway-Heath 매칭된 시야 영역이므로
# 'field' 접미사로 y축(행, 섹터 자체)과 구분 — 180도 반전 미적용처럼 보이는 것 방지.
draw(axg, gcl, [GCL_LAB[s] for s in secs], [GCL_LAB[s] + ' field' for s in secs],
     'GCIPL sector $\\times$ VF region (6$\\times$6)')
axg.set_ylabel('GCIPL sector (macula)', fontsize=9)

cax = fig.add_subplot(gs[0, 2])
fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, label='Spearman $\\rho$',
             ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.55])

# 그림 내부 제목 없음 — LaTeX 캡션이 제목 역할(중복 제거). 패널별 ax.set_title(RNFL/GCIPL
# 구분)은 캡션이 "Left:/Right:"로 참조하는 식별자라 유지.
fig.tight_layout()
out_dir = _ROOT / 'paper' / 'manuscript' / 'figures'
for ext in ('pdf', 'png'):
    fig.savefig(out_dir / f'sf_heatmap.{ext}', bbox_inches='tight', dpi=300)
print('Fig 2 저장:', out_dir / 'sf_heatmap.pdf')

# 진단 출력
diag_r = [rnfl[i][i] for i in range(4)]; off_r = [rnfl[i][j] for i in range(4) for j in range(4) if i != j]
diag_g = [gcl[i][i] for i in range(6)]; off_g = [gcl[i][j] for i in range(6) for j in range(6) if i != j]
print(f'RNFL diag_mean={sum(diag_r)/4:.4f} off_mean={sum(off_r)/12:.4f} Δ={sum(diag_r)/4-sum(off_r)/12:+.4f}')
print(f'GCIPL diag_mean={sum(diag_g)/6:.4f} off_mean={sum(off_g)/30:.4f} Δ={sum(diag_g)/6-sum(off_g)/30:+.4f}')
