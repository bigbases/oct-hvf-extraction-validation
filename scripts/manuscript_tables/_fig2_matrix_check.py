# -*- coding: utf-8 -*-
"""187 코호트 RNFL 4x4 + GCIPL 6x6 full matrix 계산 + 대각이 행/열 최댓값인지 진단."""
import csv, sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))
import phase5_2_sector_vf as p52

D_KEYS = set()
for row in csv.reader(open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig')):
    if row == ['patient_id', 'eye', 'vf_date']: continue
    D_KEYS.add(tuple(row))

rows = p52.load_data()
rows = [r for r in rows if (r['patient_id'], r['eye'], r['vf_date']) in D_KEYS]
print('n =', len(rows))
pids = [r['_pid'] for r in rows]

def rho(xv, yv, grp=None):
    b = p52.corr_block(xv, yv, 'x', groups=grp)
    return b.get('spearman', {}).get('rho')

# ── RNFL 4x4 ──
quads = list(p52.RNFL_QUAD_VF.keys())  # q_s, q_i, q_t, q_n
print('\n=== RNFL 4x4 (row=structure quadrant, col=VF region) ===')
print('        ' + ''.join(f'{q:>8s}' for q in quads))
rnfl_mat = {}
for rq in quads:
    xv = [r[f'_rnfl_q{rq[-1]}'] for r in rows]
    line = f'{rq:8s}'
    rnfl_mat[rq] = {}
    for cq in quads:
        _, opts = None, p52.RNFL_QUAD_VF[cq][0]
        yv = [p52.vf_group_mean(r['_vf'], opts) for r in rows]
        v = rho(xv, yv)
        rnfl_mat[rq][cq] = v
        mark = '*' if rq == cq else ' '
        line += f'{v:7.3f}{mark}'
    print(line)
# 대각이 행 최댓값인가
print('  대각이 행 최댓값?:', {rq: (max(rnfl_mat[rq], key=lambda c: rnfl_mat[rq][c]) == rq) for rq in quads})
print('  대각이 열 최댓값?:', {cq: (max(quads, key=lambda r: rnfl_mat[r][cq]) == cq) for cq in quads})

# ── GCIPL 6x6 ──
secs = list(p52.GCA_VF_DEF.keys())
print('\n=== GCIPL 6x6 (row=structure sector, col=VF region) ===')
print('        ' + ''.join(f'{s:>8s}' for s in secs))
gcl_mat = {}
for rs in secs:
    xv = [r[f'_gcl_{rs}'] for r in rows]
    line = f'{rs:8s}'
    gcl_mat[rs] = {}
    for cs in secs:
        yv = [p52.vf_group_mean(r['_vf'], p52.GCA_VF_PTS[cs]) for r in rows]
        v = rho(xv, yv)
        gcl_mat[rs][cs] = v
        mark = '*' if rs == cs else ' '
        line += f'{v:7.3f}{mark}'
    print(line)
print('  대각이 행 최댓값?:', {rs: (max(gcl_mat[rs], key=lambda c: gcl_mat[rs][c]) == rs) for rs in secs})
n_row_max_rnfl = sum(1 for rq in quads if max(rnfl_mat[rq], key=lambda c: rnfl_mat[rq][c]) == rq)
n_row_max_gcl = sum(1 for rs in secs if max(gcl_mat[rs], key=lambda c: gcl_mat[rs][c]) == rs)
print(f'\n요약: RNFL 대각=행최댓값 {n_row_max_rnfl}/4,  GCIPL 대각=행최댓값 {n_row_max_gcl}/6')
