# -*- coding: utf-8 -*-
"""187 코호트에서 topographic diag vs off-diag 평균 ρ (phase5_2 로직 재사용)."""
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
print('187 필터 후:', len(rows))
pids = [r['_pid'] for r in rows]

# RNFL 4분면 교차행렬
quads = list(p52.RNFL_QUAD_VF.keys())
diag_rho, off_rho = [], []
for q, (_, _) in p52.RNFL_QUAD_VF.items():
    xv = [r[f'_rnfl_q{q[-1]}'] for r in rows]
    for oq, (opts, _) in p52.RNFL_QUAD_VF.items():
        yv = [p52.vf_group_mean(r['_vf'], opts) for r in rows]
        blk = p52.corr_block(xv, yv, label='x', groups=pids if oq == q else None)
        rho = blk.get('spearman', {}).get('rho')
        if rho is None: continue
        (diag_rho if oq == q else off_rho).append(rho)
rnfl_diag = sum(diag_rho)/len(diag_rho)
rnfl_off = sum(off_rho)/len(off_rho)
print(f'RNFL: diag_mean_rho={rnfl_diag:.4f} offdiag={rnfl_off:.4f} delta={rnfl_diag-rnfl_off:+.4f}')

# GCIPL 6섹터 교차행렬
secs = list(p52.GCA_VF_DEF.keys())
diag_rho2, off_rho2 = [], []
for s in secs:
    xv = [r[f'_gcl_{s}'] for r in rows]
    for osec in secs:
        yv = [p52.vf_group_mean(r['_vf'], p52.GCA_VF_PTS[osec]) for r in rows]
        blk = p52.corr_block(xv, yv, label='x', groups=pids if osec == s else None)
        rho = blk.get('spearman', {}).get('rho')
        if rho is None: continue
        (diag_rho2 if osec == s else off_rho2).append(rho)
gcl_diag = sum(diag_rho2)/len(diag_rho2)
gcl_off = sum(off_rho2)/len(off_rho2)
print(f'GCIPL: diag_mean_rho={gcl_diag:.4f} offdiag={gcl_off:.4f} delta={gcl_diag-gcl_off:+.4f}')
