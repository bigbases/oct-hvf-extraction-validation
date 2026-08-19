# -*- coding: utf-8 -*-
"""
Fig 2 (sf_heatmap.pdf) 를 187 코호트(cohort_final_D.csv) 기준으로 재생성.
phase5_2 의 load_data 를 187 필터로 monkeypatch → 표준 계산 로직 그대로 재사용 →
results/phase5_2_sector_correlations.json 을 187 기준으로 갱신(276본은 백업) →
stage_30_figures.make_sf_heatmap 으로 PDF/PNG 재생성.
기존 스크립트 로직은 수정하지 않음(주입만).
"""
import csv, sys, shutil
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))

import phase5_2_sector_vf as p52

# 187 키 로드
D_KEYS = set()
for row in csv.reader(open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig')):
    if row == ['patient_id', 'eye', 'vf_date']: continue
    D_KEYS.add(tuple(row))
assert len(D_KEYS) == 187

# 276 JSON 백업
json_path = _ROOT / 'results' / 'phase5_2_sector_correlations.json'
if json_path.exists():
    bak = _ROOT / 'results' / '_backup_phase5_2_276_sector_correlations.json'
    if not bak.exists():
        shutil.copy(json_path, bak)
        print(f'백업: {bak.name}')

# load_data 를 187 필터로 교체
_orig_load = p52.load_data
def _load_187():
    rows = _orig_load()
    filtered = [r for r in rows if (r['patient_id'], r['eye'], r['vf_date']) in D_KEYS]
    print(f'[187 필터] {len(rows)} -> {len(filtered)} rows')
    return filtered
p52.load_data = _load_187

# phase5_2 전체 재실행 → json(187) 갱신
p52.main()

# 히트맵 재생성
import stage_30_figures as s30
out_dir = _ROOT / 'paper' / 'manuscript' / 'figures'
out_dir.mkdir(parents=True, exist_ok=True)
s30.make_sf_heatmap(out_dir)
print('Fig 2 재생성 완료:', out_dir / 'sf_heatmap.pdf')
