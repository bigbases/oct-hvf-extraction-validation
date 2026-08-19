# -*- coding: utf-8 -*-
"""
phase5_1 (전역 SF 상관) 과 phase5_4 (GCL/RNFL 결합모델) 를 187 코호트로 재생성.
_regen_fig2 와 동일 원칙: 원본 phase 스크립트 로직은 수정하지 않고 주입만.
- 187 필터: cohort_final_D.csv 의 (patient_id, eye, vf_date) 키.
- 주입: csv.DictReader 를 파일명 기준 shim 으로 교체 → analysis_master.csv 읽을 때만 187 필터.
- 백업은 호출 전 별도로 수행됨 (_backup_*_275_*.json).
"""
import csv as _csv, sys, os
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))

# 187 keys
D = set()
for row in _csv.reader(open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig')):
    if row == ['patient_id', 'eye', 'vf_date']:
        continue
    D.add(tuple(row))
assert len(D) == 187, len(D)

# csv.DictReader shim: analysis_master.csv 를 읽을 때만 187 필터
_orig = _csv.DictReader
class _FilteredDictReader:
    def __new__(cls, f, *a, **k):
        rdr = _orig(f, *a, **k)
        name = str(getattr(f, 'name', ''))
        if 'analysis_master' in name.replace('\\', '/'):
            rows = [r for r in rdr
                    if (r.get('patient_id'), r.get('eye'), r.get('vf_date')) in D]
            return iter(rows)
        return rdr
_csv.DictReader = _FilteredDictReader

# phase5_1 재생성
import phase5_1_correlations as p51
print('=== phase5_1 (187) 재생성 ===')
p51.main()

# phase5_4 재생성
import phase5_4_gcl_rnfl_compare as p54
print('\n=== phase5_4 (187) 재생성 ===')
p54.main()

# 복원
_csv.DictReader = _orig
print('\n재생성 완료 (results/phase5_1_correlations.json, results/phase5_4_gcl_rnfl_compare.json)')
