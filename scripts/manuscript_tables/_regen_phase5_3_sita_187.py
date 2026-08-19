# -*- coding: utf-8 -*-
"""phase5_3 (covariate) + sita_sensitivity 를 187 코호트로 재생성. _regen 방식(주입만)."""
import csv as _csv, sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / 'scripts'))

D = set()
for row in _csv.reader(open(_ROOT / 'cohort_final_D.csv', encoding='utf-8-sig')):
    if row == ['patient_id', 'eye', 'vf_date']:
        continue
    D.add(tuple(row))
assert len(D) == 187

_orig = _csv.DictReader
class _FDR:
    def __new__(cls, f, *a, **k):
        rdr = _orig(f, *a, **k)
        name = str(getattr(f, 'name', '')).replace('\\', '/')
        if 'analysis_master' in name:
            return iter([r for r in rdr
                         if (r.get('patient_id'), r.get('eye'), r.get('vf_date')) in D])
        return rdr
_csv.DictReader = _FDR

import phase5_3_covariate as p3
print("=== phase5_3 covariate (187) ===")
p3.main()

import phase5_sita_sensitivity as ps
print("\n=== sita_sensitivity (187) ===")
ps.main()

_csv.DictReader = _orig
print("\n완료")
