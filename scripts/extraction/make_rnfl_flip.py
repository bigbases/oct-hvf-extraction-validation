"""
rnfl_detail.csv / clockhours.csv OS→OD flip 버전 생성

OS 행만 적용:
  quadrant: rnfl_q_t ↔ rnfl_q_n
  clock hours (좌우 대칭):
    H01(12시), H07(6시) 고정
    H02↔H12, H03↔H11, H04↔H10, H05↔H09, H06↔H08
"""
import csv
from pathlib import Path

import os as _os
from pathlib import Path as _Path
# 저장소 루트는 이 파일 위치에서 유도한다(개인 경로 하드코딩 금지).
# 다른 위치에서 쓰려면 HVF_ROOT 환경변수로 덮어쓴다.
_REPO_ROOT = _os.environ.get('HVF_ROOT', str(_Path(__file__).resolve().parents[2]))


BASE = Path(_REPO_ROOT)

QUAD_SWAP = [('rnfl_q_t', 'rnfl_q_n')]
HOUR_SWAP = [('rnfl_h02','rnfl_h12'), ('rnfl_h03','rnfl_h11'),
             ('rnfl_h04','rnfl_h10'), ('rnfl_h05','rnfl_h09'),
             ('rnfl_h06','rnfl_h08')]


def flip_row(row: dict) -> dict:
    r = dict(row)
    for a, b in QUAD_SWAP + HOUR_SWAP:
        r[a], r[b] = row[b], row[a]
    return r


def process(src_name: str, dst_name: str):
    rows = list(csv.DictReader(open(BASE / src_name, encoding='utf-8-sig')))
    fieldnames = list(rows[0].keys())
    out, n_os = [], 0
    for row in rows:
        if row.get('eye', '').upper() == 'OS':
            out.append(flip_row(row))
            n_os += 1
        else:
            out.append(row)
    with open(BASE / dst_name, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out)
    print(f'{dst_name}: {len(rows)}행 저장 (OS flip {n_os}건)')


process('rnfl_detail.csv',     'rnfl_detail_flip.csv')
process('rnfl_detail_180d.csv','rnfl_detail_180d_flip.csv')
process('clockhours.csv',      'clockhours_flip.csv')
print('완료')
