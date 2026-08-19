"""
OS rows의 VF 수치 + OCT 섹터 수치를 좌우 반전한 flip 버전 CSV 생성

VF flip: 24-2 grid 각 행 내에서 값을 역순으로 (OS→OD 좌표계 통일)
  Row 0 (p01-p04):  4점 reverse
  Row 1 (p05-p10):  6점 reverse
  Row 2 (p11-p18):  8점 reverse
  Row 3 (p19-p27):  9점 reverse
  Row 4 (p28-p36):  9점 reverse
  Row 5 (p37-p44):  8점 reverse
  Row 6 (p45-p50):  6점 reverse
  Row 7 (p51-p54):  4점 reverse

OCT flip (OS rows만):
  os_s_sup_t ↔ os_s_sup_n
  os_s_inf_t ↔ os_s_inf_n
"""
import csv
from pathlib import Path

import os as _os
from pathlib import Path as _Path
# 저장소 루트는 이 파일 위치에서 유도한다(개인 경로 하드코딩 금지).
# 다른 위치에서 쓰려면 HVF_ROOT 환경변수로 덮어쓴다.
_REPO_ROOT = _os.environ.get('HVF_ROOT', str(_Path(__file__).resolve().parents[2]))


BASE = Path(_REPO_ROOT)

# 24-2 grid row별 point 구간 (1-indexed label 번호)
VF_ROWS = [
    list(range(1, 5)),    # row 0: p01-p04
    list(range(5, 11)),   # row 1: p05-p10
    list(range(11, 19)),  # row 2: p11-p18
    list(range(19, 28)),  # row 3: p19-p27
    list(range(28, 37)),  # row 4: p28-p36
    list(range(37, 45)),  # row 5: p37-p44
    list(range(45, 51)),  # row 6: p45-p50
    list(range(51, 55)),  # row 7: p51-p54
]

VF_COLS = [f'p{i:02d}' for i in range(1, 55)]


def flip_vf(row_dict: dict) -> dict:
    """OS row의 VF 값을 각 grid row 내에서 reverse."""
    vals = [row_dict.get(col) for col in VF_COLS]  # 54개, None 가능
    flipped = list(vals)
    for group in VF_ROWS:
        idxs = [i - 1 for i in group]  # 0-indexed
        sub = [vals[i] for i in idxs]
        sub_rev = list(reversed(sub))
        for idx, val in zip(idxs, sub_rev):
            flipped[idx] = val
    result = dict(row_dict)
    for i, col in enumerate(VF_COLS):
        result[col] = flipped[i]
    return result


def flip_oct(row_dict: dict) -> dict:
    """OS row의 OCT temporal↔nasal 섹터 swap."""
    result = dict(row_dict)
    # sup: temporal↔nasal swap
    result['os_s_sup_t'] = row_dict.get('os_s_sup_n')
    result['os_s_sup_n'] = row_dict.get('os_s_sup_t')
    # inf: temporal↔nasal swap
    result['os_s_inf_t'] = row_dict.get('os_s_inf_n')
    result['os_s_inf_n'] = row_dict.get('os_s_inf_t')
    return result


def process_file(src: Path, dst: Path):
    rows = list(csv.DictReader(open(src, encoding='utf-8-sig')))
    if not rows:
        print(f'[SKIP] {src.name}: 빈 파일')
        return

    fieldnames = list(rows[0].keys())
    out_rows = []
    n_od, n_os = 0, 0

    for row in rows:
        if row.get('eye', '').upper() == 'OS':
            row = flip_vf(row)
            row = flip_oct(row)
            n_os += 1
        else:
            n_od += 1
        out_rows.append(row)

    with open(dst, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    print(f'저장: {dst.name}  (OD={n_od}, OS={n_os} → flip 적용)')


TARGETS = [
    ('ml_final_90d.csv',            'ml_final_90d_flip.csv'),
    ('ml_final_90d_excl_empty.csv', 'ml_final_90d_excl_empty_flip.csv'),
    ('ml_final_180d.csv',           'ml_final_180d_flip.csv'),
    ('ml_final_180d_excl_empty.csv','ml_final_180d_excl_empty_flip.csv'),
]

for src_name, dst_name in TARGETS:
    src = BASE / src_name
    dst = BASE / dst_name
    if not src.exists():
        print(f'[SKIP] {src_name}: 파일 없음')
        continue
    process_file(src, dst)

print('\n완료.')
