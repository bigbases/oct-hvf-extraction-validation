"""
SFA threshold_grid 이미지에서 50개 수치 추출 (하이브리드 OCR)
1. 고정 격자 좌표 (크기별 캘리브레이션 완료)
2. 전체 이미지 word-level OCR → 격자 위치 매핑 (빠름)
3. 미검출 셀만 개별 크롭 OCR (보완)
"""
import os, csv, sys, re
import pytesseract
from PIL import Image, ImageFilter, ImageEnhance
from pathlib import Path
import statistics

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
pytesseract.pytesseract.tesseract_cmd = r'C:\Program Files\Tesseract-OCR\tesseract.exe'

BASE    = _REPO_ROOT
ML_CSV  = os.path.join(BASE, 'ml_dataset.csv')
OUT_CSV = os.path.join(BASE, 'sfa_thresholds.csv')
DATASET_FINAL_CSV = os.path.join(BASE, 'kcc', 'dataset_final.csv')  # 원본 raw JPG 경로(vf_path) 조회용

sys.path.insert(0, os.path.join(BASE, 'src'))
from hvf.config import get as _cfg  # noqa: E402

import os as _os
from pathlib import Path as _Path
# 저장소 루트는 이 파일 위치에서 유도한다(개인 경로 하드코딩 금지).
# 다른 위치에서 쓰려면 HVF_ROOT 환경변수로 덮어쓴다.
_REPO_ROOT = _os.environ.get('HVF_ROOT', str(_Path(__file__).resolve().parents[2]))


LAYOUT_OD = {0:[2,3,4,5], 1:[1,2,3,4,5,6], 2:[0,1,2,3,4,5,6,7],
             3:[0,1,2,3,4,5,6,7,8], 4:[0,1,2,3,4,5,6,7,8],
             5:[0,1,2,3,4,5,6,7], 6:[1,2,3,4,5,6], 7:[2,3,4,5]}
LAYOUT_OS = {0:[2,3,4,5], 1:[1,2,3,4,5,6], 2:[0,1,2,3,4,5,6,7],
             3:[0,1,2,3,4,5,6,7,8], 4:[0,1,2,3,4,5,6,7,8],
             5:[0,1,2,3,4,5,6,7], 6:[1,2,3,4,5,6], 7:[2,3,4,5]}
POINT_LABELS = [f'p{i+1:02d}' for i in range(54)]

def make_gp(layout):
    return [(r,c) for r in range(8) for c in layout[r]]

GRID_POINTS_OD = make_gp(LAYOUT_OD)
GRID_POINTS_OS = make_gp(LAYOUT_OS)

# 크기별 고정 격자 중심 (검증 완료)
FIXED_GRIDS = {
    (385, 335): {
        'row': [59, 96, 132, 169, 206, 242, 279, 316],
        'col': [28, 64, 101, 137, 174, 212, 248, 285, 322],
        'half': 20,
        'tol': 16,   # 매핑 허용 거리 (px)
    },
    (513, 446): {
        'row': [80, 128, 177, 226, 274, 324, 372, 422],
        'col': [37, 85, 134, 183, 232, 281, 331, 378, 427],
        'half': 27,
        'tol': 22,
    },
}


def get_grid_info(img_w, img_h):
    """이미지 크기 → (row_c, col_c, half, tol). 모르는 크기면 비례 추정."""
    size = (img_w, img_h)
    if size in FIXED_GRIDS:
        g = FIXED_GRIDS[size]
        return g['row'], g['col'], g['half'], g['tol']
    # 비례 추정 (385×335 기준)
    scale = img_w / 385
    ref = FIXED_GRIDS[(385, 335)]
    row_c = [int(v * scale) for v in ref['row']]
    col_c = [int(v * scale) for v in ref['col']]
    half  = int(ref['half'] * scale)
    tol   = int(ref['tol']  * scale)
    return row_c, col_c, half, tol


def preprocess_full(img: Image.Image, scale: int = 2) -> Image.Image:
    img2 = img.resize((img.width*scale, img.height*scale), Image.LANCZOS)
    img2 = ImageEnhance.Contrast(img2).enhance(1.8)
    img2 = img2.filter(ImageFilter.SHARPEN)
    return img2


def full_image_ocr(img: Image.Image, scale: int = 2):
    """
    전체 이미지 word-level OCR.
    반환: list of (cx, cy, txt, conf)  ← 원본 좌표
    """
    img2 = preprocess_full(img, scale)
    cfg = r'--psm 11 --oem 1 -c tessedit_char_whitelist=0123456789<'
    data = pytesseract.image_to_data(img2, config=cfg, output_type=pytesseract.Output.DICT)

    results = []
    for txt, conf, x, y, w, h in zip(
            data['text'], data['conf'],
            data['left'], data['top'], data['width'], data['height']):
        txt = txt.strip()
        if not txt or int(conf) < 20:
            continue
        if not re.match(r'^<?[0-9]+$', txt):
            continue
        if w // scale > 55:   # 너무 넓으면 두 셀 합침
            continue
        cx = (x + w//2) // scale
        cy = (y + h//2) // scale
        results.append((cx, cy, txt, int(conf)))
    return results


def parse_value(txt: str):
    txt = txt.strip().replace(' ', '')
    if not txt:
        return None
    if '<' in txt:
        return -1
    m = re.search(r'\d+', txt)
    if not m:
        return None
    v = int(m.group())
    if v < 0:
        return -1
    if v > 50:
        return None
    return v


def cell_ocr(img: Image.Image, cx: int, cy: int, half: int) -> int | None:
    """단일 셀 크롭 OCR — psm7 → psm8 순으로 시도."""
    x0 = max(0, cx - half)
    y0 = max(0, cy - half)
    x1 = min(img.width,  cx + half)
    y1 = min(img.height, cy + half)
    if x1 <= x0 or y1 <= y0:
        return None

    crop = img.crop((x0, y0, x1, y1))
    scale = 3
    crop2 = crop.resize((crop.width*scale, crop.height*scale), Image.LANCZOS)
    crop2 = ImageEnhance.Contrast(crop2).enhance(2.5)
    crop2 = crop2.filter(ImageFilter.SHARPEN).filter(ImageFilter.SHARPEN)

    for psm in (7, 8, 6):
        cfg = f'--psm {psm} --oem 1 -c tessedit_char_whitelist=0123456789<'
        txt = pytesseract.image_to_string(crop2, config=cfg).strip()
        val = parse_value(txt)
        if val is not None:
            return val
    return None


def map_to_grid(detections, row_c, col_c, tol):
    """
    감지된 (cx, cy, txt) 목록을 격자 (row, col) → val 매핑.
    각 셀에서 가장 가까운 감지값 사용.
    """
    grid = {}
    for cx, cy, txt, conf in detections:
        # 가장 가까운 row/col 찾기
        best_r = min(range(len(row_c)), key=lambda i: abs(cy - row_c[i]))
        best_c = min(range(len(col_c)), key=lambda i: abs(cx - col_c[i]))
        dy = abs(cy - row_c[best_r])
        dx = abs(cx - col_c[best_c])
        if dy <= tol and dx <= tol:
            key = (best_r, best_c)
            # conf 높은 것 우선
            if key not in grid or conf > grid[key][1]:
                val = parse_value(txt)
                if val is not None:
                    grid[key] = (val, conf)
    return {k: v[0] for k, v in grid.items()}


def ocr_threshold(img_path: str, eye: str = 'OD', debug: bool = False) -> dict:
    eye = eye.upper()
    grid_points = GRID_POINTS_OD if eye == 'OD' else GRID_POINTS_OS
    empty = {lbl: None for lbl in POINT_LABELS}
    empty.update({'n_detected': 0, 'n_missing': 54})

    img = Image.open(img_path).convert('L')
    row_c, col_c, half, tol = get_grid_info(img.width, img.height)

    # Pass 1: 전체 이미지 OCR
    detections = full_image_ocr(img)
    grid = map_to_grid(detections, row_c, col_c, tol)

    if debug:
        print(f'  {img.width}x{img.height}  전체OCR감지={len(detections)}  격자매핑={len(grid)}')

    # Pass 2: 미검출 셀 개별 크롭 OCR
    result = {}
    n_detected = n_missing = 0

    for i, (r, c) in enumerate(grid_points):
        lbl = POINT_LABELS[i]
        if (r, c) in grid:
            result[lbl] = grid[(r, c)]
            n_detected += 1
        else:
            cy = row_c[r] if r < len(row_c) else None
            cx = col_c[c] if c < len(col_c) else None
            val = cell_ocr(img, cx, cy, half) if cx is not None and cy is not None else None
            result[lbl] = val
            if val is not None:
                n_detected += 1
                if debug:
                    print(f'  {lbl} 크롭OCR: ({cx},{cy}) → {val}')
            else:
                n_missing += 1

    result['n_detected'] = n_detected
    result['n_missing']  = n_missing
    return result


# ── row3/4(p19-p27, p28-p36) OD 전용 보정 패스 ──────────────────
# 2026-07-23 발견: 표준 grid_crop이 row3/4(맹점 인접 9점 두 행)의 col0(p19/p28)을
# 왼쪽 경계 밖으로 잘라낸다 — OD 전용(OS는 문제 없음, exact-match 100% 확인됨).
# src/hvf/ocr_vf.py._extract_row34_od 와 동일 로직을 이 파일의 기존 헬퍼
# (full_image_ocr/map_to_grid/cell_ocr/get_grid_info)로 이식. 이 스크립트는
# threshold_grid.png(이미 크롭됨)만 받으므로, 패치는 원본 SFA JPG(kcc/
# dataset_final.csv의 vf_path)를 별도로 열어 크롭창을 왼쪽으로 37px 이동해
# 재추출한다. 9개 값을 전부 얻었을 때만 표준 결과를 덮어쓰고, 실패하면
# 표준 결과(기존 값)를 그대로 둔다(회귀 방지 원칙).
_ROW34_OD_SHIFT_PX = 37  # col_c 평균 간격(~36.75px)에서 유도, 좌측 여유 확보


def extract_row34_od(orig_img_path: str) -> dict:
    """row3/4 전용 재추출(원본 SFA JPG 입력). 반환: {'row3': {p19..p27} or {}, 'row4': {p28..p36} or {}}."""
    img = Image.open(orig_img_path)
    W, H = img.width, img.height
    gx1, gy1, gx2, gy2 = _cfg('ocr', 'sfa', 'grid_crop')
    shift = _ROW34_OD_SHIFT_PX
    x1 = max(0, int(W * gx1) - shift)
    x2 = int(W * gx2) - shift
    y1 = int(H * gy1)
    y2 = int(H * gy2)
    shifted_img = img.crop((x1, y1, x2, y2)).convert('L')

    row_c, col_c, half, tol = get_grid_info(shifted_img.width, shifted_img.height)
    detections = full_image_ocr(shifted_img)
    grid = map_to_grid(detections, row_c, col_c, tol)

    out = {'row3': {}, 'row4': {}}
    for row_idx, base_p, key in [(3, 19, 'row3'), (4, 28, 'row4')]:
        vals = []
        for c in range(9):
            v = grid.get((row_idx, c))
            if v is None:
                v = cell_ocr(shifted_img, col_c[c], row_c[row_idx], half)
            vals.append(v)
        if all(v is not None for v in vals):
            out[key] = {f'p{base_p + i:02d}': v for i, v in enumerate(vals)}
    return out


def run_test(n=3, debug=False):
    ds = list(csv.DictReader(open(ML_CSV, encoding='utf-8-sig')))
    seen = set()
    for row in ds:
        key = row['sfa_dir']
        if not key or key in seen:
            continue
        seen.add(key)
        img_path = os.path.join(row['sfa_dir'], 'threshold_grid.png')
        if not os.path.exists(img_path):
            continue

        eye = row.get('eye', 'OD')
        res = ocr_threshold(img_path, eye=eye, debug=debug)
        print(f"\n[{row['patient_id']} {eye}]  detected={res['n_detected']}  missing={res['n_missing']}")

        gp = make_gp(LAYOUT_OD if eye == 'OD' else LAYOUT_OS)
        grid_display = [['   ']*8 for _ in range(8)]
        for i, (r, c) in enumerate(gp):
            v = res[POINT_LABELS[i]]
            grid_display[r][c] = ' <0' if v == -1 else (f'{v:3d}' if v is not None else '  ?')
        for r in range(8):
            print('  ' + ' '.join(grid_display[r]))

        if len(seen) >= n:
            break


if __name__ == '__main__':
    import sys as _sys

    if '--test' in _sys.argv:
        args = _sys.argv[_sys.argv.index('--test'):]
        n = int(args[1]) if len(args) > 1 and args[1].isdigit() else 3
        debug = '--debug' in _sys.argv
        run_test(n, debug=debug)
        _sys.exit(0)

    BATCH_LIMIT = int(_sys.argv[1]) if len(_sys.argv) > 1 and _sys.argv[1].isdigit() else None

    ds = list(csv.DictReader(open(ML_CSV, encoding='utf-8-sig')))
    vf_path_lookup = {(r['patient_id'], r['eye'], r['vf_date']): r['vf_path']
                       for r in csv.DictReader(open(DATASET_FINAL_CSV, encoding='utf-8-sig'))}

    fields = ['patient_id','eye','vf_date'] + POINT_LABELS + ['n_detected','n_missing']

    # 2026-07-23: row3/4 OD 패치로 눈당 OCR 호출량이 늘어 장시간 실행 시
    # Windows 리소스 고갈(STATUS_DLL_INIT_FAILED) 위험 — ocr_oct_values.py와
    # 동일하게 done_keys 이어쓰기 + 배치 한도 방식으로 재시작 가능하게 함.
    done_keys = set()
    if os.path.exists(OUT_CSV):
        for ex in csv.DictReader(open(OUT_CSV, encoding='utf-8-sig')):
            done_keys.add((ex['patient_id'], ex['eye'], ex['vf_date']))
        print(f'기존 {len(done_keys)}건 → 이어서 진행', flush=True)

    write_header = not os.path.exists(OUT_CSV) or len(done_keys) == 0
    f_out = open(OUT_CSV, 'a' if done_keys else 'w', newline='', encoding='utf-8-sig')
    writer = csv.DictWriter(f_out, fieldnames=fields, extrasaction='ignore')
    if write_header:
        writer.writeheader()

    seen = set()
    ok = err = 0
    total_missing = 0
    n_patch_applied = 0

    for i, row in enumerate(ds):
        sfa_dir = row['sfa_dir']
        key = (row['patient_id'], row['eye'], row['vf_date'])
        if not sfa_dir or key in seen or key in done_keys:
            continue
        seen.add(key)

        img_path = os.path.join(sfa_dir, 'threshold_grid.png')
        if not os.path.exists(img_path):
            err += 1
            continue

        try:
            eye = row.get('eye', 'OD')
            res = ocr_threshold(img_path, eye=eye)

            # row3/4 OD 전용 보정 패스
            if eye.upper() == 'OD':
                orig_path = vf_path_lookup.get(key)
                if orig_path and os.path.exists(orig_path):
                    patch = extract_row34_od(orig_path)
                    applied = {'row3': False, 'row4': False}
                    if patch['row3']:
                        res.update(patch['row3']); applied['row3'] = True
                    if patch['row4']:
                        res.update(patch['row4']); applied['row4'] = True
                    if applied['row3'] or applied['row4']:
                        n_patch_applied += 1

            res['patient_id'] = row['patient_id']
            res['eye']        = row['eye']
            res['vf_date']    = row['vf_date']
            # n_detected/n_missing 재계산(패치로 값이 바뀌었을 수 있으므로)
            res['n_detected'] = sum(1 for lbl in POINT_LABELS if res.get(lbl) is not None)
            res['n_missing']  = 54 - res['n_detected']
            writer.writerow(res)
            f_out.flush()
            total_missing += res['n_missing']
            ok += 1
        except Exception as e:
            err += 1
            print(f'ERR {row["patient_id"]}: {e}')

        if (ok + err) % 50 == 0:
            print(f'  {ok+err}/{len(ds)-len(done_keys)}(신규분) | ok:{ok} err:{err} | row34패치적용:{n_patch_applied}', flush=True)

        if BATCH_LIMIT is not None and (ok + err) >= BATCH_LIMIT:
            print(f'\n=== 배치 한도({BATCH_LIMIT}) 도달 — 정상 종료(래퍼가 이어서 재호출) ===', flush=True)
            break

    f_out.close()

    print(f'\n=== 완료(이번 실행분) ===')
    print(f'성공: {ok}  에러: {err}')
    print(f'row3/4 OD 패치 적용: {n_patch_applied}건')
    avg = total_missing / ok if ok else 0
    print(f'평균 미검출(이번 실행분): {avg:.1f}/54')
    print(f'저장: {OUT_CSV}')
