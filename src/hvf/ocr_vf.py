"""
VF(SFA) 역치 OCR + test_pattern 추출 모듈.

기존 ocr_threshold.py 에서 이관.
변경점:
  - 파일 경로 하드코딩 제거 → config 로드
  - extract_vf(sfa_jpg_path, eye) 추가: 원본 SFA JPG에서 격자 크롭 + 헤더 추출
  - extract_test_pattern(banner_img) 추가: 배너에서 "24-2" / "10-2" 판별
  OCR 로직(get_grid_info / full_image_ocr / cell_ocr / map_to_grid)은 무변경.
"""
import re
import pytesseract
from PIL import Image, ImageFilter, ImageEnhance

from .config import configure_tesseract, get as _cfg

# ── 레이아웃 상수 ─────────────────────────────────────────────
LAYOUT_OD = {
    0: [2, 3, 4, 5],
    1: [1, 2, 3, 4, 5, 6],
    2: [0, 1, 2, 3, 4, 5, 6, 7],
    3: [0, 1, 2, 3, 4, 5, 6, 7, 8],   # 수정본: col 0 포함 (9점)
    4: [0, 1, 2, 3, 4, 5, 6, 7, 8],   # 수정본: col 0 포함 (9점)
    5: [0, 1, 2, 3, 4, 5, 6, 7],
    6: [1, 2, 3, 4, 5, 6],
    7: [2, 3, 4, 5],
}
LAYOUT_OS = {
    0: [2, 3, 4, 5],
    1: [1, 2, 3, 4, 5, 6],
    2: [0, 1, 2, 3, 4, 5, 6, 7],
    3: [0, 1, 2, 3, 4, 5, 6, 7, 8],
    4: [0, 1, 2, 3, 4, 5, 6, 7, 8],
    5: [0, 1, 2, 3, 4, 5, 6, 7],
    6: [1, 2, 3, 4, 5, 6],
    7: [2, 3, 4, 5],
}

POINT_LABELS = [f'p{i + 1:02d}' for i in range(54)]


def make_gp(layout):
    return [(r, c) for r in range(8) for c in layout[r]]


GRID_POINTS_OD = make_gp(LAYOUT_OD)
GRID_POINTS_OS = make_gp(LAYOUT_OS)

# ── 크기별 고정 격자 중심 (검증 완료, ocr_threshold.py 원본과 동일) ────────
FIXED_GRIDS = {
    (385, 335): {
        'row':  [59, 96, 132, 169, 206, 242, 279, 316],
        'col':  [28, 64, 101, 137, 174, 212, 248, 285, 322],
        'half': 20,
        'tol':  16,
    },
    (513, 446): {
        'row':  [80, 128, 177, 226, 274, 324, 372, 422],
        'col':  [37,  85, 134, 183, 232, 281, 331, 378, 427],
        'half': 27,
        'tol':  22,
    },
}


# ── 격자 정보 조회 ─────────────────────────────────────────────
def get_grid_info(img_w, img_h):
    size = (img_w, img_h)
    if size in FIXED_GRIDS:
        g = FIXED_GRIDS[size]
        return g['row'], g['col'], g['half'], g['tol']
    scale   = img_w / 385
    ref     = FIXED_GRIDS[(385, 335)]
    row_c   = [int(v * scale) for v in ref['row']]
    col_c   = [int(v * scale) for v in ref['col']]
    half    = int(ref['half'] * scale)
    tol     = int(ref['tol'] * scale)
    return row_c, col_c, half, tol


# ── OCR 유틸 (원본과 동일) ────────────────────────────────────
def _preprocess_full(img: Image.Image, scale: int = 2) -> Image.Image:
    img2 = img.resize((img.width * scale, img.height * scale), Image.LANCZOS)
    img2 = ImageEnhance.Contrast(img2).enhance(1.8)
    img2 = img2.filter(ImageFilter.SHARPEN)
    return img2


def _full_image_ocr(img: Image.Image, scale: int = 2):
    img2 = _preprocess_full(img, scale)
    cfg  = r'--psm 11 --oem 1 -c tessedit_char_whitelist=0123456789<'
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
        if w // scale > 55:
            continue
        cx = (x + w // 2) // scale
        cy = (y + h // 2) // scale
        results.append((cx, cy, txt, int(conf)))
    return results


def _parse_value(txt: str):
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


def _cell_ocr(img: Image.Image, cx: int, cy: int, half: int):
    x0 = max(0, cx - half); y0 = max(0, cy - half)
    x1 = min(img.width, cx + half); y1 = min(img.height, cy + half)
    if x1 <= x0 or y1 <= y0:
        return None
    crop  = img.crop((x0, y0, x1, y1))
    scale = 3
    crop2 = crop.resize((crop.width * scale, crop.height * scale), Image.LANCZOS)
    crop2 = ImageEnhance.Contrast(crop2).enhance(2.5)
    crop2 = crop2.filter(ImageFilter.SHARPEN).filter(ImageFilter.SHARPEN)
    for psm in (7, 8, 6):
        cfg = f'--psm {psm} --oem 1 -c tessedit_char_whitelist=0123456789<'
        txt = pytesseract.image_to_string(crop2, config=cfg).strip()
        val = _parse_value(txt)
        if val is not None:
            return val
    return None


def _map_to_grid(detections, row_c, col_c, tol):
    grid = {}
    for cx, cy, txt, conf in detections:
        best_r = min(range(len(row_c)), key=lambda i: abs(cy - row_c[i]))
        best_c = min(range(len(col_c)), key=lambda i: abs(cx - col_c[i]))
        dy = abs(cy - row_c[best_r]); dx = abs(cx - col_c[best_c])
        if dy <= tol and dx <= tol:
            key = (best_r, best_c)
            if key not in grid or conf > grid[key][1]:
                val = _parse_value(txt)
                if val is not None:
                    grid[key] = (val, conf)
    return {k: v[0] for k, v in grid.items()}


# ── 역치 OCR 코어 (기존 ocr_threshold 동일, img 객체 받음) ───────────────
def ocr_threshold_img(img: Image.Image, eye: str = 'OD') -> dict:
    """역치격자 PIL Image → {p01..p54, n_detected, n_missing}. OCR 로직 무변경."""
    eye          = eye.upper()
    grid_points  = GRID_POINTS_OD if eye == 'OD' else GRID_POINTS_OS
    empty        = {lbl: None for lbl in POINT_LABELS}
    empty.update({'n_detected': 0, 'n_missing': 54})

    img_l        = img.convert('L')
    row_c, col_c, half, tol = get_grid_info(img_l.width, img_l.height)

    detections   = _full_image_ocr(img_l)
    grid         = _map_to_grid(detections, row_c, col_c, tol)

    result       = {}
    n_detected   = n_missing = 0
    for i, (r, c) in enumerate(grid_points):
        lbl = POINT_LABELS[i]
        if (r, c) in grid:
            result[lbl] = grid[(r, c)]
            n_detected += 1
        else:
            cy  = row_c[r] if r < len(row_c) else None
            cx  = col_c[c] if c < len(col_c) else None
            val = _cell_ocr(img_l, cx, cy, half) if cx is not None else None
            result[lbl] = val
            if val is not None:
                n_detected += 1
            else:
                n_missing += 1

    result['n_detected'] = n_detected
    result['n_missing']  = n_missing
    return result


# ── 신규: test_pattern 추출 ───────────────────────────────────
def extract_test_pattern(banner_img: Image.Image) -> str:
    """
    SFA 헤더 배너 크롭에서 "24-2" / "10-2" 추출.
    반환: '24-2' | '10-2' | 'unknown'
    """
    img2 = banner_img.resize(
        (banner_img.width * 3, banner_img.height * 3), Image.LANCZOS
    ).convert('L')
    img2 = ImageEnhance.Contrast(img2).enhance(2.0)
    cfg  = '--psm 7 --oem 1 -c tessedit_char_whitelist=Central0123456789-ThresholdTest '
    txt  = pytesseract.image_to_string(img2, config=cfg)
    m    = re.search(r'(24|10)-2', txt)
    if m:
        return f'{m.group(1)}-2'
    # fallback: looser whitelist
    cfg2 = '--psm 6 --oem 1'
    txt2 = pytesseract.image_to_string(img2, config=cfg2)
    m2   = re.search(r'(24|10)-2', txt2)
    return f'{m2.group(1)}-2' if m2 else 'unknown'


# ── 신규: row3/4(p19-p27, p28-p36) OD 전용 보정 패스 ───────────
# 2026-07-23 발견: 표준 grid_crop이 row3/4(맹점 인접 9점 두 행)의 첫 값
# (col0, p19/p28)을 왼쪽 경계 밖으로 잘라낸다 — OD 전용(OS는 p19 exact-match
# 100%로 문제 없음 확인됨). 크롭 창을 한 컬럼 폭(~37px, col_c 간격 실측
# 36-38px 균일 확인)만큼 왼쪽으로 옮긴 별도 패스로 이 두 행만 재추출한다.
# 표준 col_c/grid_crop/다른 행은 무변경 — 이 패스가 9개 값을 전부 유효하게
# 얻었을 때만 표준 결과를 덮어쓰고, 실패하면 표준 결과(기존 8개)를 그대로
# 둔다(회귀 방지 원칙). 성공 여부는 result['row34_od_patch']에 기록한다.
_ROW34_OD_SHIFT_PX = 37  # col_c 평균 간격(~36.75px)에서 유도, 좌측 여유 확보


def _extract_row34_od(img: Image.Image, W: int, H: int) -> dict:
    """row3/4 전용 재추출. 반환: {'row3': {p19..p27} or {}, 'row4': {p28..p36} or {}}."""
    gx1, gy1, gx2, gy2 = _cfg('ocr', 'sfa', 'grid_crop')
    shift = _ROW34_OD_SHIFT_PX
    x1 = max(0, int(W * gx1) - shift)
    x2 = int(W * gx2) - shift
    y1 = int(H * gy1)
    y2 = int(H * gy2)
    shifted_img = img.crop((x1, y1, x2, y2)).convert('L')

    row_c, col_c, half, tol = get_grid_info(shifted_img.width, shifted_img.height)
    detections = _full_image_ocr(shifted_img)
    grid = _map_to_grid(detections, row_c, col_c, tol)

    out = {'row3': {}, 'row4': {}}
    for row_idx, base_p, key in [(3, 19, 'row3'), (4, 28, 'row4')]:
        vals = []
        for c in range(9):
            v = grid.get((row_idx, c))
            if v is None:
                # 표준 ocr_threshold_img과 동일한 셀 단위 폴백 (전체스캔 미검출 시)
                v = _cell_ocr(shifted_img, col_c[c], row_c[row_idx], half)
            vals.append(v)
        if all(v is not None for v in vals):
            out[key] = {f'p{base_p + i:02d}': v for i, v in enumerate(vals)}
    return out


# ── 신규: 원본 SFA JPG → 역치 + test_pattern ──────────────────
def extract_vf(sfa_jpg_path: str, eye: str) -> dict:
    """
    원본 SFA JPG 경로 + eye → {test_pattern, p01..p54, n_detected, n_missing}.

    이미지 내 격자 크롭은 config ocr.sfa.grid_crop 비율 사용.
    test_pattern은 config ocr.sfa.banner_crop 비율 사용.
    """
    configure_tesseract()

    img = Image.open(sfa_jpg_path)
    W, H = img.width, img.height

    # 역치격자 크롭
    gx1, gy1, gx2, gy2 = _cfg('ocr', 'sfa', 'grid_crop')
    grid_img = img.crop((int(W * gx1), int(H * gy1), int(W * gx2), int(H * gy2)))

    # 테스트 패턴 배너 크롭
    bx1, by1, bx2, by2 = _cfg('ocr', 'sfa', 'banner_crop')
    banner_img = img.crop((int(W * bx1), int(H * by1), int(W * bx2), int(H * by2)))

    result               = ocr_threshold_img(grid_img, eye)
    result['test_pattern'] = extract_test_pattern(banner_img)

    # row3/4 OD 전용 보정 패스 (2026-07-23)
    if eye.upper() == 'OD':
        patch = _extract_row34_od(img, W, H)
        applied = {'row3': False, 'row4': False}
        if patch['row3']:
            result.update(patch['row3'])
            applied['row3'] = True
        if patch['row4']:
            result.update(patch['row4'])
            applied['row4'] = True
        result['row34_od_patch'] = applied

    return result
