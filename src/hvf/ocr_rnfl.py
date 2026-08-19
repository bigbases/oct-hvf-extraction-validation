"""
RNFL quadrant + clock-hour OCR 모듈.

기존 ocr_rnfl_detail.py 에서 이관.
변경점:
  - parse_quadrants(path) → parse_quadrants(img_rgb): 파일 경로 대신 PIL Image
  - parse_clockhours(path) → parse_clockhours(img_rgb): 동일
  - extract_rnfl(rnfl_jpg_path, patient_id, eye, vf_date, oct_date) 추가
  OCR 후보 수집·할당 로직(collect_candidates / multi_stage / assign)은 무변경.
"""
import re, math
import pytesseract
from PIL import Image, ImageEnhance, ImageOps

from .config import configure_tesseract, get as _cfg
from .ocr_oct import ocr_cell

_CFG = '--psm 11 --oem 1 -c tessedit_char_whitelist=0123456789'


# ── OCR 후보 수집 (원본 동일) ─────────────────────────────────
def _collect_candidates(pil_l, conf_thr=10):
    ch_lo, ch_hi = _cfg('constants', 'plausibility', 'rnfl_clockhour')
    data = pytesseract.image_to_data(pil_l, config=_CFG, output_type=pytesseract.Output.DICT)
    out  = []
    for i, txt in enumerate(data['text']):
        m = re.search(r'\d+', txt.strip())
        if not m or data['conf'][i] < conf_thr:
            continue
        v = int(m.group())
        if v < ch_lo or v > ch_hi:
            continue
        cx = data['left'][i] + data['width'][i] // 2
        cy = data['top'][i]  + data['height'][i] // 2
        out.append((v, cx, cy))
    return out


def _multi_stage_candidates(img_rgb, scale=2):
    big  = img_rgb.resize((img_rgb.width * scale, img_rgb.height * scale), Image.LANCZOS)
    gray = big.convert('L')
    all_cands = []
    for contrast in (3.0, 2.0, 4.0, 1.5, 6.0):
        enh = ImageEnhance.Contrast(gray).enhance(contrast)
        all_cands.extend(_collect_candidates(enh, conf_thr=10))
    inv = ImageOps.invert(ImageEnhance.Contrast(gray).enhance(3.0))
    all_cands.extend(_collect_candidates(inv, conf_thr=10))

    dedup = []
    for v, cx, cy in all_cands:
        if not any(v == v2 and abs(cx - cx2) < 15 and abs(cy - cy2) < 15
                   for v2, cx2, cy2 in dedup):
            dedup.append((v, cx, cy))
    return [(v, cx // scale, cy // scale) for v, cx, cy in dedup]


# ── Quadrants (원본 로직 유지, img_rgb 받음) ─────────────────
def parse_quadrants(img_rgb: Image.Image) -> dict:
    keys  = ['S', 'T', 'I', 'N']
    empty = {f'{e}_{k}': None for e in ('od', 'os') for k in keys}
    if img_rgb is None:
        return empty
    W, H  = img_rgb.width, img_rgb.height
    # 2026-07-16 재보정 (v2, 셀단위 ocr_cell로 교체):
    # v1(좌표만 수정 + 전역 _multi_stage_candidates 최근접매칭)은 좌표는 맞았지만
    # 검출 신뢰도가 낮아 T/N 미검출·획소실 오독(58→36, 94→54 등)이 남았음.
    # 실측으로 크롭이 이웃 요소(범례/B-scan)를 물지 않는 것도 확인함 — 원인은
    # 전역 PSM11 탐지 신뢰도. GCA/summary와 동일한 셀단위 다중임계값 재시도
    # (ocr_cell)로 교체. 좌표 자체는 v1 실측값(od_S/os_S frac_x=0.430/0.866,
    # T/N행 frac_y=0.413, I행 frac_y=0.736) 그대로 사용.
    targets = {
        'od_S': (0.430, 0.088), 'od_T': (0.317, 0.413),
        'od_N': (0.542, 0.413), 'od_I': (0.430, 0.736),
        'os_S': (0.866, 0.088), 'os_N': (0.752, 0.413),
        'os_T': (0.978, 0.413), 'os_I': (0.866, 0.736),
    }
    # 2026-07-16 진단(미해결, 좌표 아님): I 필드 잔여오류(74%대)는 크롭 위치 확인
    # 결과 정상(박스가 정답 숫자 중앙에 정확히 위치, 잘림/이웃겹침 없음)이었음에도
    # 발생 — 순수 OCR 엔진 인식 실패. 깨끗한 2자리 이진화 이미지에 PSM 7/8/13 +
    # 임계값 4단계(70/100/120/150) 전부 테스트해도 오독(예: 실제 "98"을 PSM8은
    # "9"만, PSM7/13은 "45"/"46"/"9"로 반환). taxonomy="engine recognition
    # failure"(paper/manuscript/section/ch3_methods.tex) — 대응은 pipeline 수정이
    # 아니라 텍스트검출 기반 엔진(예: PaddleOCR) 대체가 필요.
    quad_lo, quad_hi = _cfg('constants', 'plausibility', 'rnfl_quadrant')
    result = dict(empty)
    for key, (xr, yr) in targets.items():
        cx, cy = int(W * xr), int(H * yr)
        result[key] = ocr_cell(img_rgb, cx, cy, hw=28, hh=18, min_val=quad_lo, max_val=quad_hi)
    return result


# ── Clock-hours (원본 로직 유지, img_rgb 받음) ───────────────
def parse_clockhours(img_rgb: Image.Image) -> dict:
    empty = {f'{e}_h{n:02d}': None for e in ('od', 'os') for n in range(1, 13)}
    if img_rgb is None:
        return empty
    W, H  = img_rgb.width, img_rgb.height
    cands = _multi_stage_candidates(img_rgb)

    cx_od = W * 0.206; cx_os = W * 0.817; cy_c = H * 0.522
    r_min = min(W, H) * 0.10; r_max = min(W, H) * 0.55

    cands_od, cands_os = [], []
    for v, cx, cy in cands:
        d_od = math.sqrt((cx - cx_od) ** 2 + (cy - cy_c) ** 2)
        d_os = math.sqrt((cx - cx_os) ** 2 + (cy - cy_c) ** 2)
        if d_od < d_os and r_min <= d_od <= r_max:
            ang = math.degrees(math.atan2(-(cy - cy_c), cx - cx_od))
            cands_od.append((v, ang, d_od))
        elif d_os < d_od and r_min <= d_os <= r_max:
            ang = math.degrees(math.atan2(-(cy - cy_c), cx - cx_os))
            cands_os.append((v, ang, d_os))

    def assign(cands_eye):
        target = {f'h{n:02d}': 90 - n * 30 for n in range(1, 13)}
        max_d  = 25
        pairs  = []
        for key, tgt in target.items():
            for idx, (v, ang, _) in enumerate(cands_eye):
                diff = abs((ang - tgt + 180) % 360 - 180)
                if diff < max_d:
                    pairs.append((diff, key, idx, v))
        pairs.sort()
        res = {}; used = set()
        for diff, key, idx, v in pairs:
            if key in res or idx in used:
                continue
            res[key] = v; used.add(idx)
        for n in range(1, 13):
            res.setdefault(f'h{n:02d}', None)
        return res

    od_h = assign(cands_od); os_h = assign(cands_os)
    return {**{f'od_{k}': v for k, v in od_h.items()},
            **{f'os_{k}': v for k, v in os_h.items()}}


# ── CV 룰 (원본 동일) ─────────────────────────────────────────
def cross_validate_rnfl(eye: str, quad: dict, clk: dict, avg_rnfl) -> list:
    flags = []
    qv    = [quad[f'{eye}_{d}'] for d in 'STIN']
    cv    = [clk[f'{eye}_h{n:02d}'] for n in range(1, 13)]
    qv_ok = [x for x in qv if x is not None]
    cv_ok = [x for x in cv if x is not None]

    # 2026-07-16: 절대오차 10 -> 5 로 강화(RNFL 50-120 범위 대비 10은 너무 느슨해
    # file2 OD의 오차 8.25(mean 65.75 vs avg 74) 같은 실오류를 놓쳤음).
    if avg_rnfl and len(qv_ok) >= 2:
        m = sum(qv_ok) / len(qv_ok)
        if abs(m - avg_rnfl) > 5:
            flags.append(f'CV_quad({len(qv_ok)}of4): mean={m:.1f} vs avg_rnfl={avg_rnfl}')

    if avg_rnfl and len(cv_ok) >= 6:
        m = sum(cv_ok) / len(cv_ok)
        if abs(m - avg_rnfl) > 5:
            flags.append(f'CV_clk({len(cv_ok)}of12): mean={m:.1f} vs avg_rnfl={avg_rnfl}')

    if quad[f'{eye}_S']:
        s_clk = [clk[f'{eye}_h{n:02d}'] for n in (10, 11, 12, 1, 2)]
        s_ok  = [x for x in s_clk if x is not None]
        if len(s_ok) >= 3:
            m = sum(s_ok) / len(s_ok)
            if abs(m - quad[f'{eye}_S']) > 15:
                flags.append(f'CV_S: quad_S={quad[f"{eye}_S"]} vs clk_top={m:.1f}')

    if quad[f'{eye}_I']:
        i_clk = [clk[f'{eye}_h{n:02d}'] for n in (4, 5, 6, 7, 8)]
        i_ok  = [x for x in i_clk if x is not None]
        if len(i_ok) >= 3:
            m = sum(i_ok) / len(i_ok)
            if abs(m - quad[f'{eye}_I']) > 15:
                flags.append(f'CV_I: quad_I={quad[f"{eye}_I"]} vs clk_bot={m:.1f}')

    return flags


RNFL_FIELDS = (
    ['patient_id', 'eye', 'vf_date', 'oct_date'] +
    [f'rnfl_q_{d.lower()}' for d in ['S', 'T', 'I', 'N']] +
    [f'rnfl_h{n:02d}' for n in range(1, 13)] +
    ['cv_flags', 'notes']
)


# ── 신규: 원본 RNFL JPG → rnfl_detail 행 ──────────────────────
def extract_rnfl(rnfl_jpg_path: str, avg_rnfl_od: float | None,
                 avg_rnfl_os: float | None,
                 patient_id: str, eye: str,
                 vf_date: str, oct_date: str) -> dict:
    """
    원본 RNFL JPG → {patient_id, eye, vf_date, oct_date, rnfl_q_*, rnfl_h*, cv_flags, notes}.
    avg_rnfl_* 는 cross-validation용 (oct_values_raw에서 전달).
    """
    configure_tesseract()

    row   = {k: '' for k in RNFL_FIELDS}
    row.update({'patient_id': patient_id, 'eye': eye,
                'vf_date': vf_date, 'oct_date': oct_date})
    notes = []

    if not rnfl_jpg_path:
        notes.append('no_rnfl_path'); row['notes'] = '|'.join(notes); return row

    try:
        rnfl_img = Image.open(rnfl_jpg_path).convert('RGB')
    except Exception as e:
        notes.append(f'open_error:{str(e)[:40]}'); row['notes'] = '|'.join(notes); return row

    # 영역 크롭
    def _cr(key_path):
        ratio = _cfg(*key_path)
        W, H  = rnfl_img.width, rnfl_img.height
        x1r, y1r, x2r, y2r = ratio
        return rnfl_img.crop((int(W*x1r), int(H*y1r), int(W*x2r), int(H*y2r)))

    quad_img = _cr(('ocr', 'rnfl', 'quadrants'))
    clk_img  = _cr(('ocr', 'rnfl', 'clockhours'))

    q = parse_quadrants(quad_img)
    c = parse_clockhours(clk_img)

    eye_l    = eye.lower()
    avg_rnfl = avg_rnfl_od if eye == 'OD' else avg_rnfl_os
    cv_flags = cross_validate_rnfl(eye_l, q, c, avg_rnfl)

    # 누락 카운트 경고
    n_q_miss = sum(1 for d in ['S', 'T', 'I', 'N'] if q[f'{eye_l}_{d}'] is None)
    n_c_miss = sum(1 for n in range(1, 13) if c[f'{eye_l}_h{n:02d}'] is None)
    if n_q_miss >= 2: notes.append(f'quad_miss={n_q_miss}/4')
    if n_c_miss >= 3: notes.append(f'clk_miss={n_c_miss}/12')

    row['cv_flags'] = '|'.join(cv_flags)
    row['notes']    = '|'.join(notes)
    for d in ['S', 'T', 'I', 'N']:
        val = q[f'{eye_l}_{d}']
        row[f'rnfl_q_{d.lower()}'] = val if val is not None else ''
    for n in range(1, 13):
        val = c[f'{eye_l}_h{n:02d}']
        row[f'rnfl_h{n:02d}'] = val if val is not None else ''

    return row
