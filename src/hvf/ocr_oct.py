"""
OCT GCA/RNFL 수치 OCR + signal_strength 추출 모듈.

기존 ocr_oct_values.py 에서 이관.
변경점:
  - parse_*(path) → parse_*(img_rgb): 파일 경로 대신 PIL Image 객체 수신
  - _crop_region(img, ratio) 헬퍼 추가: 원본 JPG에서 영역 인라인 크롭
  - extract_signal_strength(header_img) 추가: GCA/RNFL 헤더 "신호 강도: X/10 X/10"
  - extract_oct(gca_jpg, rnfl_jpg) 추가: 원본 JPG 두 장 → 전체 추출
  OCR 셀 추출 로직(ocr_cell / cross_validate)은 무변경.
"""
import re, math
import pytesseract
from PIL import Image, ImageEnhance, ImageOps, ImageChops

from .config import configure_tesseract, get as _cfg

# ── Sanity check 범위 ──────────────────────────────────────────
# 2026-07-23: 하드코딩 값(avg_rnfl 40-130 등)이 GT(수동검수본) 관측범위 밖의
# 진짜 값을 폐기하고 있었음이 확인되어 config.plausibility로 교체.
# (이 dict는 _normalize_row 최종 게이트 — parse_* 함수 내부의 1차 min_val/
# max_val과 반드시 동일한 값을 참조해야 함. 두 레이어가 다른 값을 쓰면
# 1차를 고쳐도 여기서 다시 걸러짐 — 2026-07-23 실측으로 발견.)
def _ranges():
    p = lambda k: tuple(x if x is not None else float('-inf') if i == 0 else float('inf')
                        for i, x in enumerate(_cfg('constants', 'plausibility', k)))
    return {
        'avg_rnfl': p('rnfl_avg'),
        'min_gcl':  p('min_gcl'),
        'avg_gcl':  p('avg_gcl'),
        'vert_cd':  p('vert_cd'),
        '_s_':      p('gca_sector'),
    }

_OCR_CFG8 = '--psm 8 --oem 1 -c tessedit_char_whitelist=0123456789.'
_OCR_CFG7 = '--psm 7 --oem 1 -c tessedit_char_whitelist=0123456789.'
_OCR_CFG6 = '--psm 6 --oem 1 -c tessedit_char_whitelist=0123456789.'
_PAD = 10


def apply_sanity(key, val, flagged=None):
    if val is None:
        return None
    for pattern, (lo, hi) in _ranges().items():
        if pattern in key:
            if val < lo or val > hi:
                if flagged is not None:
                    flagged.append((key, val, lo, hi))
                return None
    return val


def _pad(pil_l):
    out = Image.new('L', (pil_l.width + 2 * _PAD, pil_l.height + 2 * _PAD), 255)
    out.paste(pil_l, (_PAD, _PAD))
    return out


def _rgb_binary(crop_rgb, thr=70):
    r, g, b = crop_rgb.split()
    mx = ImageChops.lighter(ImageChops.lighter(r, g), b)
    return mx.point(lambda x: 0 if x < thr else 255, 'L')


def _try_val(pil_img, cfg, min_val=None, max_val=None):
    txt = pytesseract.image_to_string(pil_img, config=cfg).strip()
    m   = re.search(r'\d+\.?\d*', txt)
    if not m:
        return None
    v = float(m.group())
    if min_val is not None and v < min_val:
        return None
    if max_val is not None and v > max_val:
        return None
    return v


# ── OCR 셀 추출 (원본 ocr_oct_values.py ocr_cell 동일) ────────
def ocr_cell(img_rgb, cx, cy, hw, hh, hh_below=None, psm8_first=False,
             min_val=None, max_val=None):
    hb = hh_below if hh_below is not None else hh
    x0, y0 = max(0, cx - hw), max(0, cy - hh)
    x1, y1 = min(img_rgb.width, cx + hw), min(img_rgb.height, cy + hb)
    crop     = img_rgb.crop((x0, y0, x1, y1))
    big_rgb  = crop.resize((crop.width * 4, crop.height * 4), Image.LANCZOS)
    big_gray = big_rgb.convert('L')
    kw = dict(min_val=min_val, max_val=max_val)
    s1 = ((_OCR_CFG8,), (_OCR_CFG7,)) if psm8_first else ((_OCR_CFG7,), (_OCR_CFG8,))

    for (cfg,) in s1:
        for thr in (70, 100, 120, 150):
            v = _try_val(_pad(_rgb_binary(big_rgb, thr=thr)), cfg, **kw)
            if v is not None:
                return v

    for contrast in (3.0, 2.0, 4.0, 1.5):
        c = _pad(ImageEnhance.Contrast(big_gray).enhance(contrast))
        for (cfg,) in s1:
            v = _try_val(c, cfg, **kw)
            if v is not None:
                return v
        v = _try_val(c, _OCR_CFG6, **kw)
        if v is not None:
            return v

    inv = _pad(ImageOps.invert(_rgb_binary(big_rgb, thr=60)))
    for (cfg,) in s1:
        v = _try_val(inv, cfg, **kw)
        if v is not None:
            return v

    for contrast in (3.0, 4.0, 2.0):
        ig = _pad(ImageOps.invert(ImageEnhance.Contrast(big_gray).enhance(contrast)))
        for (cfg,) in s1:
            v = _try_val(ig, cfg, **kw)
            if v is not None:
                return v

    return None


# ── 영역 크롭 헬퍼 ────────────────────────────────────────────
def _crop_region(img: Image.Image, ratio: list) -> Image.Image:
    """비율 [x1,y1,x2,y2] → 픽셀 크롭. RGB 유지."""
    x1r, y1r, x2r, y2r = ratio
    W, H = img.width, img.height
    return img.crop((int(W * x1r), int(H * y1r), int(W * x2r), int(H * y2r)))


# ── 1. GCL 두께표 (원본 로직 유지, img_rgb 받음) ─────────────
def parse_gcl_table(img_rgb: Image.Image) -> dict:
    empty = {k: None for k in ['avg_gcl_od', 'avg_gcl_os', 'min_gcl_od', 'min_gcl_os']}
    if img_rgb is None:
        return empty
    W, H = img_rgb.width, img_rgb.height
    avg_lo, avg_hi = _cfg('constants', 'plausibility', 'avg_gcl')
    min_lo, min_hi = _cfg('constants', 'plausibility', 'min_gcl')
    cells = {
        'avg_gcl_od': (int(W * 0.64), int(H * 0.42), 45, 22, avg_lo, avg_hi),
        'avg_gcl_os': (int(W * 0.85), int(H * 0.42), 45, 22, avg_lo, avg_hi),
        'min_gcl_od': (int(W * 0.64), int(H * 0.81), 45, 22, min_lo, min_hi),
        'min_gcl_os': (int(W * 0.85), int(H * 0.81), 45, 22, min_lo, min_hi),
    }
    result = {}
    for key, (cx, cy, hw, hh, lo, hi) in cells.items():
        x0, y0 = max(0, cx - hw), max(0, cy - hh)
        x1, y1 = min(W, cx + hw), min(H, cy + hh)
        cell = img_rgb.crop((x0, y0, x1, y1))
        big  = cell.resize((cell.width * 4, cell.height * 4), Image.LANCZOS)
        # 1차: 원래 방식 (grayscale+contrast, PSM6→PSM8) — 흰/회색 배경에 강함
        v = None
        for contrast in (3.0, 2.0, 4.0, 1.5):
            c = ImageEnhance.Contrast(big.convert('L')).enhance(contrast)
            v = _try_val(c, _OCR_CFG6, min_val=lo, max_val=hi)
            if v is not None:
                break
            v = _try_val(c, _OCR_CFG8, min_val=lo, max_val=hi)
            if v is not None:
                break
        # 2차: _rgb_binary 폴백 — 녹색/색상 배경에 강함
        if v is None:
            v = ocr_cell(img_rgb, cx, cy, hw, hh, min_val=lo, max_val=hi)
        # 3차: 밝은 녹색 배경용 고threshold — 패딩 없이
        if v is None:
            for thr in (170, 190):
                bi = _rgb_binary(big, thr=thr)
                v = _try_val(bi, _OCR_CFG8, min_val=lo, max_val=hi)
                if v is not None:
                    break
        result[key] = int(v) if v is not None else None
    return result


# ── 2. GCA 섹터 (원본 로직 유지, img_rgb 받음) ───────────────
def parse_sectors(img_rgb: Image.Image, eye: str = 'OD') -> dict:
    keys  = ['s_sup', 's_sup_t', 's_inf_t', 's_inf', 's_inf_n', 's_sup_n']
    empty = {f'{eye.lower()}_{k}': None for k in keys}
    if img_rgb is None:
        return empty
    W, H = img_rgb.width, img_rgb.height
    cx_c, cy_c = W // 2, H // 2
    r = min(W, H) * 0.35
    # 2026-07-16 최종 결론 (하루 종일 검증, 두 번 뒤집힘 — 근거 남김):
    # GCA 6섹터 파이차트는 스키마틱 오버레이로, B-scan/두께맵(사진)과 달리 T 라벨
    # 위치가 양안 동일(30도)함을 대표 리포트 육안 대조로 확인함(두께맵 시신경유두
    # 위치로 정한 해부학적 temporal이 30도 값과 일치, eye 분기 불필요).
    # GT(oct_tabular_90d.csv)는 row의 eye가 OD/OS냐에 따라 os_s_* 를 다르게 기록
    # (같은 리포트인데 OD행·OS행의 os_s_sup_t가 서로 다름) — 이게 지금까지의 모든
    # 불일치의 실제 원인. 이 identical-angle 코드 + OS행에서만 스왑하는 채점
    # 로직(scripts/verify_gca_os_mirroring.py)을 함께 쓰면 GT와 ~83% 일치
    # (OD행 그대로=83.5%, OS행 스왑=82.8%; paper 82.7%와 부합).
    # eye 분기(mirroring)를 코드에 넣으면 오히려 OD행 대비 5.4%로 떨어짐 — 실측 확인.
    angles_deg = [90, 30, 330, 270, 210, 150]
    sec_lo, sec_hi = _cfg('constants', 'plausibility', 'gca_sector')
    result = {}
    for key, ang in zip(keys, angles_deg):
        rad = math.radians(ang)
        px  = int(cx_c + r * math.cos(rad))
        py  = int(cy_c - r * math.sin(rad))
        hw, hh = (55, 20) if ang in (90, 270) else (55, 30)
        # TODO(미해결, 2026-07-16 진단): s_inf(270도) 정확도가 sup(90도)보다 훨씬
        # 낮음(od_s_sup 89.0% vs od_s_inf 32.2%). 색상/자릿수 기각(직접측정: green
        # 배경에서도 sup 99.1% vs inf 41.9%). 크롭 hw를 좁히면 일부 케이스는 복구되나
        # (예: hw=40에서 76 정확 추출) 다른 케이스는 hw 0~55 전 범위에서 실패 —
        # 이건 폭이 아니라 중심좌표(px,py) 자체가 숫자를 벗어난 것으로 추정됨.
        # 고정 오프셋으로는 리포트별 렌더링 변동을 못 따라감. 다음 단계: 색/경계
        # 기반 동적 크롭(예: 웨지 색상 세그멘테이션으로 숫자 중심 재탐지) 필요.
        # 결과: results/gca_os_mirroring_verification.json, taxonomy=
        # "geometric localization failure"(paper/manuscript/section/ch3_methods.tex).
        v = ocr_cell(img_rgb, px, py, hw, hh, psm8_first=True, min_val=sec_lo, max_val=sec_hi)
        result[f'{eye.lower()}_{key}'] = int(v) if v is not None else None
    return result


# ── 3. RNFL 요약표 (원본 로직 유지, img_rgb 받음) ────────────
def parse_rnfl_summary(img_rgb: Image.Image) -> dict:
    empty = {'od_avg_rnfl': None, 'os_avg_rnfl': None,
             'od_vert_cd': None,  'os_vert_cd': None}
    if img_rgb is None:
        return empty
    W, H = img_rgb.width, img_rgb.height
    # 실측 보정값 (560×348 크롭 기준):
    #   od_avg_rnfl: xr=0.46, yr=0.18  (OD 행 중심)
    #   os_avg_rnfl: xr=0.73, yr=0.16  (OS 셀이 OD보다 약 2% 위; hw=28로 제한)
    #   vert_cd: yr=0.77  (수직C/D 6번째 행; od_vert_cd=0.66 실측 확인)
    # (cx, cy, hw, hh, min_val, max_val)
    rnfl_lo, rnfl_hi = _cfg('constants', 'plausibility', 'rnfl_avg')
    cd_lo, cd_hi = _cfg('constants', 'plausibility', 'vert_cd')
    cells = {
        'od_avg_rnfl': (int(W * 0.46), int(H * 0.18), 40, 16, rnfl_lo, rnfl_hi),
        'os_avg_rnfl': (int(W * 0.73), int(H * 0.16), 28, 18, rnfl_lo, rnfl_hi),
        'od_vert_cd':  (int(W * 0.49), int(H * 0.77), 28, 16, cd_lo, cd_hi),
        'os_vert_cd':  (int(W * 0.76), int(H * 0.77), 28, 16, cd_lo, cd_hi),
    }
    result = {}
    for key, (cx, cy, hw, hh, min_val, max_val) in cells.items():
        v = ocr_cell(img_rgb, cx, cy, hw, hh, min_val=min_val, max_val=max_val)
        result[key] = v
    return result


# ── 4. Quadrants CV (원본 로직 유지, img_rgb 받음) ───────────
def parse_quadrants_cv(img_rgb: Image.Image) -> dict:
    empty = {'od_S': None, 'od_T': None, 'od_I': None, 'od_N': None,
             'os_S': None, 'os_T': None, 'os_I': None, 'os_N': None}
    if img_rgb is None:
        return empty
    W, H = img_rgb.width, img_rgb.height
    # 2026-07-16 재보정: 기존 좌표(od_S 0.21/0.13 등)가 실제 숫자 위치와 안 맞아
    # od_S/od_T/od_N/os_S/os_N/os_T가 전부 B-scan 영역을 검색하고 있었음(=결측
    # 54-65%의 주 원인). dw_images 8건 Tesseract 실측 좌표(od_S/os_S frac_x=
    # 0.430/0.866 등)로 대체하고, 전역 PSM11 탐지+최근접매칭 대신 GCA/summary와
    # 동일한 셀단위 ocr_cell(다중 임계값 재시도)로 교체 — 더 견고함.
    # od_I/os_I: 이 파이차트 옆에는 수치 라벨 자체가 인쇄되지 않음(문자 I만 있음) —
    # 좌표 문제가 아니라 이 크롭 영역에 값이 없는 것. summary 표 등 다른 위치 확인 필요.
    targets = {
        'od_S': (0.430, 0.137), 'od_T': (0.317, 0.639), 'od_N': (0.542, 0.639),
        'os_S': (0.866, 0.137), 'os_N': (0.752, 0.639), 'os_T': (0.978, 0.639),
    }
    result = dict(empty)
    for key, (xr, yr) in targets.items():
        cx, cy = int(W * xr), int(H * yr)
        result[key] = ocr_cell(img_rgb, cx, cy, hw=24, hh=14, min_val=30, max_val=250)
    return result


# ── 5. Clock hours CV (원본 로직 유지, img_rgb 받음) ─────────
def parse_clockhours_cv(img_rgb: Image.Image) -> dict:
    empty_od = {f'od_h{i + 1:02d}': None for i in range(12)}
    empty_os = {f'os_h{i + 1:02d}': None for i in range(12)}
    if img_rgb is None:
        return {**empty_od, **empty_os}
    W, H = img_rgb.width, img_rgb.height
    cx_od = int(W * 0.206); cx_os = int(W * 0.817); cy = int(H * 0.522)
    r_min, r_max = 60, 135
    enh = ImageEnhance.Contrast(img_rgb.convert('L')).enhance(6.0)
    cfg = '--psm 11 --oem 1 -c tessedit_char_whitelist=0123456789'
    data = pytesseract.image_to_data(enh, config=cfg, output_type=pytesseract.Output.DICT)
    cands_od, cands_os = [], []
    for i, txt in enumerate(data['text']):
        m = re.search(r'\d+', txt.strip())
        if not m or data['conf'][i] < 20:
            continue
        v = int(m.group())
        if v < 5 or v > 300:
            continue
        bx = data['left'][i] + data['width'][i] // 2
        by = data['top'][i] + data['height'][i] // 2
        d_od = math.sqrt((bx - cx_od) ** 2 + (by - cy) ** 2)
        d_os = math.sqrt((bx - cx_os) ** 2 + (by - cy) ** 2)
        if d_od < d_os and r_min <= d_od <= r_max:
            cands_od.append((v, math.degrees(math.atan2(-(by - cy), bx - cx_od)), d_od))
        elif d_os < d_od and r_min <= d_os <= r_max:
            cands_os.append((v, math.degrees(math.atan2(-(by - cy), bx - cx_os)), d_os))

    def assign(cands):
        target = {f'h{n:02d}': 90 - n * 30 for n in range(1, 13)}
        res = {}; used = set()
        for key, tgt in target.items():
            best_v, best_d, best_i = None, 20, -1
            for idx, (v, ang, _) in enumerate(cands):
                if idx in used:
                    continue
                diff = abs((ang - tgt + 180) % 360 - 180)
                if diff < best_d:
                    best_d, best_v, best_i = diff, v, idx
            if best_i >= 0:
                res[key] = best_v; used.add(best_i)
            else:
                res[key] = None
        return res

    od_h = assign(cands_od); os_h = assign(cands_os)
    return {**{f'od_{k}': v for k, v in od_h.items()},
            **{f'os_{k}': v for k, v in os_h.items()}}


# ── 신규: signal_strength 추출 ────────────────────────────────
def extract_signal_strength(header_img: Image.Image) -> dict:
    """
    GCA 또는 RNFL 헤더 크롭에서 "신호 강도: OD_val/10  OS_val/10" 추출.
    반환: {'ss_od': int|None, 'ss_os': int|None}
    """
    if header_img is None:
        return {'ss_od': None, 'ss_os': None}
    img2 = header_img.resize(
        (header_img.width * 3, header_img.height * 3), Image.LANCZOS
    ).convert('L')
    img2 = ImageEnhance.Contrast(img2).enhance(2.5)
    cfg  = '--psm 6 --oem 1 -c tessedit_char_whitelist=0123456789/'
    txt  = pytesseract.image_to_string(img2, config=cfg)
    vals = re.findall(r'(\d{1,2})/10', txt)
    return {
        'ss_od': int(vals[0]) if len(vals) >= 1 else None,
        'ss_os': int(vals[1]) if len(vals) >= 2 else None,
    }


# ── Cross-validation (원본 동일) ──────────────────────────────
_CV_THR = 5


def cross_validate(norm_row, quad, clk) -> list:
    flags = []

    def _mean(vals):
        v = [x for x in vals if x is not None]
        return sum(v) / len(v) if v else None

    sec_keys = ['sup', 'sup_t', 'inf_t', 'inf', 'inf_n', 'sup_n']
    for eye in ('od', 'os'):
        avg_gcl = norm_row.get(f'avg_gcl_{eye}')
        secs    = [norm_row.get(f'{eye}_s_{k}') for k in sec_keys]
        m_sec   = _mean(secs)
        if avg_gcl is not None and m_sec is not None and len([x for x in secs if x]) >= 3:
            diff = abs(avg_gcl - m_sec)
            if diff > _CV_THR:
                flags.append(f'CV1_{eye}: avg_gcl={avg_gcl} vs sec_mean={m_sec:.1f} (diff={diff:.1f})')

        avg_rnfl = norm_row.get(f'{eye}_avg_rnfl')
        q_vals   = [quad.get(f'{eye}_{d}') for d in ('S', 'T', 'I', 'N')]
        m_quad   = _mean(q_vals)
        if avg_rnfl is not None and m_quad is not None and len([x for x in q_vals if x]) >= 2:
            diff = abs(avg_rnfl - m_quad)
            if diff > _CV_THR:
                flags.append(f'CV2_{eye}: avg_rnfl={avg_rnfl} vs quad_mean={m_quad:.1f} (diff={diff:.1f})')

        min_gcl = norm_row.get(f'min_gcl_{eye}')
        if min_gcl is not None and m_sec is not None:
            min_sec = min((x for x in secs if x is not None), default=None)
            if min_sec is not None and min_gcl > min_sec + _CV_THR:
                flags.append(f'CV3_{eye}: min_gcl_table={min_gcl} > min_sec={min_sec}')

        clk_vals = [clk.get(f'{eye}_h{n:02d}') for n in range(1, 13)]
        m_clk    = _mean(clk_vals)
        n_clk    = len([x for x in clk_vals if x])
        if avg_rnfl is not None and m_clk is not None and n_clk >= 6:
            diff = abs(avg_rnfl - m_clk)
            if diff > _CV_THR:
                flags.append(f'CV4_{eye}: avg_rnfl={avg_rnfl} vs clk_mean={m_clk:.1f} (diff={diff:.1f},n={n_clk})')

    return flags


def _normalize_row(pid, eye, vals, quad, clk, flagged_log=None):
    if flagged_log is None:
        flagged_log = []
    r = {'patient_id': pid, 'eye': eye}
    for k in FIELDS:
        if k in ('patient_id', 'eye', 'cv_flags'):
            continue
        r[k] = apply_sanity(k, vals.get(k), flagged_log)
    r['cv_flags'] = '|'.join(cross_validate(r, quad, clk))
    return r


FIELDS = (
    ['patient_id', 'eye'] +
    ['avg_gcl_od', 'avg_gcl_os', 'min_gcl_od', 'min_gcl_os'] +
    [f'od_s_{k}' for k in ['sup', 'sup_t', 'inf_t', 'inf', 'inf_n', 'sup_n']] +
    [f'os_s_{k}' for k in ['sup', 'sup_t', 'inf_t', 'inf', 'inf_n', 'sup_n']] +
    ['od_avg_rnfl', 'os_avg_rnfl', 'od_vert_cd', 'os_vert_cd', 'cv_flags']
)


# ── 신규: 원본 GCA + RNFL JPG → 전체 OCT 수치 + signal_strength ──
def extract_oct(gca_jpg_path: str, rnfl_jpg_path: str,
                patient_id: str, eye: str) -> dict:
    """
    원본 GCA JPG + RNFL JPG 경로 → oct_values 행 + signal_strength 행.
    반환: {'oct': dict(FIELDS + ss_*), 'ss_gca': {'ss_od', 'ss_os'}, 'ss_rnfl': {...}}
    """
    configure_tesseract()

    # GCA 열기
    gca_img  = Image.open(gca_jpg_path).convert('RGB') if gca_jpg_path else None
    rnfl_img = Image.open(rnfl_jpg_path).convert('RGB') if rnfl_jpg_path else None

    # 영역 크롭 (GCA)
    gcl_crop   = _crop_region(gca_img,  _cfg('ocr', 'gca', 'gcl_table'))   if gca_img  else None
    od_s_crop  = _crop_region(gca_img,  _cfg('ocr', 'gca', 'od_sectors'))  if gca_img  else None
    os_s_crop  = _crop_region(gca_img,  _cfg('ocr', 'gca', 'os_sectors'))  if gca_img  else None
    gca_ss_crop= _crop_region(gca_img,  _cfg('ocr', 'signal_strength', 'header_crop')) if gca_img else None

    # 영역 크롭 (RNFL)
    rnfl_sum_crop = _crop_region(rnfl_img, _cfg('ocr', 'rnfl', 'summary'))    if rnfl_img else None
    quad_crop     = _crop_region(rnfl_img, _cfg('ocr', 'rnfl', 'quadrants'))  if rnfl_img else None
    clk_crop      = _crop_region(rnfl_img, _cfg('ocr', 'rnfl', 'clockhours')) if rnfl_img else None
    rnfl_ss_crop  = _crop_region(rnfl_img, _cfg('ocr', 'signal_strength', 'header_crop')) if rnfl_img else None

    # OCR 수치 추출
    gcl    = parse_gcl_table(gcl_crop)
    od_s   = parse_sectors(od_s_crop, 'OD')
    os_s   = parse_sectors(os_s_crop, 'OS')
    rnfl_s = parse_rnfl_summary(rnfl_sum_crop)
    quad   = parse_quadrants_cv(quad_crop)
    clk    = parse_clockhours_cv(clk_crop)

    # signal strength
    ss_gca  = extract_signal_strength(gca_ss_crop)
    ss_rnfl = extract_signal_strength(rnfl_ss_crop)

    vals = {**gcl, **od_s, **os_s, **rnfl_s}
    oct_row = _normalize_row(patient_id, eye, vals, quad, clk)
    oct_row['ss_gca_od']  = ss_gca['ss_od']
    oct_row['ss_gca_os']  = ss_gca['ss_os']
    oct_row['ss_rnfl_od'] = ss_rnfl['ss_od']
    oct_row['ss_rnfl_os'] = ss_rnfl['ss_os']

    return oct_row
