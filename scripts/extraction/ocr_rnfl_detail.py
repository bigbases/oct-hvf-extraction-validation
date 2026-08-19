"""
RNFL Quadrant + Clock-hour OCR (개선판)
- 기존 parse_quadrants_cv / parse_clockhours_cv 좌표·임계값 재조정
- 다단계 contrast / conf 완화 / 매칭 허용범위 확대
- 결과: rnfl_detail.csv 저장 (patient_id, eye, vf_date, oct_date, 16개 수치, cv_flags, notes)
- 결측은 빈 셀로 두고 사용자가 수동 보완
"""
import csv, os, re, math, sys, time
import pytesseract
from PIL import Image, ImageEnhance, ImageOps

import os as _os
from pathlib import Path as _Path
# 저장소 루트는 이 파일 위치에서 유도한다(개인 경로 하드코딩 금지).
# 다른 위치에서 쓰려면 HVF_ROOT 환경변수로 덮어쓴다.
_REPO_ROOT = _os.environ.get('HVF_ROOT', str(_Path(__file__).resolve().parents[2]))


sys.stdout.reconfigure(encoding='utf-8', errors='replace')
pytesseract.pytesseract.tesseract_cmd = r'C:\Program Files\Tesseract-OCR\tesseract.exe'

BASE = _REPO_ROOT
ML_FULL = os.path.join(BASE, 'ml_dataset.csv')
ML_90D  = os.path.join(BASE, 'ml_final_90d.csv')
ML_180D = os.path.join(BASE, 'ml_final_180d.csv')
OCT     = os.path.join(BASE, 'oct_values.csv')
OUT_90  = os.path.join(BASE, 'rnfl_detail.csv')
OUT_180 = os.path.join(BASE, 'rnfl_detail_180d.csv')
DATASET_FINAL_CSV = os.path.join(BASE, 'kcc', 'dataset_final.csv')  # 원본 raw RNFL 이미지 경로(rnfl_path) 조회용

# 2026-07-23: 프로덕션 quadrants.png(cirrus_out/by_case)를 만든 생성 스크립트가
# git 이력에 없음 — 유실 확정. 앵커 텍스트("NA","95%" 범례) 위치를 원본 이미지와
# 기존 quadrants.png에서 각각 OCR로 찾아 역산한 결과:
#   x1=590/1700(=0.3471, config.yaml의 0.2000과 다름 — config 불신 원칙에 따라
#   직접 재측정한 값 채택), y1=1429/2200(=0.6495, config와 일치),
#   x2=1225/1700(=0.7206).
# y2=1630/2200(=0.7409)는 회귀 전 원래값 그대로 유지 — S/T/N은 이 크롭만 쓴다.
# (2026-07-23 1차 시도: y2를 0.758로 늘려 I를 같은 크롭에서 잡으려 했으나
#  max_d(=max(W,H)*0.10)가 커진 이미지 크기에 비례해 커지면서 greedy 최근접
#  매칭이 T/N을 오배정 — rnfl_q_s 97.8%→87.9%, q_t 96.7%→68.8%, q_n 93.3%→
#  78.2%로 회귀 확인 후 되돌림. I는 아래 RNFL_I_ONLY_CROP으로 완전히 분리.)
RNFL_QUADRANTS_CROP = (590/1700, 1429/2200, 1225/1700, 1630/2200)

# I(하부)분면 전용 크롭 — 위 크롭과 y범위가 절대 겹치지 않음(1630/2200에서 바로
# 시작). I값은 y≈0.752-0.753 실측, 시계시간 첫 줄은 y≈0.766~부터라 y2=1672/2200
# (=0.76)이 안전 마진. x1/x2는 위 크롭과 동일(OD/OS 두 원 다 포함, 폭 그대로).
RNFL_I_ONLY_CROP = (590/1700, 1630/2200, 1225/1700, 1672/2200)


def crop_quadrants_fresh(orig_rnfl_path):
    """원본 RNFL 리포트 JPG에서 quadrants(S/T/N) 영역을 재크롭 — 회귀 전 원래 경계."""
    if not orig_rnfl_path or not os.path.exists(orig_rnfl_path):
        return None
    img = Image.open(orig_rnfl_path)
    W, H = img.size
    x1, y1, x2, y2 = RNFL_QUADRANTS_CROP
    return img.crop((int(W*x1), int(H*y1), int(W*x2), int(H*y2)))


def crop_i_only_fresh(orig_rnfl_path):
    """I분면 전용 — S/T/N 크롭과 겹치지 않는 별도 영역."""
    if not orig_rnfl_path or not os.path.exists(orig_rnfl_path):
        return None
    img = Image.open(orig_rnfl_path)
    W, H = img.size
    x1, y1, x2, y2 = RNFL_I_ONLY_CROP
    return img.crop((int(W*x1), int(H*y1), int(W*x2), int(H*y2)))


def parse_i_only(img_or_path):
    """I분면 전용 크롭에서 od_I/os_I 값만 추출. S/T/N 타겟이 아예 없어 간섭 불가."""
    empty = {'od_I': None, 'os_I': None}
    if img_or_path is None:
        return empty
    img = Image.open(img_or_path).convert('RGB') if isinstance(img_or_path, str) else img_or_path.convert('RGB')
    W, H = img.width, img.height
    cands = multi_stage_candidates(img)
    # 메인 quadrants 크롭에서 쓰던 od_I/os_I 상대 x비율(0.21/0.79) 재사용,
    # y는 이 작은 크롭의 세로 중앙 부근(I 값이 크롭 상단~중앙에 위치, 실측 기반)
    targets = {'od_I': (W*0.21, H*0.45), 'os_I': (W*0.79, H*0.45)}
    max_d = max(W, H) * 0.5  # 크롭 자체가 작고 타겟이 2개뿐이라 넉넉히
    result = {}
    used = set()
    pairs = []
    for key, (tx, ty) in targets.items():
        for idx, (v, cx, cy) in enumerate(cands):
            d = math.sqrt((cx - tx) ** 2 + (cy - ty) ** 2)
            if d < max_d:
                pairs.append((d, key, idx, v))
    pairs.sort()
    for d, key, idx, v in pairs:
        if key in result or idx in used:
            continue
        result[key] = v
        used.add(idx)
    for key in targets:
        result.setdefault(key, None)
    return result


# ── OCR 후보 수집 (다단계 fallback) ──────────────────────────
_CFG = '--psm 11 --oem 1 -c tessedit_char_whitelist=0123456789'

def collect_candidates(pil_l, conf_thr=10):
    """주어진 grayscale 이미지에서 (val, cx, cy) 후보 리스트 반환"""
    data = pytesseract.image_to_data(pil_l, config=_CFG, output_type=pytesseract.Output.DICT)
    out = []
    for i, txt in enumerate(data['text']):
        m = re.search(r'\d+', txt.strip())
        if not m or data['conf'][i] < conf_thr:
            continue
        v = int(m.group())
        # 2026-07-23: quadrant(GT lo=6)/clock-hour(GT lo=4) 공용 후보 수집기라
        # 더 낮은 쪽(clock-hour) 기준으로 느슨하게 잡음 — config plausibility 참고.
        if v < 4 or v > 300:
            continue
        cx = data['left'][i] + data['width'][i] // 2
        cy = data['top'][i]  + data['height'][i] // 2
        out.append((v, cx, cy))
    return out


def multi_stage_candidates(img_rgb, scale=2):
    """다단계 contrast로 후보 수집 → 좌표 dedup (같은 위치는 1번만)"""
    big = img_rgb.resize((img_rgb.width*scale, img_rgb.height*scale), Image.LANCZOS)
    gray = big.convert('L')
    all_cands = []
    for contrast in (3.0, 2.0, 4.0, 1.5, 6.0):
        enh = ImageEnhance.Contrast(gray).enhance(contrast)
        all_cands.extend(collect_candidates(enh, conf_thr=10))
    # 반전 1회
    inv = ImageOps.invert(ImageEnhance.Contrast(gray).enhance(3.0))
    all_cands.extend(collect_candidates(inv, conf_thr=10))

    # 좌표 dedup (15px 이내 같은 값이면 중복)
    dedup = []
    for v, cx, cy in all_cands:
        is_dup = False
        for v2, cx2, cy2 in dedup:
            if v == v2 and abs(cx-cx2) < 15 and abs(cy-cy2) < 15:
                is_dup = True
                break
        if not is_dup:
            dedup.append((v, cx, cy))
    # scale 보정 → 원본 좌표계로
    return [(v, cx//scale, cy//scale) for v, cx, cy in dedup]


# ── 1. Quadrants (4개 × 2안 = 8개) ──────────────────────────
def parse_quadrants(path_or_img):
    keys = ['S', 'T', 'I', 'N']
    empty = {f'{e}_{k}': None for e in ('od', 'os') for k in keys}
    if path_or_img is None:
        return empty
    if isinstance(path_or_img, str):
        if not path_or_img or not os.path.exists(path_or_img):
            return empty
        img = Image.open(path_or_img).convert('RGB')
    else:
        img = path_or_img.convert('RGB')
    W, H = img.width, img.height
    cands = multi_stage_candidates(img)

    # OD/OS 원의 중심 (대략 W*0.21 / W*0.79, H*0.5)
    # quadrants.png에서 수치 위치를 더 넓게 잡음
    targets = {
        'od_S': (W*0.21,  H*0.13),
        'od_T': (W*0.040, H*0.55),
        'od_I': (W*0.21,  H*0.92),
        'od_N': (W*0.378, H*0.55),
        'os_S': (W*0.79,  H*0.13),
        'os_N': (W*0.660, H*0.55),
        'os_I': (W*0.79,  H*0.92),
        'os_T': (W*0.985, H*0.55),
    }

    # 거리 매칭 (허용범위 W의 8% 확대 = ~50px)
    max_d = max(W, H) * 0.10
    result = {}
    used = set()
    # 가까운 순으로 그리디 할당
    pairs = []
    for key, (tx, ty) in targets.items():
        for idx, (v, cx, cy) in enumerate(cands):
            d = math.sqrt((cx-tx)**2 + (cy-ty)**2)
            if d < max_d:
                pairs.append((d, key, idx, v))
    pairs.sort()
    for d, key, idx, v in pairs:
        if key in result or idx in used:
            continue
        result[key] = v
        used.add(idx)
    for key in targets:
        result.setdefault(key, None)
    return result


# ── 2. Clock-hours (12개 × 2안 = 24개) ──────────────────────
def parse_clockhours(path):
    empty = {f'{e}_h{n:02d}': None for e in ('od', 'os') for n in range(1, 13)}
    if not path or not os.path.exists(path):
        return empty
    img = Image.open(path).convert('RGB')
    W, H = img.width, img.height
    cands = multi_stage_candidates(img)

    cx_od = W * 0.206
    cx_os = W * 0.817
    cy_c  = H * 0.522
    r_min = min(W, H) * 0.10
    r_max = min(W, H) * 0.55

    cands_od, cands_os = [], []
    for v, cx, cy in cands:
        d_od = math.sqrt((cx-cx_od)**2 + (cy-cy_c)**2)
        d_os = math.sqrt((cx-cx_os)**2 + (cy-cy_c)**2)
        if d_od < d_os:
            if r_min <= d_od <= r_max:
                ang = math.degrees(math.atan2(-(cy-cy_c), cx-cx_od))
                cands_od.append((v, ang, d_od))
        else:
            if r_min <= d_os <= r_max:
                ang = math.degrees(math.atan2(-(cy-cy_c), cx-cx_os))
                cands_os.append((v, ang, d_os))

    def assign(cands_eye):
        target = {f'h{n:02d}': 90 - n*30 for n in range(1, 13)}
        # 12시(h12)는 90도, 3시(h03)는 0도, 6시(h06)는 -90도, 9시(h09)는 180/-180도
        # 위 target은 h01=60, h02=30, h03=0, h04=-30, h05=-60, h06=-90, h07=-120(=240), h08=-150, h09=180, h10=150, h11=120, h12=90
        max_d = 25  # 각도 허용범위 (30도 슬라이스 ÷ 2 + 여유)
        pairs = []
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

    od_h = assign(cands_od)
    os_h = assign(cands_os)
    return {**{f'od_{k}': v for k, v in od_h.items()},
            **{f'os_{k}': v for k, v in os_h.items()}}


# ── 3. CV 룰 ─────────────────────────────────────────────────
def cross_validate(eye, quad, clk, avg_rnfl):
    flags = []
    qv = [quad[f'{eye}_{d}'] for d in 'STIN']
    cv = [clk[f'{eye}_h{n:02d}'] for n in range(1, 13)]
    qv_ok = [x for x in qv if x is not None]
    cv_ok = [x for x in cv if x is not None]

    if avg_rnfl and len(qv_ok) >= 2:
        m = sum(qv_ok)/len(qv_ok)
        if abs(m - avg_rnfl) > 10:
            flags.append(f'CV_quad({len(qv_ok)}of4): mean={m:.1f} vs avg_rnfl={avg_rnfl} (diff={abs(m-avg_rnfl):.1f})')

    if avg_rnfl and len(cv_ok) >= 6:
        m = sum(cv_ok)/len(cv_ok)
        if abs(m - avg_rnfl) > 10:
            flags.append(f'CV_clk({len(cv_ok)}of12): mean={m:.1f} vs avg_rnfl={avg_rnfl} (diff={abs(m-avg_rnfl):.1f})')

    # CV3: S quadrant ≈ mean(h10, h11, h12, h01, h02) (위쪽 5시간)
    if quad[f'{eye}_S']:
        s_clk = [clk[f'{eye}_h{n:02d}'] for n in (10, 11, 12, 1, 2)]
        s_clk_ok = [x for x in s_clk if x is not None]
        if len(s_clk_ok) >= 3:
            m = sum(s_clk_ok)/len(s_clk_ok)
            if abs(m - quad[f'{eye}_S']) > 15:
                flags.append(f'CV_S: quad_S={quad[f"{eye}_S"]} vs clk_top_mean={m:.1f}')

    if quad[f'{eye}_I']:
        i_clk = [clk[f'{eye}_h{n:02d}'] for n in (4, 5, 6, 7, 8)]
        i_clk_ok = [x for x in i_clk if x is not None]
        if len(i_clk_ok) >= 3:
            m = sum(i_clk_ok)/len(i_clk_ok)
            if abs(m - quad[f'{eye}_I']) > 15:
                flags.append(f'CV_I: quad_I={quad[f"{eye}_I"]} vs clk_bot_mean={m:.1f}')

    return flags


# ── 4. Run ───────────────────────────────────────────────────
def run(pairs_csv, mode='verify', out_csv=None):
    """mode: 'verify' (샘플 10개 출력) or 'full' (전체 csv 저장)"""
    ml_full = list(csv.DictReader(open(ML_FULL, encoding='utf-8-sig')))
    rnfl_map = {(r['patient_id'], r['eye'], r['vf_date']): r['rnfl_dir'] for r in ml_full}
    rnfl_path_map = {(r['patient_id'], r['eye'], r['vf_date']): r['rnfl_path']
                      for r in csv.DictReader(open(DATASET_FINAL_CSV, encoding='utf-8-sig'))}

    oct_idx = {}
    for r in csv.DictReader(open(OCT, encoding='utf-8-sig')):
        oct_idx[(r['patient_id'], r['eye'])] = r

    pair_rows = list(csv.DictReader(open(pairs_csv, encoding='utf-8-sig')))

    if mode == 'verify':
        pair_rows = pair_rows[:10]

    QKEYS = ['S', 'T', 'I', 'N']
    HKEYS = [f'h{n:02d}' for n in range(1, 13)]
    FIELDS = (
        ['patient_id', 'eye', 'vf_date', 'oct_date'] +
        [f'rnfl_q_{d.lower()}' for d in QKEYS] +
        [f'rnfl_{h}' for h in HKEYS] +
        ['cv_flags', 'notes']
    )

    # path별 캐시 (양안 통합 리포트라 OD에서 OS 결과 같이 추출됨)
    path_cache = {}

    out_rows = []
    t0 = time.time()
    for i, r in enumerate(pair_rows, 1):
        pid, eye, vfd, octd = r['patient_id'], r['eye'], r['vf_date'], r['oct_date']
        rnfl_dir = rnfl_map.get((pid, eye, vfd), '')
        notes = []

        if not rnfl_dir:
            notes.append('no_rnfl_dir')
            row = {k: '' for k in FIELDS}
            row.update({'patient_id': pid, 'eye': eye, 'vf_date': vfd, 'oct_date': octd,
                        'notes': '|'.join(notes)})
            out_rows.append(row)
            continue

        if rnfl_dir not in path_cache:
            clk_path  = os.path.join(rnfl_dir, 'clockhours.png')
            t1 = time.time()
            try:
                # S/T/N: 정적 quadrants.png(cirrus_out/by_case) 그대로 사용 —
                # 이게 진짜 "회귀 전 원래 상태"임(2026-07-23 검증: 원본에서 재크롭
                # 시 10샘플 S/T/N 73.7/26.3/63.2%로 정적파일 100/89.5/89.5%보다
                # 뚜렷이 나쁨 — 정적 파일이 진짜 기준, 손대지 않음).
                q = parse_quadrants(os.path.join(rnfl_dir, 'quadrants.png'))

                # I(하부)분면만 원본 이미지에서 별도 크롭+전용함수로 추가 추출해
                # 병합(2026-07-23, 옵션2 — S/T/N 경로와 완전히 분리, 간섭 불가).
                orig_path = rnfl_path_map.get((pid, eye, vfd))
                i_crop = crop_i_only_fresh(orig_path) if orig_path else None
                if i_crop is not None:
                    i_vals = parse_i_only(i_crop)
                    q['od_I'] = i_vals['od_I']
                    q['os_I'] = i_vals['os_I']

                c = parse_clockhours(clk_path)
                path_cache[rnfl_dir] = (q, c, None)
            except Exception as e:
                path_cache[rnfl_dir] = (None, None, str(e))
            if mode == 'verify':
                print(f'[{i}] {pid} {eye} ({time.time()-t1:.1f}s)')

        q, c, err = path_cache[rnfl_dir]
        if err:
            notes.append(f'ocr_error:{err[:40]}')
            row = {k: '' for k in FIELDS}
            row.update({'patient_id': pid, 'eye': eye, 'vf_date': vfd, 'oct_date': octd,
                        'notes': '|'.join(notes)})
            out_rows.append(row)
            continue

        eye_l = eye.lower()
        avg = None
        oct_r = oct_idx.get((pid, eye), {})
        try:
            avg = float(oct_r.get(f'{eye_l}_avg_rnfl', '') or 0) or None
        except: pass

        cv_flags = cross_validate(eye_l, q, c, avg)

        # 누락 카운트
        n_q_miss = sum(1 for d in QKEYS if q[f'{eye_l}_{d}'] is None)
        n_c_miss = sum(1 for h in HKEYS if c[f'{eye_l}_{h}'] is None)
        if n_q_miss >= 2: notes.append(f'quad_miss={n_q_miss}/4')
        if n_c_miss >= 3: notes.append(f'clk_miss={n_c_miss}/12')

        row = {
            'patient_id': pid, 'eye': eye, 'vf_date': vfd, 'oct_date': octd,
            'cv_flags': '|'.join(cv_flags),
            'notes': '|'.join(notes),
        }
        for d in QKEYS:
            row[f'rnfl_q_{d.lower()}'] = q[f'{eye_l}_{d}'] if q[f'{eye_l}_{d}'] is not None else ''
        for h in HKEYS:
            row[f'rnfl_{h}'] = c[f'{eye_l}_{h}'] if c[f'{eye_l}_{h}'] is not None else ''
        out_rows.append(row)

    if mode == 'full':
        target = out_csv or OUT_90
        with open(target, 'w', newline='', encoding='utf-8-sig') as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(out_rows)
        print(f'\n저장: {target} ({len(out_rows)}행)')

    # 요약
    print(f'\n=== {mode} 요약 ({len(out_rows)} 페어, {time.time()-t0:.1f}s) ===')
    q_total = c_total = 0
    q_filled = c_filled = 0
    n_cv_flagged = 0
    n_notes = 0
    for r in out_rows:
        for d in 'STIN':
            q_total += 1
            if r[f'rnfl_q_{d.lower()}'] != '':
                q_filled += 1
        for n in range(1, 13):
            c_total += 1
            if r[f'rnfl_h{n:02d}'] != '':
                c_filled += 1
        if r['cv_flags']: n_cv_flagged += 1
        if r['notes']:    n_notes += 1
    print(f'Quadrant 충원율  : {q_filled}/{q_total} ({q_filled/q_total*100:.1f}%)')
    print(f'Clock-hour 충원율: {c_filled}/{c_total} ({c_filled/c_total*100:.1f}%)')
    print(f'CV flag 발생     : {n_cv_flagged}/{len(out_rows)}')
    print(f'notes 발생       : {n_notes}/{len(out_rows)}')
    return out_rows


if __name__ == '__main__':
    mode  = sys.argv[1] if len(sys.argv) > 1 else 'verify'
    which = sys.argv[2] if len(sys.argv) > 2 else '90d'
    if which == '180d':
        run(ML_180D, mode, OUT_180)
    else:
        run(ML_90D, mode, OUT_90)
