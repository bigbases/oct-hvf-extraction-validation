"""
scripts/extract_demographics.py
GCA 리포트 헤더 OCR → data/demographics.csv

dob: re.findall 기반 — 검사날짜(2026)가 먼저 나오는 버그 수정 (1920-2010 필터)
sex: Tesseract kor+eng로 한국어 추출 불가 → oct_inventory_merged 기존값만 사용

안전장치:
  - 기존 dob 있는 환자(105명): dob 무변경 (existing or ocr로 기존 우선)
  - dob 결측 환자(39명)만 헤더 OCR
  - 회귀 10개 OK 샘플: OCR dob vs 기존값 비교
  - sex: Unknown/''/이외 → '' (분석 시 결측 처리)
"""
import csv, json, re, datetime, sys, os
from pathlib import Path
from collections import defaultdict

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))

import yaml
import pytesseract
from PIL import Image

cfg = yaml.safe_load(open(_ROOT / 'config' / 'params.yaml', encoding='utf-8'))
pytesseract.pytesseract.tesseract_cmd = cfg['paths']['tesseract']

_OCR_CFG = '--psm 6'
_LANG    = 'kor+eng'
_HDR_Y   = 0.20
_HDR_X   = 0.60
_SCALE   = 2

# re.findall: 전체 순회 후 1920-2010만 취득 (검사날짜 2026 버그 수정)
_DOB_RE = re.compile(r'(\d{4})\s*[-./]\s*(\d{1,2})\s*[-./]\s*(\d{1,2})')

_SEX_MAP = {
    'Male': 'Male', 'Female': 'Female',
    'male': 'Male', 'female': 'Female',
    '남성': 'Male', '여성': 'Female',
    '남': 'Male',   '여': 'Female',
}
_VALID_SEX = {'Male', 'Female'}


def _read_csv(path):
    return list(csv.DictReader(open(path, encoding='utf-8-sig')))


def _ocr_header(img_path):
    img = Image.open(img_path).convert('RGB')
    W, H = img.width, img.height
    hdr = img.crop((0, 0, int(W * _HDR_X), int(H * _HDR_Y)))
    big = hdr.resize((hdr.width * _SCALE, hdr.height * _SCALE), Image.LANCZOS)
    return pytesseract.image_to_string(big, config=_OCR_CFG, lang=_LANG)


def _parse_dob(text):
    for y_s, mo_s, d_s in _DOB_RE.findall(text):
        y, mo, d = int(y_s), int(mo_s), int(d_s)
        if 1920 <= y <= 2010 and 1 <= mo <= 12 and 1 <= d <= 31:
            return f'{y:04d}-{mo:02d}-{d:02d}'
    return None


def _norm_sex(s):
    s = (s or '').strip()
    return _SEX_MAP.get(s) or _SEX_MAP.get(s.lower()) or ''


def main():
    from hvf.config import data_dir, results_dir
    from hvf.registry import sha256_file

    data = data_dir()
    res  = results_dir()

    vf_rows  = _read_csv(_ROOT / 'data' / 'vf_results.csv')
    inv_rows = _read_csv(_ROOT / 'oct_inventory_merged.csv')
    ds_rows  = _read_csv(_ROOT / 'kcc' / 'dataset_final.csv')

    cohort_pats = sorted({r['patient_id'] for r in vf_rows})

    # patient별 기존 demo
    inv_map = {}
    for r in inv_rows:
        pid = r['patient_id']
        if pid not in inv_map:
            inv_map[pid] = {
                'dob':  r.get('dob', '').strip(),
                'sex':  _norm_sex(r.get('sex', '')),
                'conf': r.get('ocr_confidence', '').strip(),
                'src':  r.get('source', '').strip(),
            }

    # patient별 GCA 이미지
    gca_map = {}
    for r in ds_rows:
        pid = r['patient_id']
        if pid not in gca_map:
            gca = r.get('gca_path', '')
            if gca and os.path.exists(gca):
                gca_map[pid] = gca

    # ---- 회귀 검증: 10개 OK 샘플 (dob만) ----
    ok_pids = [
        pid for pid in cohort_pats
        if inv_map.get(pid, {}).get('conf') == 'OK'
        and inv_map.get(pid, {}).get('dob')
        and pid in gca_map
    ][:10]

    print('=== 회귀 검증: 10개 OK 샘플 dob OCR vs 기존값 ===')
    reg_pass = 0
    for pid in ok_pids:
        text    = _ocr_header(gca_map[pid])
        ocr_dob = _parse_dob(text)
        exist   = inv_map[pid]['dob']
        match   = (ocr_dob == exist)
        if match:
            reg_pass += 1
    print(f'  dob 회귀: {reg_pass}/10 일치')
    print()

    # ---- dob-missing 39명만 OCR ----
    dob_missing = [pid for pid in cohort_pats if not inv_map.get(pid, {}).get('dob')]
    print(f'dob OCR 대상: {len(dob_missing)}명')

    by_src = defaultdict(list)   # 'dw' / 'dw+zeiss' / 'no_image'
    dob_ocr_out = {}             # pid → extracted dob (or None)

    for pid in dob_missing:
        inv_r    = inv_map.get(pid, {})
        src_type = 'dw' if inv_r.get('src') == 'dw' else 'dw+zeiss'

        if pid not in gca_map:
            by_src['no_image'].append({'pid': pid, 'status': 'no_image'})
            dob_ocr_out[pid] = None
            continue

        text    = _ocr_header(gca_map[pid])
        ocr_dob = _parse_dob(text)
        status  = 'recovered' if ocr_dob else 'fail'
        by_src[src_type].append({'pid': pid, 'status': status, 'ocr_dob': ocr_dob})
        dob_ocr_out[pid] = ocr_dob

    # ---- 유형별 복구율 ----
    print()
    print('=== 유형별 dob 복구율 ===')
    for src_type in ('dw', 'dw+zeiss', 'no_image'):
        results = by_src[src_type]
        if not results:
            continue
        n_rec  = sum(1 for r in results if r['status'] == 'recovered')
        n_fail = len(results) - n_rec
        print(f'  {src_type} ({len(results)}명): 복구={n_rec}  실패={n_fail}')
        for r in results:
            if r['status'] == 'fail':
                print(f'    fail: pid={r["pid"]}')

    # ---- demographics.csv 생성 ----
    # sex: 기존 oct_inventory_merged 값만 사용 (Unknown → '' 로 정규화)
    # note: sex OCR 시도했으나 Tesseract kor+eng에서 한국어 텍스트 추출 불가
    out_rows = []
    for pid in cohort_pats:
        inv_r = inv_map.get(pid, {})
        existing_dob = inv_r.get('dob', '')
        existing_sex = inv_r.get('sex', '')  # 이미 _norm_sex 적용됨

        # dob: 기존 → OCR 순
        final_dob = existing_dob or dob_ocr_out.get(pid) or ''
        # sex: 기존값만 (Unknown/'' → incomplete로 flag, 변경 없음)
        final_sex = existing_sex if existing_sex in _VALID_SEX else ''

        demo_source = 'oct_inventory'
        if pid in dob_missing:
            demo_source = 'header_ocr' if final_dob else 'missing'

        dob_ok = bool(final_dob)
        sex_ok = final_sex in _VALID_SEX
        if dob_ok and sex_ok:
            flag = 'OK'
        elif dob_ok:
            flag = 'dob_only'
        elif sex_ok:
            flag = 'sex_only'
        else:
            flag = 'incomplete'

        out_rows.append({
            'patient_id': pid,
            'dob':         final_dob,
            'sex':         final_sex,
            'demo_source': demo_source,
            'demo_flag':   flag,
        })

    out_path = data / 'demographics.csv'
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['patient_id', 'dob', 'sex',
                                          'demo_source', 'demo_flag'])
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    # ---- 최종 요약 ----
    n_both     = sum(1 for r in out_rows if r['demo_flag'] == 'OK')
    n_dob_only = sum(1 for r in out_rows if r['demo_flag'] == 'dob_only')
    n_sex_only = sum(1 for r in out_rows if r['demo_flag'] == 'sex_only')
    n_neither  = sum(1 for r in out_rows if r['demo_flag'] == 'incomplete')

    print()
    print('=== 최종 demographics.csv (144명) ===')
    print(f'  완전 (dob+sex): {n_both}명')
    print(f'  dob only:       {n_dob_only}  (sex not recoverable - Korean OCR limit)')
    print(f'  sex만:          {n_sex_only}명')
    print(f'  둘다 없음:      {n_neither}명')
    print()
    # 눈 수 환산 (평균 276/144=1.917눈/환자)
    # 실제 눈수는 vf_results에서 계산
    vf_rows2 = _read_csv(_ROOT / 'data' / 'vf_results.csv')
    pid_to_eyes = defaultdict(int)
    for r in vf_rows2:
        pid_to_eyes[r['patient_id']] += 1
    eyes_both     = sum(pid_to_eyes[r['patient_id']] for r in out_rows if r['demo_flag'] == 'OK')
    eyes_dob_only = sum(pid_to_eyes[r['patient_id']] for r in out_rows if r['demo_flag'] == 'dob_only')
    print(f'  공변량 분석 가능 눈수:')
    print(f'    age+sex 보정: {eyes_both}눈 / 276 ({100*eyes_both/276:.1f}%)')
    print(f'    age만 보정:   {eyes_both+eyes_dob_only}눈 / 276 ({100*(eyes_both+eyes_dob_only)/276:.1f}%)')

    incomplete_pids = [r['patient_id'] for r in out_rows if r['demo_flag'] == 'incomplete']
    if incomplete_pids:
        print(f'  둘다 없음 PIDs ({len(incomplete_pids)}명): {incomplete_pids}')
    print(f'  출력: {out_path}')

    # Provenance
    sha  = sha256_file(str(out_path))
    prov = {
        'generated': datetime.date.today().isoformat(),
        'output':    'data/demographics.csv',
        'sha256':    sha,
        'n_patients': len(out_rows),
        'stats': {
            'dob_and_sex': n_both,
            'dob_only':    n_dob_only,
            'sex_only':    n_sex_only,
            'incomplete':  n_neither,
        },
        'eyes_stats': {
            'age_sex_adjusted': eyes_both,
            'age_only_adjusted': eyes_both + eyes_dob_only,
            'total': 276,
        },
        'regression': {
            'n_samples': 10,
            'dob_match': reg_pass,
            'note': 'OK 환자 10개 샘플 dob OCR vs 기존값 비교',
        },
        'method': {
            'dob': 're.findall 전체 순회 후 1920-2010 필터 (검사날짜 버그 수정)',
            'sex': 'oct_inventory_merged 기존값만 — Korean OCR 불가 (Tesseract kor+eng 한계)',
        },
        'dob_recovery': {
            src: {
                'n': len(by_src[src]),
                'recovered': sum(1 for r in by_src[src] if r['status'] == 'recovered'),
            }
            for src in ('dw', 'dw+zeiss', 'no_image') if by_src[src]
        },
    }
    prov_path = res / 'demographics_provenance.json'
    json.dump(prov, open(prov_path, 'w', encoding='utf-8'),
              indent=2, ensure_ascii=False)
    print(f'  봉인: {prov_path}')


if __name__ == '__main__':
    main()
