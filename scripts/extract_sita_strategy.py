"""
scripts/extract_sita_strategy.py
원본 SFA(시야) 리포트 헤더에서 SITA 전략 + SW 버전을 추출한다.

배경:
  phase_b_images / cirrus_out 의 threshold_grid.png 는 비식별 크롭이라 헤더가 없다.
  전략(SITA Standard/Fast/Faster)과 버전은 원본 리포트 헤더에만 존재하므로
  Downloads\시야검사 트리의 원본 JPG를 다시 찾아 OCR 한다.

매핑:
  ml_dataset.csv 의 sfa_dir basename = "{pid}__{origstem}" 형태이고,
  origstem = 원본 리포트 파일명(stem). Downloads\시야검사 전체를 os.walk 하여
  stem->path 인덱스를 만든 뒤 (새 폴더 등 하위 포함) origstem 으로 원본을 찾는다.

출력:
  A) results/sita_strategy_distribution.json  — 집계(비 PHI)
  B) data/sita_per_eye.csv                    — 눈별 라벨(PHI: patient_id/date 포함,
                                                data/ 는 .gitignore 로 커밋 차단됨)

재현성/보안:
  - 입력 해시(ml_dataset.csv)·리포트 개수 로그
  - PHI(환자ID) 콘솔 출력 금지 — 진행률/집계만 출력
"""
import csv, os, re, sys, io, json, hashlib
from pathlib import Path
from collections import Counter

import pytesseract
from PIL import Image

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))
from hvf.config import configure_tesseract

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
configure_tesseract()

ROOT    = str(_ROOT)
# 원본(비크롭) 리포트 아카이브 — IRB 제한 원본 이미지라 레포에 포함되지 않음.
# 로컬 환경변수로 지정(미설정 시 아래 기본값은 예시 플레이스홀더).
DLROOT  = os.environ.get('HVF_RAW_VF_ARCHIVE', str(_ROOT / 'raw_vf_archive_not_included'))
ML_CSV  = os.path.join(ROOT, 'ml_dataset.csv')
MASTER  = os.path.join(ROOT, 'data', 'analysis_master.csv')
OUT_DIST = os.path.join(ROOT, 'results', 'sita_strategy_distribution.json')
OUT_PER  = os.path.join(ROOT, 'data', 'sita_per_eye.csv')


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def build_stem_index(root):
    """Downloads\시야검사 전체 트리에서 *.jpg stem -> full path (마지막 승자)."""
    idx = {}
    for dp, _dn, fn in os.walk(root):
        for f in fn:
            if f.lower().endswith('.jpg'):
                idx[f[:-4]] = os.path.join(dp, f)
    return idx


def build_origstem_map(ml_csv):
    """(patient_id, eye, vf_date) -> origstem  (ml_dataset sfa_dir basename 파싱)."""
    mp = {}
    for r in csv.DictReader(open(ml_csv, encoding='utf-8-sig')):
        sd = r.get('sfa_dir', '')
        if sd and '__' in os.path.basename(sd):
            _pid_part, rest = os.path.basename(sd).split('__', 1)
            mp[(r['patient_id'], r['eye'], r['vf_date'])] = rest
    return mp


def extract_header(img_path):
    """상단 42% + 하단 12%(2배 확대) OCR → (strategy, all_strats, version, hfa)."""
    im = Image.open(img_path).convert('RGB')
    W, H = im.size
    top = im.crop((0, 0, W, int(H * 0.42)))
    bot = im.crop((0, int(H * 0.88), W, H))
    bot = bot.resize((bot.width * 2, bot.height * 2), Image.LANCZOS)
    txt_top = pytesseract.image_to_string(top)
    txt_bot = pytesseract.image_to_string(bot)
    txt_all = txt_top + '\n' + txt_bot

    def _norm(s):
        s = s.lower()
        if 'stand' in s:  return 'SITA Standard'
        if 'faster' in s: return 'SITA Faster'
        if 'fast' in s:   return 'SITA Fast'
        return None

    strat = None
    m = re.search(r'SITA\s*(Stand?ard|Fast(?:er)?)', txt_all, re.I)
    if m:
        strat = _norm(m.group(1))
    all_strats = [_norm(s) for s in re.findall(r'SITA\s*(Stand?ard|Fast(?:er)?)', txt_all, re.I)]
    all_strats = [s for s in all_strats if s]

    ver = None
    mv = re.search(r'Version\s*([\d.]+)', txt_all, re.I)
    if mv:
        ver = mv.group(1)

    hfa = None
    mh = re.search(r'(HFA\s*3|HFA3|840-\d+)', txt_all, re.I)
    if mh:
        hfa = re.sub(r'\s+', '', mh.group(1))

    return strat, all_strats, ver, hfa


def main():
    print('=' * 60)
    print('SITA 전략 + SW 버전 추출')
    print('=' * 60)
    print(f'입력 ml_dataset SHA256: {sha256_file(ML_CSV)}')

    stem_index = build_stem_index(DLROOT)
    print(f'원본 리포트 stem 인덱스: {len(stem_index)}개 (.jpg)')

    origmap = build_origstem_map(ML_CSV)
    print(f'ml_dataset sfa_dir 매핑: {len(origmap)}개')

    cohort = list(csv.DictReader(open(MASTER, encoding='utf-8')))
    print(f'코호트: {len(cohort)}안')

    strat_counter = Counter()
    ver_counter   = Counter()
    hfa_counter   = Counter()
    multi_strat   = 0
    miss_map      = 0
    miss_file     = 0
    unreadable    = 0
    per_eye_rows  = []

    for i, r in enumerate(cohort):
        pid, eye, date = r['patient_id'], r['eye'], r['vf_date']
        key = (pid, eye, date)
        origstem = origmap.get(key)

        label = 'UNDETERMINED'
        if not origstem:
            miss_map += 1
        else:
            fpath = stem_index.get(origstem)
            if not fpath:
                miss_file += 1
            else:
                try:
                    strat, all_strats, ver, hfa = extract_header(fpath)
                    if len(set(all_strats)) > 1:
                        multi_strat += 1
                    if strat:
                        label = strat
                        strat_counter[strat] += 1
                    else:
                        unreadable += 1
                    if ver: ver_counter[ver] += 1
                    if hfa: hfa_counter[hfa] += 1
                except Exception:
                    miss_file += 1

        per_eye_rows.append({'patient_id': pid, 'eye': eye,
                             'vf_date': date, 'sita_strategy': label})

        if (i + 1) % 50 == 0:
            print(f'  진행 {i + 1}/{len(cohort)} ...', flush=True)

    n_mapped = len(cohort) - miss_map - miss_file

    # ── 눈별 CSV (PHI → data/) ────────────────────────────────
    with open(OUT_PER, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=['patient_id', 'eye', 'vf_date', 'sita_strategy'])
        w.writeheader()
        w.writerows(per_eye_rows)

    # ── 집계 JSON (비 PHI → results/) ─────────────────────────
    dist = {
        'generated': '2026-07-09',
        'input_sha256': {'ml_dataset': sha256_file(ML_CSV),
                         'analysis_master': sha256_file(MASTER)},
        'n_cohort': len(cohort),
        'n_readable': n_mapped - unreadable,
        'sita_strategy_distribution': dict(strat_counter),
        'unreadable': unreadable,
        'file_not_found': miss_file,
        'not_mapped': miss_map,
        'multi_strategy_reports': multi_strat,
        'sw_version_distribution': dict(ver_counter),
        'hfa_model_distribution': dict(hfa_counter),
        'method': ('OCR of original SFA report header (top 42% + bottom 12% x2); '
                   'mapped via ml_dataset sfa_dir basename origstem; '
                   'full 시야검사 tree os.walk index'),
        'per_eye_file': 'data/sita_per_eye.csv (PHI, gitignored)',
    }
    json.dump(dist, open(OUT_DIST, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)

    # ── 콘솔 요약 (비 PHI) ────────────────────────────────────
    print()
    print(f'SITA 전략 분포 (판독 {n_mapped - unreadable}/{len(cohort)}):')
    for k, v in strat_counter.most_common():
        print(f'  {k:16s}: {v}')
    print(f'  {"UNREADABLE":16s}: {unreadable}')
    print(f'  {"FILE NOT FOUND":16s}: {miss_file}')
    print(f'  {"NOT MAPPED":16s}: {miss_map}')
    print(f'  (복수 전략 표기 리포트: {multi_strat})')
    print()
    print('SW 버전 분포:')
    for k, v in ver_counter.most_common():
        print(f'  Version {k}: {v}')
    print('HFA 모델 분포:')
    for k, v in hfa_counter.most_common():
        print(f'  {k}: {v}')
    print()
    print(f'저장 A(집계): {OUT_DIST}')
    print(f'저장 B(눈별): {OUT_PER}  [PHI, gitignored]')


if __name__ == '__main__':
    main()
