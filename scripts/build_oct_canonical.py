"""
STEP 4 — oct_canonical.csv 빌드

입력:
  oct_tabular_90d.csv       (루트 — GCL+RNFL 수동검수본, primary)
  rnfl_detail.csv           (루트 — RNFL fallback, 4눈용)
  kcc/dataset_final.csv     (코호트 인덱스 276눈)
  data/signal_strength.csv  (SS, patient+oct_date 키)

출력:
  data/oct_canonical.csv    (276 rows)
  results/step4_provenance.json (SHA256 봉인)

규칙:
  - GCL/RNFL 소스: oct_tabular_90d (수동검수, OS T-N swap 수정, inferior 복원, q_i 복원)
  - RNFL fallback: rnfl_detail.csv (oct_tabular_90d에 없는 4눈)
  - GCL avg/min BVAL 미적용 (수동검수 완료), 섹터 BVAL>=200 → None (극단값 안전장치)
  - RNFL BVAL>=200 → None
  - rnfl_q_i 포함 (oct_tabular_90d에 정상값 확인)
  - 코호트 276눈만 (nearest-1, GCA+RNFL ±90d)
  - SS: patient+oct_date 키 (eye 없음)
"""
import csv, json, hashlib, datetime, sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(_ROOT / 'src'))

from hvf.config import data_dir, results_dir
from hvf.registry import sha256_file

BVAL_THR = 200

GCL_SECTOR_COLS = [
    'od_s_sup','od_s_sup_t','od_s_inf_t','od_s_inf','od_s_inf_n','od_s_sup_n',
    'os_s_sup','os_s_sup_t','os_s_inf_t','os_s_inf','os_s_inf_n','os_s_sup_n',
]
RNFL_COLS = [
    'rnfl_q_s','rnfl_q_t','rnfl_q_i','rnfl_q_n',   # q_i 복원 (수동검수본에 정상값)
    'rnfl_h01','rnfl_h02','rnfl_h03','rnfl_h04','rnfl_h05','rnfl_h06',
    'rnfl_h07','rnfl_h08','rnfl_h09','rnfl_h10','rnfl_h11','rnfl_h12',
]
SS_COLS = ['ss_gca_od','ss_gca_os','ss_rnfl_od','ss_rnfl_os']

OUT_COLS = (
    ['patient_id','eye','oct_date'] +
    ['avg_gcl_od','avg_gcl_os','min_gcl_od','min_gcl_os'] +
    GCL_SECTOR_COLS +
    RNFL_COLS +
    SS_COLS
)


def _read_csv(path):
    return list(csv.DictReader(open(path, encoding='utf-8-sig')))


def _to_int_or_none(val, bval_thr=None):
    if val in ('', None, 'None', 'nan'):
        return None
    try:
        v = float(val)
        if bval_thr is not None and v >= bval_thr:
            return None
        return int(v)
    except (ValueError, TypeError):
        return None


def _to_float_or_none(val):
    if val in ('', None, 'None', 'nan'):
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def main():
    data  = data_dir()
    res   = results_dir()

    cohort  = _read_csv(_ROOT / 'kcc' / 'dataset_final.csv')
    tab_oct = _read_csv(_ROOT / 'oct_tabular_90d.csv')    # GCL+RNFL 수동검수본 (primary)
    rnfl_fb = _read_csv(_ROOT / 'rnfl_detail.csv')         # RNFL fallback (4눈용)
    ss_raw  = _read_csv(data / 'signal_strength.csv')

    # 인덱스 빌드
    oct_map  = {(r['patient_id'], r['eye'], r['oct_date']): r for r in tab_oct}
    rnfl_map = {(r['patient_id'], r['eye'], r['oct_date']): r for r in rnfl_fb}
    ss_map   = {(r['patient_id'], r['oct_date']): r for r in ss_raw}

    # 276 코호트 키: vf_results.csv anchor → dataset_final에서 oct_date 취득
    # (일부 vf_date가 date-remap으로 불일치 → gap_days로 중복 해소)
    from collections import defaultdict
    vf_rows  = list(csv.DictReader(open(_ROOT / 'data' / 'vf_results.csv', encoding='utf-8-sig')))
    vf_map   = {(r['patient_id'], r['eye']): r for r in vf_rows}
    vf_eyes  = set(vf_map.keys())

    cands = defaultdict(list)
    for r in cohort:
        k = (r['patient_id'], r['eye'])
        if k in vf_eyes and r.get('has_rnfl','0') == '1':
            cands[k].append(r)

    cohort_keys = []
    for k, rows in cands.items():
        vf_r = vf_map[k]
        if len(rows) == 1:
            best = rows[0]
        else:
            match = [r for r in rows if r.get('gap_days','') == vf_r.get('gap_days','')]
            if len(match) == 1:
                best = match[0]
            else:
                try:
                    best = min(rows, key=lambda r: abs(int(r.get('gap_days', 999))))
                except Exception:
                    best = rows[0]
        cohort_keys.append((best['patient_id'], best['eye'], best['oct_date']))

    print(f'코호트 키: {len(cohort_keys)}개')

    rows_out = []
    missing_oct = 0
    missing_rnfl = 0

    for pid, eye, oct_date in cohort_keys:
        row = {'patient_id': pid, 'eye': eye, 'oct_date': oct_date}

        # GCL (avg/min/섹터) — oct_tabular_90d 수동검수본
        oct_r = oct_map.get((pid, eye, oct_date), {})
        for col in ['avg_gcl_od','avg_gcl_os','min_gcl_od','min_gcl_os']:
            row[col] = _to_int_or_none(oct_r.get(col,''))
        for col in GCL_SECTOR_COLS:
            row[col] = _to_int_or_none(oct_r.get(col,''), bval_thr=BVAL_THR)
        if not oct_r:
            missing_oct += 1

        # RNFL (q_s/q_t/q_i/q_n + h01-h12) — oct_tabular_90d primary, rnfl_detail fallback
        if oct_r:
            rnfl_r = oct_r           # oct_tabular_90d에 RNFL도 포함
        else:
            rnfl_r = rnfl_map.get((pid, eye, oct_date), {})
        for col in RNFL_COLS:
            row[col] = _to_int_or_none(rnfl_r.get(col,''), bval_thr=BVAL_THR)
        if not rnfl_r:
            missing_rnfl += 1

        # Signal Strength (patient+oct_date 키)
        ss_r = ss_map.get((pid, oct_date), {})
        for col in SS_COLS:
            row[col] = _to_int_or_none(ss_r.get(col,''))

        rows_out.append(row)

    out_path = data / 'oct_canonical.csv'
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)

    print(f'출력: {out_path} ({len(rows_out)}행)')
    print(f'OCT raw 미매칭: {missing_oct}')
    print(f'RNFL 미매칭: {missing_rnfl}')

    # None율 요약
    for col in ['avg_gcl_od','avg_gcl_os','min_gcl_od','min_gcl_os',
                'rnfl_q_s','rnfl_q_t','rnfl_q_n']:
        n_none = sum(1 for r in rows_out if r[col] is None)
        print(f'  {col}: None {n_none}/{len(rows_out)} ({100*n_none/len(rows_out):.1f}%)')

    # SHA-256 봉인
    sha = sha256_file(str(out_path))
    prov = {
        'generated': datetime.date.today().isoformat(),
        'output': 'data/oct_canonical.csv',
        'sha256': sha,
        'n_eyes': len(rows_out),
        'n_patients': len({r['patient_id'] for r in rows_out}),
        'sources': {
            'GCL': 'oct_tabular_90d.csv (수동검수본, OS T-N swap 수정, inferior 복원)',
            'RNFL': 'oct_tabular_90d.csv primary + rnfl_detail.csv fallback (4눈)',
            'signal_strength': 'data/signal_strength.csv',
        },
        'excluded': [
            'od_avg_rnfl, os_avg_rnfl, od_vert_cd, os_vert_cd (DL논문용 보존)',
        ],
        'restored': [
            'rnfl_q_i (수동검수본에 정상값 확인, 이전 구 파이프라인 100% None 오류 수정)',
            'od_s_inf, os_s_inf inferior 섹터 (3-stage OCR 33-36% None → 수동검수 3% None)',
            'OS 섹터 temporal↔nasal swap 수정 (ocr_oct_values.py 레이블 오류 수정)',
        ],
        'BVAL_treatment': f'GCL 섹터·RNFL >={BVAL_THR} -> None (안전장치, avg/min은 수동검수 완료)',
        'accuracy_ref': 'results/oct_extraction_accuracy.json (Phase 3 정확도 봉인)',
    }
    prov_path = res / 'step4_provenance.json'
    json.dump(prov, open(prov_path, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print(f'봉인: {prov_path}')


if __name__ == '__main__':
    main()
