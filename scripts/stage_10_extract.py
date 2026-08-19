"""
STEP 1 — 원본 이미지에서 모든 수치 재추출.

출력:
  data/sfa_thresholds_raw.csv   (340 rows, p01-p54 + test_pattern)
  data/oct_values_raw.csv       (340 rows, GCL/RNFL/sector/ss)
  data/rnfl_detail_raw.csv      (340 rows, quadrant/clockhour)
  data/signal_strength.csv      (340 rows, ss_gca/ss_rnfl × OD/OS)

검증:
  --verify-determinism  : 첫 10행을 두 번 실행해 SHA-256 일치 확인
  --count-check         : 278 눈(±90d·GCA·RNFL 필수) 재현 여부 보고

PHI 규칙:
  로그에 patient_id·날짜 원문 출력 금지 → 진행률만 번호로 출력.
  입력 파일 SHA-256 → registry 기록.
"""
import csv, os, sys, hashlib, json, datetime, argparse, logging
from pathlib import Path

# ── PYTHONPATH: src/ 추가 (scripts/에서 직접 실행 시) ──────────
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(_ROOT / 'src'))

from hvf.config import load_config, get as cfg, data_dir, results_dir, set_seed, configure_tesseract
from hvf.registry import record_result, sha256_file
from hvf.ocr_vf   import extract_vf, POINT_LABELS
from hvf.ocr_oct  import extract_oct, FIELDS as OCT_FIELDS
from hvf.ocr_rnfl import extract_rnfl, RNFL_FIELDS

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S')
log = logging.getLogger(__name__)

# ── 출력 파일명 ───────────────────────────────────────────────
OUT_VF   = 'sfa_thresholds_raw.csv'
OUT_OCT  = 'oct_values_raw.csv'
OUT_RNFL = 'rnfl_detail_raw.csv'
OUT_SS   = 'signal_strength.csv'

VF_FIELDS = (['patient_id', 'eye', 'vf_date', 'test_pattern'] +
              POINT_LABELS + ['n_detected', 'n_missing'])

OCT_RAW_FIELDS = (
    ['patient_id', 'eye', 'oct_date'] +
    [f for f in OCT_FIELDS if f not in ('patient_id', 'eye')] +
    ['ss_gca_od', 'ss_gca_os', 'ss_rnfl_od', 'ss_rnfl_os']
)

SS_FIELDS = ['patient_id', 'oct_date', 'ss_gca_od', 'ss_gca_os', 'ss_rnfl_od', 'ss_rnfl_os']


def _load_cohort():
    """dataset_final.csv 로드 → 리스트 반환."""
    path = _ROOT / cfg('cohort', 'dataset_csv')
    return list(csv.DictReader(open(path, encoding='utf-8-sig')))


def _sha256_csv(path: Path) -> str:
    return sha256_file(str(path))


def _run_extraction(rows, *, limit=None):
    """
    rows: dataset_final 행 목록.
    limit: None이면 전체, 정수면 첫 N행만.

    반환: (vf_rows, oct_rows, rnfl_rows, ss_rows)
    """
    if limit is not None:
        rows = rows[:limit]

    vf_rows = []; oct_rows = []; rnfl_rows = []; ss_rows = []
    n = len(rows)

    for i, r in enumerate(rows, 1):
        pid     = r['patient_id']
        eye     = r['eye']
        vf_date = r['vf_date']
        oct_date= r['oct_date']
        vf_path = r['vf_path']
        gca_path= r['gca_path']
        rnfl_path=r.get('rnfl_path', '')

        if i % 20 == 0 or i == n:
            log.info('[%d/%d] 처리 중 ...', i, n)

        # ── VF ─────────────────────────────────────────────────
        try:
            vf_res = extract_vf(vf_path, eye)
        except Exception as e:
            log.warning('[%d] VF 실패: %s', i, str(e)[:60])
            vf_res = {lbl: None for lbl in POINT_LABELS}
            vf_res.update({'test_pattern': 'error', 'n_detected': 0, 'n_missing': 54})

        vf_row = {'patient_id': pid, 'eye': eye, 'vf_date': vf_date}
        vf_row.update({k: vf_res.get(k) for k in POINT_LABELS + ['test_pattern', 'n_detected', 'n_missing']})
        vf_rows.append(vf_row)

        # ── OCT ────────────────────────────────────────────────
        try:
            oct_row = extract_oct(gca_path, rnfl_path, pid, eye)
        except Exception as e:
            log.warning('[%d] OCT 실패: %s', i, str(e)[:60])
            oct_row = {'patient_id': pid, 'eye': eye, 'cv_flags': f'error:{str(e)[:40]}',
                       'ss_gca_od': None, 'ss_gca_os': None, 'ss_rnfl_od': None, 'ss_rnfl_os': None}

        oct_row['oct_date'] = oct_date
        oct_rows.append(oct_row)

        # ── RNFL detail ────────────────────────────────────────
        avg_od = oct_row.get('od_avg_rnfl')
        avg_os = oct_row.get('os_avg_rnfl')
        try:
            avg_od_f = float(avg_od) if avg_od is not None else None
            avg_os_f = float(avg_os) if avg_os is not None else None
        except (TypeError, ValueError):
            avg_od_f = avg_os_f = None

        try:
            rnfl_row = extract_rnfl(rnfl_path, avg_od_f, avg_os_f,
                                    pid, eye, vf_date, oct_date)
        except Exception as e:
            log.warning('[%d] RNFL 실패: %s', i, str(e)[:60])
            rnfl_row = {k: '' for k in RNFL_FIELDS}
            rnfl_row.update({'patient_id': pid, 'eye': eye,
                             'vf_date': vf_date, 'oct_date': oct_date,
                             'notes': f'error:{str(e)[:40]}'})
        rnfl_rows.append(rnfl_row)

        # ── signal_strength 분리 행 ─────────────────────────────
        ss_rows.append({
            'patient_id': pid,
            'oct_date':   oct_date,
            'ss_gca_od':  oct_row.get('ss_gca_od'),
            'ss_gca_os':  oct_row.get('ss_gca_os'),
            'ss_rnfl_od': oct_row.get('ss_rnfl_od'),
            'ss_rnfl_os': oct_row.get('ss_rnfl_os'),
        })

    return vf_rows, oct_rows, rnfl_rows, ss_rows


def _write_csv(path: Path, fields: list, rows: list):
    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def _cohort_count_check(vf_rows):
    """
    278 눈 재현 검증: 24-2 필터 → join dataset_final → ±90d → nearest-1.
    Fail-loud: 278이 아니면 RuntimeError.
    """
    ds_path = _ROOT / cfg('cohort', 'dataset_csv')
    ds_map  = {}
    for r in csv.DictReader(open(ds_path, encoding='utf-8-sig')):
        try:
            g = int(r['gap_days'])
        except (ValueError, TypeError):
            g = 9999
        has_gca  = r.get('has_gca', '').strip()  in ('1', 'Y', 'True', 'true')
        has_rnfl = r.get('has_rnfl', '').strip() in ('1', 'Y', 'True', 'true')
        ds_map[(r['patient_id'], r['eye'], r['vf_date'])] = {
            'gap': g, 'has_gca': has_gca, 'has_rnfl': has_rnfl
        }

    rnfl_required = cfg('constants', 'cohort', 'rnfl_required') if \
        _try_cfg('constants', 'cohort', 'rnfl_required') else True
    # 실제로는 params.yaml constants.cohort 아닌 constants 아래 cohort 없음 → 별도 키 사용
    # params.yaml에 constants.cohort.rnfl_required 없음 → 기본 True 사용
    allow_patterns = set(cfg('constants', 'cohort', 'test_pattern_include')
                         if _try_cfg('constants', 'cohort', 'test_pattern_include')
                         else ['24-2'])

    best = {}  # (pid, eye) → nearest gap record
    for vf in vf_rows:
        tp  = vf.get('test_pattern', 'unknown')
        if tp not in allow_patterns:
            continue
        key_ds = (vf['patient_id'], vf['eye'], vf['vf_date'])
        ds_r   = ds_map.get(key_ds)
        if ds_r is None or not ds_r['has_gca']:
            continue
        if rnfl_required and not ds_r['has_rnfl']:
            continue
        gap = abs(ds_r['gap'])
        if gap > 90:
            continue
        k = (vf['patient_id'], vf['eye'])
        if k not in best or gap < abs(ds_map.get((best[k]['patient_id'], best[k]['eye'], best[k]['vf_date']), {}).get('gap', 9999)):
            best[k] = vf

    n_eyes = len(best)
    n_pat  = len(set(v['patient_id'] for v in best.values()))
    log.info('코호트 검증: 278눈/145환자 기준 → 실제 %d눈/%d환자', n_eyes, n_pat)

    # 24-2 필터 제거 건수 리포트 (fail-loud 아님, 보고만)
    n_excl_pattern = sum(1 for vf in vf_rows
                         if vf.get('test_pattern', 'unknown') not in allow_patterns)
    if n_excl_pattern > 0:
        log.warning('test_pattern 필터 제외: %d건 (24-2 아님)', n_excl_pattern)

    if n_eyes != 278:
        raise RuntimeError(
            f'[FAIL] 278 눈 재현 실패: 실제 {n_eyes}눈/{n_pat}환자. '
            f'OCR 누락·test_pattern 오류 여부 확인 필요.'
        )
    log.info('[OK] 278 눈 재현 성공.')


def _try_cfg(*keys):
    """cfg() 호출 시 KeyError 방어."""
    try:
        return cfg(*keys)
    except (KeyError, TypeError):
        return None


def main():
    parser = argparse.ArgumentParser(description='STEP 1 — 원본 이미지 재추출')
    parser.add_argument('--verify-determinism', action='store_true',
                        help='첫 10행을 두 번 실행해 SHA-256 일치 확인 후 종료')
    parser.add_argument('--count-check', action='store_true',
                        help='278 눈 재현 검증 실행 (기본: 실행, --no-count-check로 끔)')
    parser.add_argument('--no-count-check', dest='count_check', action='store_false')
    parser.set_defaults(count_check=True)
    args = parser.parse_args()

    load_config()
    set_seed()
    configure_tesseract()

    out_dir = data_dir()
    cohort  = _load_cohort()
    log.info('cohort 로드: %d행', len(cohort))

    # ── 결정론 검증 모드 ──────────────────────────────────────
    if args.verify_determinism:
        log.info('[determinism] 첫 10행 × 2회 실행 ...')
        _, r1_oct, _, r1_ss = _run_extraction(cohort, limit=10)
        _, r2_oct, _, r2_ss = _run_extraction(cohort, limit=10)
        tmp1 = out_dir / '_det_run1.csv'; tmp2 = out_dir / '_det_run2.csv'
        _write_csv(tmp1, SS_FIELDS, r1_ss)
        _write_csv(tmp2, SS_FIELDS, r2_ss)
        h1 = sha256_file(str(tmp1)); h2 = sha256_file(str(tmp2))
        tmp1.unlink(); tmp2.unlink()
        if h1 == h2:
            log.info('[OK] 결정론 검증 통과: sha256=%s', h1[:16])
        else:
            raise RuntimeError(f'[FAIL] 결정론 실패: run1={h1[:16]} vs run2={h2[:16]}')
        return

    # ── 전체 추출 ─────────────────────────────────────────────
    log.info('전체 %d행 추출 시작 ...', len(cohort))
    vf_rows, oct_rows, rnfl_rows, ss_rows = _run_extraction(cohort)

    # ── CSV 저장 (검증 전에 먼저 저장 — 검증 실패해도 데이터 보존) ──────
    paths = {
        'vf':   out_dir / OUT_VF,
        'oct':  out_dir / OUT_OCT,
        'rnfl': out_dir / OUT_RNFL,
        'ss':   out_dir / OUT_SS,
    }
    _write_csv(paths['vf'],   VF_FIELDS,       vf_rows)
    _write_csv(paths['oct'],  OCT_RAW_FIELDS,  oct_rows)
    _write_csv(paths['rnfl'], RNFL_FIELDS,     rnfl_rows)
    _write_csv(paths['ss'],   SS_FIELDS,        ss_rows)
    log.info('저장 완료: %s', [p.name for p in paths.values()])

    # ── 278 눈 재현 검증 (저장 후 검증) ──────────────────────────
    if args.count_check:
        _cohort_count_check(vf_rows)

    # ── SHA-256 → registry ────────────────────────────────────
    cohort_path = _ROOT / cfg('cohort', 'dataset_csv')
    for key, path in paths.items():
        record_result(
            f'stage10_{key}_sha256',
            value=sha256_file(str(path)),
            script=__file__,
            inputs=[str(cohort_path)],
            extra={'rows': len(vf_rows if key == 'vf' else
                              oct_rows if key == 'oct' else
                              rnfl_rows if key == 'rnfl' else ss_rows),
                   'output_file': path.name},
        )

    log.info('registry 기록 완료. '
             '결정론 검증은: python scripts/stage_10_extract.py --verify-determinism')


if __name__ == '__main__':
    main()
