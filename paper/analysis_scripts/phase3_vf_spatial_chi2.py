"""
scripts/phase3_vf_spatial_chi2.py
VF 추출 오류의 공간 분포 검정 — 수평 자오선 인접 2행 vs 나머지.

원고 §res-accuracy 의 "the two rows flanking the horizontal meridian accounted for
N of M mismatches (x% vs y%, chi2 p<...)" 문장을 만드는 유일한 출처.

2026-08-15 재현성 감사 배경: 이 문장의 수치(155/194, 2.7%)가 원고에만 있고 이를
계산하는 스크립트가 레포에 없었다. 봉인된 results/phase3_extraction_accuracy.json
의 점별 통계로부터 재계산하면 146/185(2.51% vs 0.33%)가 나온다. 출처 없는 숫자를
남기지 않기 위해 계산을 스크립트로 고정한다.

정의:
  - 대상 = 324눈 검증 풀(sfa_thresholds_raw x sfa_review_master_fixed 공통 키),
    phase3_extraction_accuracy.run_vf_accuracy() 의 'normal eyes' 와 동일 집합.
  - mismatch = |raw - gt| != 0 인 유효 쌍(양쪽 다 결측 아님). 결측은 분모/분자
    양쪽에서 제외 — exact_match_pct 의 여집합과 같은 정의.
  - 자오선 인접 2행 = HVF_COORDS 의 row 3, 4 (= p19..p36, 18점).

자기검증: 재계산한 점별 exact_match_pct 가 봉인 JSON 의 per_point_normal 과
54/54 일치하지 않으면 예외를 던진다. 즉 이 스크립트는 봉인 결과와 어긋난 채로는
숫자를 뱉지 못한다.

출력: results/phase3_vf_spatial_chi2.json
"""
import csv, json, datetime, sys
from pathlib import Path

import scipy.stats as sps

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))
sys.path.insert(0, str(_ROOT / 'scripts'))
from hvf.registry import sha256_file                                # noqa: E402
from phase3_extraction_accuracy import _f, accuracy_stats           # noqa: E402
from stage_30_figures import HVF_COORDS                             # noqa: E402

RAW_CSV = _ROOT / 'data' / 'sfa_thresholds_raw.csv'
GT_CSV  = _ROOT / 'sfa_review_master_fixed.csv'
REGISTRY = _ROOT / 'results' / 'phase3_extraction_accuracy.json'
OUT      = _ROOT / 'results' / 'phase3_vf_spatial_chi2.json'

# 수평 자오선을 사이에 두고 맞닿은 두 행. HVF_COORDS 는 8행 격자이고 자오선은
# row 3 과 row 4 사이를 지난다 — 좌표표를 여기서 다시 쓰지 않고 그대로 참조한다.
MERIDIAN_ROWS = (3, 4)


def load_pairs():
    """324눈 검증 풀의 (raw, gt) 쌍을 점별로 모은다."""
    raw_idx = {(r['patient_id'], r['eye'], r['vf_date']): r
               for r in csv.DictReader(open(RAW_CSV, encoding='utf-8-sig'))}
    gt_idx  = {(r['patient_id'], r['eye'], r['vf_date']): r
               for r in csv.DictReader(open(GT_CSV,  encoding='utf-8-sig'))}
    keys = sorted(set(raw_idx) & set(gt_idx))

    # 오프셋 눈 판정 — phase3_extraction_accuracy.run_vf_accuracy() 와 동일 규칙.
    # 2026-07-23 row3/4 패치 이후 실측 0건이지만, 규칙을 빼면 과거 데이터에서
    # 조용히 달라지므로 그대로 둔다.
    def count_match(rr, gr, pts, shift=0):
        n = 0
        for p in pts:
            rp = p + shift
            if rp < 1 or rp > 54: continue
            rv = _f(rr.get(f'p{rp:02d}') or rr.get(f'p{rp}'))
            gv = _f(gr.get(f'p{p:02d}')  or gr.get(f'p{p}'))
            if rv is not None and gv is not None and abs(rv - gv) <= 1:
                n += 1
        return n

    post19 = list(range(19, 55))
    normal, offset = [], []
    for k in keys:
        rr, gr = raw_idx[k], gt_idx[k]
        if count_match(rr, gr, post19, -1) > count_match(rr, gr, post19, 0) + 5:
            offset.append(k)
        else:
            normal.append(k)

    by_point = {}
    for p in range(1, 55):
        by_point[p] = [(_f(raw_idx[k].get(f'p{p:02d}') or raw_idx[k].get(f'p{p}')),
                        _f(gt_idx[k].get(f'p{p:02d}')  or gt_idx[k].get(f'p{p}')))
                       for k in normal]
    return by_point, len(keys), len(normal), len(offset)


def verify_against_registry(by_point):
    """점별 통계가 봉인 JSON 과 일치하는지 확인. 어긋나면 즉시 실패."""
    if not REGISTRY.exists():
        raise FileNotFoundError(f'봉인 결과가 없어 교차검증 불가: {REGISTRY}')
    sealed = json.load(open(REGISTRY, encoding='utf-8'))['vf_accuracy']['per_point_normal']
    bad = [p for p in range(1, 55)
           if accuracy_stats(by_point[p]) != sealed[f'p{p:02d}']]
    if bad:
        raise ValueError(
            f'점별 통계가 봉인 JSON 과 불일치: {bad[:10]} (총 {len(bad)}점). '
            'phase3_extraction_accuracy.py 를 먼저 재실행해 정본을 맞출 것.')
    return len(sealed)


def tally(by_point, points):
    valid = mismatch = 0
    for p in points:
        for rv, gv in by_point[p]:
            if rv is None or gv is None:
                continue
            valid += 1
            if abs(rv - gv) != 0:
                mismatch += 1
    return valid, mismatch


def main():
    by_point, n_common, n_normal, n_offset = load_pairs()
    n_checked = verify_against_registry(by_point)
    print(f'검증 풀: common={n_common}눈  normal={n_normal}  offset={n_offset}')
    print(f'봉인 JSON 점별 교차검증: {n_checked}/54 일치')

    mid = sorted(i + 1 for i, (r, _) in enumerate(HVF_COORDS) if r in MERIDIAN_ROWS)
    rest = sorted(i + 1 for i, (r, _) in enumerate(HVF_COORDS) if r not in MERIDIAN_ROWS)
    v_mid, m_mid = tally(by_point, mid)
    v_rest, m_rest = tally(by_point, rest)
    r_mid, r_rest = 100 * m_mid / v_mid, 100 * m_rest / v_rest

    print(f'\n자오선 인접 2행 = {len(mid)}점 (p{mid[0]:02d}..p{mid[-1]:02d})')
    print(f'  {"그룹":22s} {"유효점":>8s} {"불일치":>7s} {"비율%":>7s}')
    print('  ' + '-' * 48)
    print(f'  {"자오선 인접 2행":22s} {v_mid:8d} {m_mid:7d} {r_mid:7.2f}')
    print(f'  {"나머지":22s} {v_rest:8d} {m_rest:7d} {r_rest:7.2f}')
    print(f'  {"전체":22s} {v_mid+v_rest:8d} {m_mid+m_rest:7d} '
          f'{100*(m_mid+m_rest)/(v_mid+v_rest):7.2f}')

    table = [[m_mid, v_mid - m_mid], [m_rest, v_rest - m_rest]]
    chi2_y, p_y, dof, exp = sps.chi2_contingency(table, correction=True)
    chi2_n, p_n, _, _     = sps.chi2_contingency(table, correction=False)
    odds, p_f             = sps.fisher_exact(table)

    print(f'\n  chi2 (Yates 보정)  = {chi2_y:.4f}  df={dof}  P={p_y:.3e}   <- 원고 보고값')
    print(f'  chi2 (보정 없음)   = {chi2_n:.4f}  df={dof}  P={p_n:.3e}')
    print(f'  Fisher exact       = OR {odds:.4f}          P={p_f:.3e}')
    print(f'  최소 기대도수      = {exp.min():.2f}')

    out = {
        'generated': datetime.date.today().isoformat(),
        'input_sha256': {
            'sfa_thresholds_raw': sha256_file(str(RAW_CSV)),
            'sfa_review_master_fixed': sha256_file(str(GT_CSV)),
        },
        'design': {
            'pool': f'{n_normal}-eye validation pool (offset eyes: {n_offset})',
            'mismatch_definition': 'abs(raw - gt) != 0 among pairs valid in both',
            'meridian_rows': f'HVF_COORDS rows {MERIDIAN_ROWS} = p{mid[0]:02d}..p{mid[-1]:02d} ({len(mid)} points)',
            'cross_check': f'per-point stats verified identical to sealed registry ({n_checked}/54)',
        },
        'counts': {
            'meridian_adjacent': {'valid': v_mid,  'mismatch': m_mid,  'rate_pct': round(r_mid, 2)},
            'elsewhere':         {'valid': v_rest, 'mismatch': m_rest, 'rate_pct': round(r_rest, 2)},
            'total':             {'valid': v_mid + v_rest, 'mismatch': m_mid + m_rest},
        },
        'chi2_yates':    {'statistic': round(float(chi2_y), 4), 'dof': int(dof), 'p': float(p_y)},
        'chi2_uncorrected': {'statistic': round(float(chi2_n), 4), 'dof': int(dof), 'p': float(p_n)},
        'fisher_exact':  {'odds_ratio': round(float(odds), 4), 'p': float(p_f)},
        'min_expected_cell': round(float(exp.min()), 2),
    }
    json.dump(out, open(OUT, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print(f'\n결과 저장: {OUT}')


if __name__ == '__main__':
    main()
