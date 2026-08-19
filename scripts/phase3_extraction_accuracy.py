"""
scripts/phase3_extraction_accuracy.py
Phase 3: OCR 추출 정확도 정량 검증

소스:
  VF raw     : data/sfa_thresholds_raw.csv     (340행)
  VF GT      : sfa_review_master_fixed.csv     (340행, 수동검수)
  OCT raw    : data/oct_values_raw.csv         (340행)
  RNFL raw   : data/rnfl_detail_raw.csv        (340행)
  OCT/RNFL GT: oct_tabular_90d.csv             (280행, 수동검수)
  Taxonomy   : results/ocr_failure_taxonomy.json

오프셋 처리:
  54눈에서 raw p19+ 1칸 밀림 → 자동 감지(shift_improves_match) 후 재매핑

출력: results/phase3_extraction_accuracy.json
"""
import csv, json, math, datetime, sys
from pathlib import Path
from collections import defaultdict

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))
from hvf.registry import sha256_file

def _f(v, sentinel_none=('-1',)):
    """VF 값 변환. -1 → 0 dB (절대암점). None/blank → None."""
    if v in ('', None, 'nan', 'None', 'N/A'): return None
    try:
        x = float(v)
        if math.isnan(x): return None
        if x == -1: return 0.0   # 절대암점
        return x
    except (TypeError, ValueError):
        return None

def _fnum(v):
    """숫자 변환만, -1 포함 그대로."""
    if v in ('', None, 'nan', 'None', 'N/A'): return None
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except: return None

# ─── 정확도 계산 유틸 ─────────────────────────────────────────────
def accuracy_stats(pairs):
    """
    pairs = [(raw_val, gt_val), ...] — None 포함 가능
    returns dict with mae, rmse, exact, within1, within2, n_valid, n_raw_missing, n_gt_missing
    """
    n_raw_miss = sum(1 for r, g in pairs if r is None and g is not None)
    n_gt_miss  = sum(1 for r, g in pairs if g is None)
    valid = [(r, g) for r, g in pairs if r is not None and g is not None]
    n = len(valid)
    if n == 0:
        return {'n_pairs': len(pairs), 'n_valid': 0,
                'n_raw_missing': n_raw_miss, 'n_gt_missing': n_gt_miss}
    diffs  = [abs(r - g) for r, g in valid]
    sq     = [d**2 for d in diffs]
    exact  = sum(1 for d in diffs if d == 0)
    w1     = sum(1 for d in diffs if d <= 1)
    w2     = sum(1 for d in diffs if d <= 2)
    mae    = sum(diffs) / n
    rmse   = math.sqrt(sum(sq) / n)
    return {
        'n_pairs': len(pairs), 'n_valid': n,
        'n_raw_missing': n_raw_miss, 'n_gt_missing': n_gt_miss,
        'missing_rate_raw_pct': round(100 * n_raw_miss / max(len(pairs), 1), 1),
        'exact_match_pct':      round(100 * exact / n, 1),
        'within1dB_pct':        round(100 * w1   / n, 1),
        'within2dB_pct':        round(100 * w2   / n, 1),
        'mae':                  round(mae, 3),
        'rmse':                 round(rmse, 3),
    }


# ════════════════════════════════════════════════════════════════
# SECTION A: VF 추출 정확도
# ════════════════════════════════════════════════════════════════
def run_vf_accuracy():
    raw_rows = list(csv.DictReader(
        open(_ROOT / 'data' / 'sfa_thresholds_raw.csv', encoding='utf-8-sig')))
    gt_rows  = list(csv.DictReader(
        open(_ROOT / 'sfa_review_master_fixed.csv',     encoding='utf-8-sig')))

    # 인덱스: (patient_id, eye, vf_date) → row
    raw_idx = {(r['patient_id'], r['eye'], r['vf_date']): r for r in raw_rows}
    gt_idx  = {(r['patient_id'], r['eye'], r['vf_date']): r for r in gt_rows}

    common_keys = sorted(set(raw_idx) & set(gt_idx))
    only_raw    = set(raw_idx) - set(gt_idx)
    only_gt     = set(gt_idx)  - set(raw_idx)
    print(f'VF: raw={len(raw_rows)} gt={len(gt_rows)} common={len(common_keys)} '
          f'only_raw={len(only_raw)} only_gt={len(only_gt)}')

    P = list(range(1, 55))   # p1..p54

    # ── 오프셋 감지: 각 눈에서 p19-p54 1칸 shift가 일치율을 높이는지 ──
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

    offset_eyes = []
    normal_eyes = []
    for key in common_keys:
        rr, gr = raw_idx[key], gt_idx[key]
        pts_post19 = list(range(19, 55))
        m0 = count_match(rr, gr, pts_post19, shift=0)
        m1 = count_match(rr, gr, pts_post19, shift=-1)  # raw에서 1칸 당겨야
        if m1 > m0 + 5:   # shift가 명확히 더 나을 때 (5점 이상 개선)
            offset_eyes.append(key)
        else:
            normal_eyes.append(key)
    print(f'VF offset eyes: {len(offset_eyes)}/324  normal: {len(normal_eyes)}')

    def get_raw_val(rr, p, is_offset):
        rp = p - 1 if (is_offset and p >= 19) else p
        if rp < 1 or rp > 54: return None
        return _f(rr.get(f'p{rp:02d}') or rr.get(f'p{rp}'))

    def get_gt_val(gr, p):
        return _f(gr.get(f'p{p:02d}') or gr.get(f'p{p}'))

    # ── 전체 54점 정확도 (정상+재매핑 눈) ──
    all_pairs = []
    for key in common_keys:
        rr, gr = raw_idx[key], gt_idx[key]
        is_off = key in set(offset_eyes)
        for p in P:
            rv = get_raw_val(rr, p, is_off)
            gv = get_gt_val(gr, p)
            all_pairs.append((rv, gv))

    # ── 정상 눈만 ──
    norm_pairs = []
    for key in normal_eyes:
        rr, gr = raw_idx[key], gt_idx[key]
        for p in P:
            rv = _f(rr.get(f'p{p:02d}') or rr.get(f'p{p}'))
            gv = get_gt_val(gr, p)
            norm_pairs.append((rv, gv))

    # ── 오프셋 눈(재매핑 후) ──
    off_pairs = []
    for key in offset_eyes:
        rr, gr = raw_idx[key], gt_idx[key]
        for p in P:
            rv = get_raw_val(rr, p, True)
            gv = get_gt_val(gr, p)
            off_pairs.append((rv, gv))

    # ── 포인트별 정확도 (전체, normal만) ──
    pt_stats = {}
    for p in P:
        pairs_p = []
        for key in normal_eyes:
            rr, gr = raw_idx[key], gt_idx[key]
            rv = _f(rr.get(f'p{p:02d}') or rr.get(f'p{p}'))
            gv = get_gt_val(gr, p)
            pairs_p.append((rv, gv))
        pt_stats[f'p{p:02d}'] = accuracy_stats(pairs_p)

    # 누락률 높은 포인트
    high_miss = sorted(
        [(p, pt_stats[f'p{p:02d}']['missing_rate_raw_pct'])
         for p in P if pt_stats[f'p{p:02d}']['missing_rate_raw_pct'] > 0],
        key=lambda x: -x[1]
    )[:10]

    vf_result = {
        'matching': {
            'n_raw': len(raw_rows), 'n_gt': len(gt_rows),
            'n_common': len(common_keys),
            'n_offset_eyes': len(offset_eyes),
            'n_normal_eyes': len(normal_eyes),
        },
        'accuracy_all_eyes_post_remap': accuracy_stats(all_pairs),
        'accuracy_normal_eyes': accuracy_stats(norm_pairs),
        'accuracy_offset_eyes_after_remap': accuracy_stats(off_pairs),
        'per_point_normal': pt_stats,
        'high_missing_points': high_miss,
    }

    # 콘솔 요약
    print('\n── VF 정확도 요약 ──────────────────────────────────')
    for label, stat in [
        ('전체(재매핑 포함)', vf_result['accuracy_all_eyes_post_remap']),
        ('정상 눈만',         vf_result['accuracy_normal_eyes']),
        ('오프셋 눈(재매핑)', vf_result['accuracy_offset_eyes_after_remap']),
    ]:
        s = stat
        if s.get('n_valid', 0) == 0:
            print(f'  {label:20s}: n=0 (해당 없음 — 2026-07-23 row3/4 패치 후 오프셋 눈 0건)')
            continue
        print(f'  {label:20s}: exact={s.get("exact_match_pct","?"):5.1f}%  '
              f'±1dB={s.get("within1dB_pct","?"):5.1f}%  '
              f'MAE={s.get("mae","?"):.3f}  RMSE={s.get("rmse","?"):.3f}  '
              f'miss%={s.get("missing_rate_raw_pct","?")}  n={s.get("n_valid","?")}')
    if high_miss:
        print(f'  누락률 상위: {high_miss[:5]}')
    return vf_result


# ════════════════════════════════════════════════════════════════
# SECTION B: GCL 추출 정확도
# ════════════════════════════════════════════════════════════════
def run_gcl_accuracy():
    raw_rows = list(csv.DictReader(
        open(_ROOT / 'data' / 'oct_values_raw.csv', encoding='utf-8-sig')))
    gt_rows  = list(csv.DictReader(
        open(_ROOT / 'oct_tabular_90d.csv',         encoding='utf-8-sig')))

    raw_idx = {(r['patient_id'], r['eye'], r['oct_date']): r for r in raw_rows}
    gt_idx  = {(r['patient_id'], r['eye'], r['oct_date']): r for r in gt_rows}

    common = sorted(set(raw_idx) & set(gt_idx))
    print(f'\nGCL: raw={len(raw_rows)} gt={len(gt_rows)} common={len(common)}')

    GCL_COLS = [
        'avg_gcl_od', 'avg_gcl_os', 'min_gcl_od', 'min_gcl_os',
        'od_s_sup', 'od_s_sup_t', 'od_s_inf_t', 'od_s_inf', 'od_s_inf_n', 'od_s_sup_n',
        'os_s_sup', 'os_s_sup_t', 'os_s_inf_t', 'os_s_inf', 'os_s_inf_n', 'os_s_sup_n',
        'ss_gca_od', 'ss_gca_os',
    ]

    # 2026-07-17: GT(oct_tabular_90d.csv)의 os_s_{sup_t,sup_n,inf_t,inf_n} 4개
    # 컬럼은 행의 eye에 따라 다른 컨벤션으로 기록됨(coordinate_frame_convention.md
    # §2-2). 어느 쪽(OD행 vs OS행)이 raw를 그대로 쓰고 어느 쪽이 스왑해야
    # 하는지는 raw 추출 코드의 현재 컨벤션(eye-dependent mirror 적용 여부)에
    # 따라 바뀐다 — 2026-07-16 하루에도 코드가 두 번 뒤집히며 방향이 반대로
    # 바뀐 전례가 있어 하드코딩하지 않고 두 방향을 모두 계산해 GT와 더 잘 맞는
    # 쪽을 자동 채택한다(scripts/verify_gca_os_mirroring.py와 동일 원리).
    SWAP_PAIRS = [('os_s_sup_t', 'os_s_sup_n'), ('os_s_inf_t', 'os_s_inf_n')]
    SWAP_MAP = {}
    for a, b in SWAP_PAIRS:
        SWAP_MAP[a] = b
        SWAP_MAP[b] = a
    SWAP_COLS = set(SWAP_MAP)

    def _score_swap_cols(od_swap, os_swap):
        n = ex = 0
        for col in SWAP_COLS:
            for key in common:
                eye = key[1]
                swap_here = od_swap if eye == 'OD' else os_swap
                src_col = SWAP_MAP[col] if swap_here else col
                rv = _fnum(raw_idx[key].get(src_col))
                gv = _fnum(gt_idx[key].get(col))
                if rv is None or gv is None:
                    continue
                n += 1
                if rv == gv:
                    ex += 1
        return ex / max(n, 1)

    pct_a = _score_swap_cols(od_swap=False, os_swap=True)
    pct_b = _score_swap_cols(od_swap=True, os_swap=False)
    OD_SWAP, OS_SWAP = (False, True) if pct_a >= pct_b else (True, False)
    swap_convention = (
        f'auto-detected: OD_swap={OD_SWAP}, OS_swap={OS_SWAP} '
        f'(A: OD-asis/OS-swap={pct_a*100:.1f}% vs B: OD-swap/OS-asis={pct_b*100:.1f}%)'
    )
    print(f'  [GCA swap 방향 자동판별] {swap_convention}')

    def _raw_col_for(col, eye):
        """row-eye 조건부 컬럼 선택 (자동판별된 방향 적용)."""
        if col not in SWAP_MAP:
            return col
        swap_here = OD_SWAP if eye == 'OD' else OS_SWAP
        return SWAP_MAP[col] if swap_here else col

    col_stats = {}
    print(f'\n── GCL/SS 컬럼별 정확도 (row-eye 조건부 채점) ──────────')
    print(f'  {"컬럼":20s}  exact%  ±1%   MAE    miss%  n')
    print('  '+'-'*65)
    def _fmt(st):
        ex  = f'{st["exact_match_pct"]:5.1f}' if isinstance(st.get("exact_match_pct"), float) else '  N/A'
        w1  = f'{st["within1dB_pct"]:5.1f}'   if isinstance(st.get("within1dB_pct"), float)   else '  N/A'
        mae = f'{st["mae"]:6.3f}'              if isinstance(st.get("mae"), float)             else '   N/A'
        miss= f'{st.get("missing_rate_raw_pct","?")}'
        return ex, w1, mae, miss

    for col in GCL_COLS:
        pairs = []
        for key in common:
            eye = key[1]
            src_col = _raw_col_for(col, eye)
            rv = _fnum(raw_idx[key].get(src_col))
            gv = _fnum(gt_idx[key].get(col))
            pairs.append((rv, gv))
        st = accuracy_stats(pairs)
        col_stats[col] = st
        ex, w1, mae, miss = _fmt(st)
        print(f'  {col:20s}  {ex}%  {w1}%  {mae}  {miss}%  {st.get("n_valid","?")}')

    gcl_result = {
        'matching': {'n_raw': len(raw_rows), 'n_gt': len(gt_rows), 'n_common': len(common)},
        'scoring_note': (
            'os_s_{sup_t,sup_n,inf_t,inf_n}: row-eye 조건부 스왑 적용, 방향은 '
            'GT 일치율 기준 자동판별(raw 추출 코드의 미러링 상태에 따라 달라짐). '
            'coordinate_frame_convention.md 참고.'
        ),
        'swap_convention_detected': swap_convention,
        'per_column': col_stats,
    }
    return gcl_result


# ════════════════════════════════════════════════════════════════
# SECTION C: RNFL 추출 정확도
# ════════════════════════════════════════════════════════════════
def run_rnfl_accuracy():
    raw_rows = list(csv.DictReader(
        open(_ROOT / 'data' / 'rnfl_detail_raw.csv', encoding='utf-8-sig')))
    gt_rows  = list(csv.DictReader(
        open(_ROOT / 'oct_tabular_90d.csv',          encoding='utf-8-sig')))

    # oct_tabular_90d에는 vf_date도 있지만, RNFL은 oct_date로만 매칭
    raw_idx = {(r['patient_id'], r['eye'], r['oct_date']): r for r in raw_rows}
    gt_idx  = {(r['patient_id'], r['eye'], r['oct_date']): r for r in gt_rows}

    common = sorted(set(raw_idx) & set(gt_idx))
    print(f'\nRNFL: raw={len(raw_rows)} gt={len(gt_rows)} common={len(common)}')

    RNFL_COLS = (
        ['rnfl_q_s', 'rnfl_q_t', 'rnfl_q_i', 'rnfl_q_n'] +
        [f'rnfl_h{h:02d}' for h in range(1, 13)]
    )

    # 2026-07-17: Method B — GT의 rnfl_q_t/rnfl_q_n은 OS행에서 OD-normalized
    # 컨벤션으로 기록됨(coordinate_frame_convention.md; scripts/
    # verify_rnfl_method_b.py로 최초 확정). OS행에 한해 두 컬럼을 스왑해야
    # on-screen raw와 올바르게 대조됨(78.4% 재현, 미적용 시 왜곡됨).
    RNFL_SWAP = {'rnfl_q_t': 'rnfl_q_n', 'rnfl_q_n': 'rnfl_q_t'}

    # 2026-08-16: clock-hour 에도 같은 원칙을 적용한다. 지금까지 시각 12개는
    # 아무 변환 없이 대조돼 OS 139안이 통째로 오답 처리됐다(OS 11.5% vs OD 89.8%).
    # 근거: 원본 리포트 육안 대조에서 raw 가 화면 순서와 일치함을 확인했고
    # (OS 2건, h01~h11 11/12 일치), GT 는 그 역순이었다 → raw=on-screen,
    # GT=OD-normalized. 사분면 T/N 과 정확히 같은 구조.
    # 매핑은 coordinate_frame_convention.md §3.1 (수직축 고정: H12·H06 그대로,
    # H0n ↔ H(12−n)). make_rnfl_flip.py 의 H01/H07 고정 매핑은 쓰지 않는다 —
    # 문서와 한 칸 어긋나고 실측 일치율이 1.0% 였다.
    RNFL_SWAP.update({f'rnfl_h{n:02d}': f'rnfl_h{12 - n:02d}' for n in range(1, 12)})

    def _rnfl_col_for(col, eye):
        if eye == 'OS' and col in RNFL_SWAP:
            return RNFL_SWAP[col]
        return col

    col_stats = {}
    print(f'\n── RNFL 컬럼별 정확도 (Method B: OS행 T/N 스왑) ─────────')
    print(f'  {"컬럼":12s}  exact%  ±1%   MAE    miss%  n')
    print('  '+'-'*58)
    def _fmt2(st):
        ex  = f'{st["exact_match_pct"]:5.1f}' if isinstance(st.get("exact_match_pct"), float) else '  N/A'
        w1  = f'{st["within1dB_pct"]:5.1f}'   if isinstance(st.get("within1dB_pct"), float)   else '  N/A'
        mae = f'{st["mae"]:6.3f}'              if isinstance(st.get("mae"), float)             else '   N/A'
        miss= f'{st.get("missing_rate_raw_pct","?")}'
        return ex, w1, mae, miss

    for col in RNFL_COLS:
        pairs = []
        for key in common:
            eye = key[1]
            src_col = _rnfl_col_for(col, eye)
            rv = _fnum(raw_idx[key].get(src_col))
            gv = _fnum(gt_idx[key].get(col))
            pairs.append((rv, gv))
        st = accuracy_stats(pairs)
        col_stats[col] = st
        ex, w1, mae, miss = _fmt2(st)
        print(f'  {col:12s}  {ex}%  {w1}%  {mae}  {miss}%  {st.get("n_valid","?")}')

    rnfl_result = {
        'matching': {'n_raw': len(raw_rows), 'n_gt': len(gt_rows), 'n_common': len(common)},
        'scoring_note': (
            'OS행에 한해 raw(on-screen)를 GT(OD-normalized) 프레임으로 변환 후 대조. '
            'rnfl_q_t/rnfl_q_n: Method B 스왑. '
            'rnfl_h01~h12: 수직축 고정 미러(H12·H06 고정, H0n<->H(12-n)), 2026-08-16 추가. '
            'coordinate_frame_convention.md §3.1 참고.'
        ),
        'per_column': col_stats,
    }
    return rnfl_result


# ════════════════════════════════════════════════════════════════
# SECTION D: 오류 taxonomy 요약
# ════════════════════════════════════════════════════════════════
def summarize_taxonomy():
    tax = json.load(open(_ROOT / 'results' / 'ocr_failure_taxonomy.json', encoding='utf-8'))
    types = {}
    for k, v in tax.get('taxonomy', {}).items():
        types[k] = {
            'label':     v.get('label', ''),
            'fixable':   v.get('fixable'),
            'count_est': v.get('estimated_count_in_158', ''),
            'description_short': v.get('description', '')[:120],
        }
    return {
        'source': 'results/ocr_failure_taxonomy.json',
        'taxonomy_types': types,
        'summary': tax.get('summary', {}),
    }


# ════════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════════
def main():
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    print('='*65)
    print('Phase 3: OCR 추출 정확도 검증')
    print('='*65)

    results = {
        'generated': datetime.date.today().isoformat(),
        'source_sha256': {
            'sfa_raw':  sha256_file(str(_ROOT / 'data' / 'sfa_thresholds_raw.csv')),
            'sfa_gt':   sha256_file(str(_ROOT / 'sfa_review_master_fixed.csv')),
            'oct_raw':  sha256_file(str(_ROOT / 'data' / 'oct_values_raw.csv')),
            'rnfl_raw': sha256_file(str(_ROOT / 'data' / 'rnfl_detail_raw.csv')),
            'oct_gt':   sha256_file(str(_ROOT / 'oct_tabular_90d.csv')),
        },
    }

    results['vf_accuracy']    = run_vf_accuracy()
    results['gcl_accuracy']   = run_gcl_accuracy()
    results['rnfl_accuracy']  = run_rnfl_accuracy()
    results['error_taxonomy'] = summarize_taxonomy()

    # ── 최종 요약 테이블 ──
    print('\n' + '='*65)
    print('Phase 3 결과 요약 테이블')
    print('='*65)

    print('\n▶ VF OCR 정확도 (sfa_thresholds_raw vs master_fixed)')
    va = results['vf_accuracy']
    for label, key in [('전체(재매핑)', 'accuracy_all_eyes_post_remap'),
                       ('정상 눈',      'accuracy_normal_eyes'),
                       ('오프셋 눈',    'accuracy_offset_eyes_after_remap')]:
        s = va[key]
        if s.get('n_valid', 0) == 0:
            print(f'  {label:16s}  n=0 (해당 없음)')
            continue
        print(f'  {label:16s}  exact={s.get("exact_match_pct","?"):5.1f}%  '
              f'±1dB={s.get("within1dB_pct","?"):5.1f}%  '
              f'MAE={s.get("mae","?"):.2f}dB  '
              f'miss={s.get("missing_rate_raw_pct","?")}%  '
              f'n_pts={s.get("n_valid","?")}')
    print(f'  오프셋 눈 수: {va["matching"]["n_offset_eyes"]}/'
          f'{va["matching"]["n_common"]} ({100*va["matching"]["n_offset_eyes"]/max(va["matching"]["n_common"],1):.0f}%)')

    def _s(st, key, fmt='.1f'):
        v = st.get(key)
        if isinstance(v, (int, float)):
            return format(v, fmt)
        return 'N/A'

    print('\n▶ GCL OCR 정확도 (oct_values_raw vs oct_tabular_90d)')
    ga = results['gcl_accuracy']['per_column']
    key_cols = ['avg_gcl_od', 'avg_gcl_os', 'min_gcl_od', 'min_gcl_os',
                'od_s_sup_t', 'od_s_inf_t', 'od_s_inf',
                'os_s_sup_t', 'os_s_inf_t', 'os_s_inf', 'ss_gca_od', 'ss_gca_os']
    for col in key_cols:
        if col in ga:
            s = ga[col]
            print(f'  {col:20s}  exact={_s(s,"exact_match_pct"):5s}%  '
                  f'MAE={_s(s,"mae",".3f"):6s}  '
                  f'miss={_s(s,"missing_rate_raw_pct",".1f")}%  n={s.get("n_valid","?")}')

    print('\n▶ RNFL OCR 정확도 (rnfl_detail_raw vs oct_tabular_90d)')
    ra = results['rnfl_accuracy']['per_column']
    for col in ['rnfl_q_s', 'rnfl_q_t', 'rnfl_q_i', 'rnfl_q_n'] + \
               [f'rnfl_h{h:02d}' for h in range(1, 13)]:
        if col in ra:
            s = ra[col]
            print(f'  {col:12s}  exact={_s(s,"exact_match_pct"):5s}%  '
                  f'MAE={_s(s,"mae",".3f"):7s}  '
                  f'miss={_s(s,"missing_rate_raw_pct",".1f")}%  n={s.get("n_valid","?")}')

    # JSON 먼저 저장
    out = _ROOT / 'results' / 'phase3_extraction_accuracy.json'
    json.dump(results, open(out, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print(f'\n출력: {out}')

    print('\n▶ OCR 오류 Taxonomy')
    for k, v in results['error_taxonomy']['taxonomy_types'].items():
        fix = 'fixable' if v['fixable'] else 'NOT fixable'
        label = v['label'].encode('utf-8', errors='replace').decode('utf-8')
        print(f'  {k}: {v["count_est"]}  [{fix}]')
        print(f'    {label[:80]}')


if __name__ == '__main__':
    main()
