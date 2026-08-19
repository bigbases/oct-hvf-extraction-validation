"""
scripts/phase5_2_sector_vf.py
Phase 5-2: 섹터별 구조-기능 공간 매핑 상관 (Garway-Heath 기반)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
좌표 규약 (OD-unified, 검증 완료 2026-07-05):
  VF 좌측(-x) = nasal for OD  │  VF 우측(+x) = temporal for OD
  VF 상단(+y) = superior      │  VF 하단(-y) = inferior
  OS 처리: vf[OS] = flip_vf(vf[OS])  /  oct[OS] 섹터 = 그대로
  근거: analysis_master OS VF p19(x=-27) 평균 22.0 dB > p27(x=+27) 14.4 dB
        → OS 원본은 temporal=LEFT(-x)로 저장 → flip 후 temporal=RIGHT(+x)로 통일

망막→VF 매핑 원칙 (Garway-Heath 2000):
  디스크/망막 각도 α  →  VF 각도 (α + 180°) mod 360°
  (180° 회전: 상하반전 + 좌우반전)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
VF 54점 좌표 (OD/unified 기준, blind spot x=+15 제외):
  Row0 y=+21: p01(-9,21) p02(-3,21) p03(+3,21) p04(+9,21)
  Row1 y=+15: p05(-15,15) p06(-9,15) p07(-3,15) p08(+3,15) p09(+9,15) p10(+15,15)
  Row2 y=+09: p11(-21,9) p12(-15,9) p13(-9,9) p14(-3,9) p15(+3,9) p16(+9,9) p17(+15,9) p18(+21,9)
  Row3 y=+03: p19(-27,3) p20(-21,3) p21(-15,3) p22(-9,3) p23(-3,3) p24(+3,3) p25(+9,3) p26(+21,3) p27(+27,3)
  Row4 y=-03: p28(-27,-3) p29(-21,-3) p30(-15,-3) p31(-9,-3) p32(-3,-3) p33(+3,-3) p34(+9,-3) p35(+21,-3) p36(+27,-3)
  Row5 y=-09: p37(-21,-9) p38(-15,-9) p39(-9,-9) p40(-3,-9) p41(+3,-9) p42(+9,-9) p43(+15,-9) p44(+21,-9)
  Row6 y=-15: p45(-15,-15) p46(-9,-15) p47(-3,-15) p48(+3,-15) p49(+9,-15) p50(+15,-15)
  Row7 y=-21: p51(-9,-21) p52(-3,-21) p53(+3,-21) p54(+9,-21)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GCA 6섹터 → VF 매핑 (60° 슬라이스, 180° 회전):
  섹터명    망막 각도       VF 각도         VF 영역
  ─────────────────────────────────────────────────────
  sup_t   [  0°, 60°)  → [180°,240°)  = inferonasal  (하비측)
  sup     [ 60°,120°)  → [240°,300°)  = inferior     (하측 중앙)
  sup_n   [120°,180°)  → [300°,360°)  = inferotemporal (하측두)
  inf_n   [180°,240°)  → [  0°, 60°)  = superotemporal (상측두)
  inf     [240°,300°)  → [ 60°,120°)  = superior     (상측 중앙)
  inf_t   [300°,360°)  → [120°,180°)  = superonasal  (상비측)

RNFL 4사분면 → VF 매핑 (hemifield, 검증됨):
  q_s  상방 disc   → inferior VF  (p28-p54, 하측 전체)
  q_i  하방 disc   → superior VF  (p01-p27, 상측 전체)
  q_t  측두 disc   → central VF   (황반 ±10°, 16pts)
  q_n  비측 disc   → peripheral temporal VF (x≥+15°, 10pts)

RNFL clock-hour → VF 매핑 (30° 슬라이스, 180° 회전):
  disc h12=90° → VF 270° (하중앙)
  disc h09=180° → VF 0°  (측두 peripheral)
  disc h06=270° → VF 90° (상중앙)
  disc h03=0°  → VF 180° (비측 중앙, PMB)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
출력: results/phase5_2_sector_correlations.json
"""
import csv, json, math, datetime, sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import scipy.stats as sps

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))
from hvf.registry import sha256_file
from hvf.cohort import filter_rows

# ═══════════════════════════════════════════════════════════════
# 1.  VF 54점 좌표 정의
# ═══════════════════════════════════════════════════════════════

# (x°, y°) in OD-unified space  LEFT=-x=nasal, RIGHT=+x=temporal
VF_COORDS = {
    1:(-9,21),  2:(-3,21),  3:(3,21),   4:(9,21),
    5:(-15,15), 6:(-9,15),  7:(-3,15),  8:(3,15),   9:(9,15),  10:(15,15),
   11:(-21,9), 12:(-15,9), 13:(-9,9),  14:(-3,9),  15:(3,9),  16:(9,9),  17:(15,9), 18:(21,9),
   19:(-27,3), 20:(-21,3), 21:(-15,3), 22:(-9,3),  23:(-3,3), 24:(3,3),  25:(9,3),  26:(21,3), 27:(27,3),
   28:(-27,-3),29:(-21,-3),30:(-15,-3),31:(-9,-3), 32:(-3,-3),33:(3,-3), 34:(9,-3), 35:(21,-3),36:(27,-3),
   37:(-21,-9),38:(-15,-9),39:(-9,-9), 40:(-3,-9), 41:(3,-9), 42:(9,-9), 43:(15,-9),44:(21,-9),
   45:(-15,-15),46:(-9,-15),47:(-3,-15),48:(3,-15),49:(9,-15),50:(15,-15),
   51:(-9,-21),52:(-3,-21),53:(3,-21), 54:(9,-21),
}
assert len(VF_COORDS) == 54

# VF flip rows (make_flip.py VF_ROWS와 동일)
_VF_FLIP_ROWS = [
    list(range(1,5)), list(range(5,11)), list(range(11,19)),
    list(range(19,28)), list(range(28,37)), list(range(37,45)),
    list(range(45,51)), list(range(51,55)),
]

def _vf_angle(p):
    """VF point의 각도 (0°=temporal=+x, 90°=superior=+y, [0,360))"""
    x, y = VF_COORDS[p]
    a = math.degrees(math.atan2(y, x))
    return a % 360.0

def _pts_in_vf_range(lo, hi):
    """VF 각도 [lo, hi)에 해당하는 VF 점 목록"""
    result = []
    for p in range(1, 55):
        a = _vf_angle(p)
        if lo < hi:
            if lo <= a < hi:
                result.append(p)
        else:  # 0° 걸침 (hi < lo, e.g., h09: 345°-15°)
            if a >= lo or a < hi:
                result.append(p)
    return sorted(result)


# ═══════════════════════════════════════════════════════════════
# 2.  GCA 6섹터 매핑 정의
#     망막 각도 슬라이스 → VF 슬라이스 = 망막 + 180°
#     형식: {column_suffix: (vf_lo, vf_hi, description)}
# ═══════════════════════════════════════════════════════════════

GCA_VF_DEF = {
    # 망막 [0,60)   → VF [180,240)  = inferonasal
    'sup_t': (180, 240, 'inferonasal VF  (하비측)'),
    # 망막 [60,120)  → VF [240,300)  = inferior
    'sup':   (240, 300, 'inferior VF     (하측 중앙)'),
    # 망막 [120,180) → VF [300,360)  = inferotemporal
    'sup_n': (300, 360, 'inferotemporal VF (하측두)'),
    # 망막 [180,240) → VF [0,60)     = superotemporal
    'inf_n': (  0,  60, 'superotemporal VF (상측두)'),
    # 망막 [240,300) → VF [60,120)   = superior
    'inf':   ( 60, 120, 'superior VF     (상측 중앙)'),
    # 망막 [300,360) → VF [120,180)  = superonasal
    'inf_t': (120, 180, 'superonasal VF  (상비측)'),
}

# VF 점 목록 자동 계산
GCA_VF_PTS = {s: _pts_in_vf_range(lo, hi) for s, (lo, hi, _) in GCA_VF_DEF.items()}


# ═══════════════════════════════════════════════════════════════
# 3.  RNFL 사분면 매핑 정의
# ═══════════════════════════════════════════════════════════════

# 상하 hemifield (검증됨: q_s ↔ inf VF r=0.44-0.51)
_VF_SUP = list(range(1, 28))   # y>0: p01-27 (27pts)
_VF_INF = list(range(28, 55))  # y<0: p28-54 (27pts)

# 중앙 VF (eccentricity ≤ 10°, q_t papillomacular bundle)
_VF_CENTRAL = [p for p in range(1,55)
               if math.sqrt(VF_COORDS[p][0]**2 + VF_COORDS[p][1]**2) <= 10.5]

# 주변 측두 VF (x ≥ +15°, q_n nasal disc → temporal field)
_VF_PERIPH_TEMP = [p for p in range(1,55) if VF_COORDS[p][0] >= 15]

RNFL_QUAD_VF = {
    'q_s': (_VF_INF,         'inferior VF       (하측 전체 p28-54)'),
    'q_i': (_VF_SUP,         'superior VF       (상측 전체 p01-27)'),
    'q_t': (_VF_CENTRAL,     'central VF        (황반 ≤10°)'),
    'q_n': (_VF_PERIPH_TEMP, 'peripheral temporal VF (x≥15°)'),
}


# ═══════════════════════════════════════════════════════════════
# 4.  RNFL clock-hour 매핑 정의
#     disc h01=60°…h12=90°, VF = disc+180°±15°
# ═══════════════════════════════════════════════════════════════

# disc 중심 각도 (OD 해부학 polar, 0°=temporal, 반시계=양수)
_DISC_CENTER = {
    'h01': 60, 'h02': 30,  'h03':  0,  'h04': 330, 'h05': 300, 'h06': 270,
    'h07':240, 'h08':210,  'h09':180,  'h10': 150, 'h11': 120, 'h12':  90,
}

CH_VF_DEF = {}
for h, dc in _DISC_CENTER.items():
    vf_center = (dc + 180) % 360
    lo = (vf_center - 15) % 360
    hi = (vf_center + 15) % 360
    CH_VF_DEF[h] = (lo, hi, f'VF center {vf_center}°')

CH_VF_PTS = {h: _pts_in_vf_range(lo, hi) for h, (lo, hi, _) in CH_VF_DEF.items()}


# ═══════════════════════════════════════════════════════════════
# 5.  역방향: VF점 → 담당 GCA 섹터 (출력용)
# ═══════════════════════════════════════════════════════════════

def _pt_to_gca(p):
    a = _vf_angle(p)
    for s, (lo, hi, _) in GCA_VF_DEF.items():
        if lo < hi:
            if lo <= a < hi: return s
        else:
            if a >= lo or a < hi: return s
    return '?'

VF_TO_GCA = {p: _pt_to_gca(p) for p in range(1, 55)}


# ═══════════════════════════════════════════════════════════════
# 6.  매핑 대응표 출력 (검토용)
# ═══════════════════════════════════════════════════════════════

def print_mapping_tables():
    print('='*70)
    print('Phase 5-2 매핑 대응표 (OD-unified, Garway-Heath 기반)')
    print('='*70)

    print('\n▶ GCA 6섹터 → VF 그룹 (60° 슬라이스, 180° 회전)')
    print(f'  {"섹터":8s}  {"VF 각도 범위":14s}  {"VF 영역":22s}  {"VF 점 수":6s}  VF 점 목록')
    print('  '+'-'*90)
    for s, (lo, hi, desc) in GCA_VF_DEF.items():
        pts = GCA_VF_PTS[s]
        print(f'  {s:8s}  [{lo:3d}°, {hi:3d}°)    {desc:22s}  {len(pts):3d}pts  {pts}')

    print(f'\n  ※ VF점별 담당 GCA 섹터 (54점 분류)::')
    for s in ['sup_t','sup','sup_n','inf_n','inf','inf_t']:
        pts = GCA_VF_PTS[s]
        coords = [(p, VF_COORDS[p]) for p in pts]
        print(f'    {s:6s}: {[(p,f"({x:+d},{y:+d})") for p,(x,y) in coords]}')

    print('\n▶ RNFL 4사분면 → VF 그룹')
    print(f'  {"사분면":6s}  {"VF 영역":38s}  {"VF 점 수":6s}')
    print('  '+'-'*60)
    for q, (pts, desc) in RNFL_QUAD_VF.items():
        print(f'  {q:6s}  {desc:38s}  {len(pts):3d}pts  {pts[:8]}{"..." if len(pts)>8 else ""}')

    print('\n▶ RNFL clock-hour → VF 그룹 (30° 슬라이스)')
    print(f'  {"clock":6s}  {"disc 중심":10s}  {"VF 각도 범위":16s}  {"pts":4s}  VF 점 목록')
    print('  '+'-'*80)
    for h, (lo, hi, desc) in CH_VF_DEF.items():
        dc = _DISC_CENTER[h]
        vc = (dc + 180) % 360
        pts = CH_VF_PTS[h]
        print(f'  {h:6s}  disc={dc:3d}°     VF [{lo:3d}°,{hi:3d}°)  ({vc}°±15°)  '
              f'{len(pts):2d}pts  {pts}')

    print('\n▶ VF 54점 → GCA 섹터 분류 (확인용)')
    sec_map = defaultdict(list)
    for p, s in VF_TO_GCA.items():
        sec_map[s].append(p)
    for s, pts in sorted(sec_map.items()):
        print(f'  {s:6s}: {pts}')

    print('\n▶ 해부학 대응 요약 (이 관계가 성립해야 함)')
    print('  GCA sup_t (황반 측두상방 망막) → inferonasal VF (하비측)  [상하+좌우 반전]')
    print('  GCA sup   (황반 상방 망막)    → inferior VF (하측 중앙)  [상하 반전]')
    print('  GCA inf   (황반 하방 망막)    → superior VF (상측 중앙)  [상하 반전]')
    print('  GCA inf_t (황반 측두하방 망막) → superonasal VF (상비측)  [상하+좌우 반전]')
    print('  RNFL q_s  (상방 디스크)       → inferior VF  (하측 전체) [검증됨]')
    print('  RNFL q_i  (하방 디스크)       → superior VF  (상측 전체) [검증됨]')
    print('  RNFL q_t  (측두 PMB)          → central VF   (황반 ≤10°)')
    print('  RNFL q_n  (비측 디스크)       → peripheral temporal VF  (x≥15°)')


# ═══════════════════════════════════════════════════════════════
# 7.  데이터 로드 + OS VF flip
# ═══════════════════════════════════════════════════════════════

def _f(v):
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except (TypeError, ValueError):
        return None

def flip_vf(vals54):
    """OS VF: 각 row 내 역순 (make_flip.py 동일 로직)"""
    out = list(vals54)
    for grp in _VF_FLIP_ROWS:
        idxs = [i - 1 for i in grp]
        sub  = [vals54[i] for i in idxs]
        for idx, val in zip(idxs, reversed(sub)):
            out[idx] = val
    return out

def load_data(cohort_filter=True):
    """analysis_master.csv → 187눈 분석 코호트 행들(기본).

    필터를 끄면 275눈 전체가 들어와 논문 표와 어긋난다(2026-08-15 감사).
    """
    rows = list(csv.DictReader(
        open(_ROOT / 'data' / 'analysis_master.csv', encoding='utf-8')))
    if cohort_filter:
        rows = filter_rows(rows)
    for r in rows:
        # VF 54점 (OS는 flip)
        raw = [_f(r.get(f'p{i}')) for i in range(1, 55)]
        r['_vf'] = flip_vf(raw) if r['eye'] == 'OS' else raw
        r['_vf_mean'] = _f(r.get('vf_mean'))
        r['_pid'] = r['patient_id']

        # GCA eye-specific
        e = r['eye'].lower()
        for s in GCA_VF_DEF:
            col = f'{e}_s_{s}'      # od_s_sup_t / os_s_sup_t 등
            r[f'_gcl_{s}'] = _f(r.get(col))
        r['_gcl_avg'] = _f(r.get(f'avg_gcl_{e}'))
        r['_gcl_min'] = _f(r.get(f'min_gcl_{e}'))

        # RNFL (eye-independent label)
        for q in ('s', 't', 'i', 'n'):
            r[f'_rnfl_q{q}'] = _f(r.get(f'rnfl_q_{q}'))
        for h in range(1, 13):
            r[f'_rnfl_h{h:02d}'] = _f(r.get(f'rnfl_h{h:02d}'))
    return rows


# ═══════════════════════════════════════════════════════════════
# 8.  VF 그룹 평균 감도 계산
# ═══════════════════════════════════════════════════════════════

def vf_group_mean(vf54, pt_list):
    vals = [vf54[p-1] for p in pt_list if vf54[p-1] is not None]
    return sum(vals) / len(vals) if vals else None


# ═══════════════════════════════════════════════════════════════
# 9.  상관 계산 유틸
# ═══════════════════════════════════════════════════════════════

def _fisher_ci(r, n):
    if n < 4 or abs(r) >= 1.0: return None, None
    z  = 0.5 * math.log((1+r)/(1-r))
    se = 1.0 / math.sqrt(n - 3)
    lo = math.tanh(z - 1.96*se)
    hi = math.tanh(z + 1.96*se)
    return round(lo, 4), round(hi, 4)

def corr_block(xs, ys, label, groups=None):
    pairs = [(x,y,g) for x,y,g in zip(xs, ys, groups or xs)
             if x is not None and y is not None]
    if len(pairs) < 5:
        return {'label': label, 'n': len(pairs), 'note': 'too few'}
    n  = len(pairs)
    xv = [p[0] for p in pairs]
    yv = [p[1] for p in pairs]
    gv = [p[2] for p in pairs]

    pr, pp = sps.pearsonr(xv, yv)
    sr, sp = sps.spearmanr(xv, yv)
    p_ci   = _fisher_ci(float(pr), n)
    s_ci   = _fisher_ci(float(sr), n)

    res = {
        'label':   label, 'n': n,
        'pearson':  {'r':   round(float(pr), 4), 'p': round(float(pp), 4), 'ci': list(p_ci)},
        'spearman': {'rho': round(float(sr), 4), 'p': round(float(sp), 4), 'ci': list(s_ci)},
    }

    if groups is not None and n >= 10:
        try:
            import pandas as pd
            import statsmodels.formula.api as smf
            df = pd.DataFrame({'x': xv, 'y': yv, 'pid': gv})
            df['zx'] = (df['x'] - df['x'].mean()) / df['x'].std()
            df['zy'] = (df['y'] - df['y'].mean()) / df['y'].std()
            m    = smf.mixedlm('zy ~ zx', df, groups=df['pid']).fit(reml=True)
            beta = float(m.params['zx'])
            se_b = float(m.bse['zx'])
            ci   = [float(m.conf_int().loc['zx', 0]), float(m.conf_int().loc['zx', 1])]
            res['lme'] = {
                'beta': round(beta, 4), 'se': round(se_b, 4),
                'p':    round(float(m.pvalues['zx']), 4),
                'ci':   [round(c, 4) for c in ci],
            }
        except Exception as e:
            res['lme'] = {'error': str(e)}
    return res


# ═══════════════════════════════════════════════════════════════
# 10.  FDR 보정 (Benjamini-Hochberg)
# ═══════════════════════════════════════════════════════════════

def fdr_bh(p_values):
    """p값 리스트 → q값 리스트 (BH)"""
    n   = len(p_values)
    idx = sorted(range(n), key=lambda i: p_values[i])
    q   = [None] * n
    min_q = 1.0
    for rank, i in enumerate(reversed(idx), 1):
        q_val = p_values[i] * n / (n - rank + 1)
        min_q = min(min_q, q_val)
        q[i]  = round(min(1.0, min_q), 5)
    return q


# ═══════════════════════════════════════════════════════════════
# 11.  메인 분석
# ═══════════════════════════════════════════════════════════════

def main():
    print_mapping_tables()

    print('\n\n' + '='*70)
    print('Phase 5-2 섹터별 구조-기능 상관 계산')
    print('='*70)

    rows = load_data()
    od   = [r for r in rows if r['eye'] == 'OD']
    os_  = [r for r in rows if r['eye'] == 'OS']
    all_ = rows
    pids = [r['_pid'] for r in all_]

    results = {
        'generated': datetime.date.today().isoformat(),
        'design': {
            'n_total': len(rows), 'n_od': len(od), 'n_os': len(os_),
            'vf_flip': 'OS 행 역순 (make_flip.py 동일)',
            'oct_flip': 'OS 섹터 그대로 (oct_tabular_90d 해부학 방향)',
            'mapping': 'Garway-Heath 2000 기반, 180° 회전, 60° 슬라이스 (GCA)',
            'vf_coords': {str(p): VF_COORDS[p] for p in range(1,55)},
        },
        'gcl_sector':     {},
        'rnfl_quad':      {},
        'rnfl_clock':     {},
        'cross_sector':   {},
    }

    # --- GCA 6섹터 × VF 그룹 (대각 + 교차) ---
    print('\n── GCA 섹터 상관 ──────────────────────────────────────')
    cross_gcl = {}
    for sec in GCA_VF_DEF:
        vf_pts     = GCA_VF_PTS[sec]
        struct_col = f'_gcl_{sec}'

        x_all = [r[struct_col] for r in all_]
        vf_corr = {}
        for other_sec in GCA_VF_DEF:
            yv = [vf_group_mean(r['_vf'], GCA_VF_PTS[other_sec]) for r in all_]
            tag = 'DIAG' if other_sec == sec else 'off'
            blk = corr_block(x_all, yv,
                             label=f'gcl_{sec} x vf_{other_sec}',
                             groups=pids if other_sec == sec else None)
            vf_corr[other_sec] = blk
        cross_gcl[sec] = vf_corr

        diag = cross_gcl[sec][sec]
        rho  = diag.get('spearman', {}).get('rho', '?')
        lme_b = diag.get('lme', {}).get('beta', '')
        lme_s = f'  LME β={lme_b:.3f}' if isinstance(lme_b, float) else ''
        print(f'  gcl_{sec:6s} ↔ vf_{sec:6s}: ρ={rho:6.3f}{lme_s}  (n={diag.get("n")})')
        results['gcl_sector'][sec] = vf_corr

    # --- RNFL 사분면 (대각 + 교차, GCA와 동일 구조) ---
    print('\n── RNFL 사분면 상관 (교차행렬) ──────────────────────────')
    quads = list(RNFL_QUAD_VF.keys())
    cross_rnfl = {}
    for q, (_, desc) in RNFL_QUAD_VF.items():
        xcol = f'_rnfl_q{q[-1]}'
        xv   = [r[xcol] for r in all_]
        vf_corr = {}
        for other_q, (other_pts, _) in RNFL_QUAD_VF.items():
            yv  = [vf_group_mean(r['_vf'], other_pts) for r in all_]
            blk = corr_block(xv, yv,
                             label=f'rnfl_{q} x vf_{other_q}',
                             groups=pids if other_q == q else None)
            vf_corr[other_q] = blk
        cross_rnfl[q] = vf_corr
        results['rnfl_quad'][q] = vf_corr[q]  # 대각(대응 VF)만 기존 키로 유지 — 하위호환

        diag = vf_corr[q]
        rho  = diag.get('spearman', {}).get('rho', '?')
        lme_b = diag.get('lme', {}).get('beta', '')
        lme_s = f'  LME β={lme_b:.3f}' if isinstance(lme_b, float) else ''
        print(f'  rnfl_{q} ↔ {desc[:30]:30s}: ρ={rho:6.3f}{lme_s}  (n={diag.get("n")})')

    print('\n── RNFL 사분면 교차상관 행렬 (Spearman ρ) ──────────────')
    header = f'  {"rnfl \\ vf":10s}' + ''.join(f'{q:8s}' for q in quads)
    print(header)
    for r_q in quads:
        row_str = f'  {r_q:8s}  '
        for c_q in quads:
            blk = cross_rnfl[r_q][c_q]
            rho = blk.get('spearman', {}).get('rho', float('nan'))
            marker = '*' if r_q == c_q else ' '
            row_str += f'{rho:7.3f}{marker}'
        print(row_str)
    print('  (* = 대각: 해부학 대응 사분면)')

    rnfl_diag_rhos    = [cross_rnfl[q][q].get('spearman', {}).get('rho', 0) for q in quads]
    rnfl_offdiag_rhos = [cross_rnfl[r][c].get('spearman', {}).get('rho', 0)
                         for r in quads for c in quads if r != c]
    results['cross_rnfl'] = {
        q: {oq: cross_rnfl[q][oq] for oq in quads} for q in quads
    }
    results['cross_rnfl']['diag_mean_rho']    = round(sum(rnfl_diag_rhos)/len(rnfl_diag_rhos), 4)
    results['cross_rnfl']['offdiag_mean_rho'] = round(sum(rnfl_offdiag_rhos)/len(rnfl_offdiag_rhos), 4)
    print(f'\n  대각(대응 사분면) 평균 ρ = {results["cross_rnfl"]["diag_mean_rho"]}')
    print(f'  비대각(비대응)   평균 ρ = {results["cross_rnfl"]["offdiag_mean_rho"]}')

    # --- RNFL clock-hour ---
    print('\n── RNFL clock-hour 상관 (개별) ──────────────────────────')
    for h in [f'h{i:02d}' for i in range(1, 13)]:
        xcol = f'_rnfl_{h}'
        vfp  = CH_VF_PTS[h]
        xv   = [r[xcol] for r in all_]
        yv   = [vf_group_mean(r['_vf'], vfp) for r in all_]
        blk  = corr_block(xv, yv, label=f'rnfl_{h} x vf_ch{h}', groups=pids)
        results['rnfl_clock'][h] = blk
        rho  = blk.get('spearman', {}).get('rho', '?')
        n    = blk.get('n', 0)
        lme_b = blk.get('lme', {}).get('beta', '')
        lme_s = f'  β={lme_b:.3f}' if isinstance(lme_b, float) else ''
        vc   = ((_DISC_CENTER[h] + 180) % 360)
        print(f'  {h}: VF {vc:3d}°  ρ={rho:6.3f}{lme_s}  (n={n}, vfpts={len(vfp)})')

    # --- FDR 보정 ---
    # GCA 대각 (6개) + RNFL 사분면 (4개) + clock-hour (12개) = 22 tests
    fdr_items = []
    for sec in GCA_VF_DEF:
        blk = results['gcl_sector'][sec][sec]
        fdr_items.append(('gcl_' + sec, blk.get('spearman', {}).get('p', 1.0)))
    for q in RNFL_QUAD_VF:
        blk = results['rnfl_quad'][q]
        fdr_items.append(('rnfl_' + q, blk.get('spearman', {}).get('p', 1.0)))
    for h in [f'h{i:02d}' for i in range(1,13)]:
        blk = results['rnfl_clock'][h]
        fdr_items.append(('rnfl_' + h, blk.get('spearman', {}).get('p', 1.0)))

    labels, pvals = zip(*fdr_items)
    qvals = fdr_bh(list(pvals))
    fdr_out = {lbl: {'p_raw': round(p, 5), 'q_bh': q}
               for lbl, p, q in zip(labels, pvals, qvals)}
    results['fdr_bh'] = fdr_out

    print('\n── FDR 보정 요약 ────────────────────────────────────────')
    print(f'  {"섹터":20s}  {"p_raw":8s}  {"q_BH":8s}  sig(q<0.05)')
    for lbl, p, q in zip(labels, pvals, qvals):
        sig = '***' if q < 0.001 else ('**' if q < 0.01 else ('*' if q < 0.05 else ''))
        print(f'  {lbl:20s}  {p:8.4f}  {q:8.4f}  {sig}')

    # --- Cross-sector 비교 (대각 vs 비대각) ---
    print('\n── GCA 섹터 교차상관 행렬 (Spearman ρ) ─────────────────')
    secs = list(GCA_VF_DEF.keys())
    header = f'  {"gcl \\ vf":10s}' + ''.join(f'{s:8s}' for s in secs)
    print(header)
    for r_sec in secs:
        row_str = f'  gcl_{r_sec:6s}'
        for c_sec in secs:
            blk = cross_gcl[r_sec][c_sec]
            rho = blk.get('spearman', {}).get('rho', float('nan'))
            marker = '*' if r_sec == c_sec else ' '
            row_str += f'{rho:7.3f}{marker}'
        print(row_str)
    print('  (* = 대각: 해부학 대응 섹터)')

    diag_rhos    = [cross_gcl[s][s].get('spearman',{}).get('rho', 0) for s in secs]
    offdiag_rhos = [cross_gcl[r][c].get('spearman',{}).get('rho', 0)
                    for r in secs for c in secs if r != c]
    results['cross_sector']['diag_mean_rho']    = round(sum(diag_rhos)/len(diag_rhos), 4)
    results['cross_sector']['offdiag_mean_rho'] = round(sum(offdiag_rhos)/len(offdiag_rhos), 4)
    print(f'\n  대각(대응 섹터) 평균 ρ = {results["cross_sector"]["diag_mean_rho"]}')
    print(f'  비대각(비대응)  평균 ρ = {results["cross_sector"]["offdiag_mean_rho"]}')
    print('  ※ Garway-Heath 매핑이 맞으면 대각 > 비대각')

    # --- 저장 ---
    out_path = _ROOT / 'results' / 'phase5_2_sector_correlations.json'
    results['input_sha256'] = sha256_file(str(_ROOT / 'data' / 'analysis_master.csv'))
    json.dump(results, open(out_path, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
    print(f'\n출력: {out_path}')


if __name__ == '__main__':
    main()
