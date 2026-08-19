"""STAGE 30 — 데이터셋/registry → 그림 (paper/figures).

그림에 찍히는 수치는 registry(results/*.json) 값을 읽어 쓴다 (그림-표 불일치 방지).
현재 구현: Fig. VF-grid-heatmap (24-2 점별 결손률/정확도), Fig. SF-heatmap
(구조-기능 상관 비대칭 2-패널). 문서 내 등장 순서(ch4_results.tex에서 VF-grid를
먼저 \\input/\\includegraphics)로 LaTeX가 Figure 2/3을 자동 번호 매김.

VF-grid-heatmap: 24-2 격자(HVF_COORDS, kcc/analysis_ab.py의 검증된 배치 재사용)
위에 phase3_extraction_accuracy.json의 vf_accuracy.per_point_normal 점별 결손률/
exact-match를 얹는다. p27/p36만 뚜렷한 결손률(47.5%)로 보여 "추출 실패가 아니라
맹점 인접의 구조적 부재"를 그림으로 증명 — 기여2(검증 프로토콜+오류 taxonomy)의
시각물이 없었던 구멍을 메움.

SF-heatmap 구조 (구조-기능 상관 비대칭 2-패널):
  왼쪽  — RNFL 사분면 대응쌍(matched-pair)만, 1열 히트맵 (4행: S/T/I/N)
  오른쪽 — GCIPL 6섹터 교차상관 풀 매트릭스 (6x6)

비대칭인 이유(2026-07-20 확정, ch5_discussion.tex Limitations 참고):
  RNFL의 4개 사분면-매핑 VF 영역은 균등분할이 아니라 크기가 다른 개별 ROI
  (27/27/12/10점) — 풀 4x4 매트릭스를 그리면 좁은 중심영역(q_t, 12점)이
  전반적 중증도와 정의상 얽혀 대각선 특이성이 아티팩트로 흐려짐. GCIPL 6섹터는
  60도씩 균등 분할이라 6x6 대각/비대각 비교가 공정함. 두 구조가 다르므로
  시각 형식도 다르게 — RNFL은 검증된 대응쌍만, GCIPL은 풀 매트릭스.

  버리지 않은 전체 RNFL 4x4 교차행렬은 results/phase5_2_sector_correlations.json
  의 'cross_rnfl' 키에 보존되어 있다 (근거: 위 문서의 Limitations 문단).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))

from hvf.config import set_seed, results_dir

STAGE = "stage_30_figures"

# GCIPL 6섹터 표시 순서 (시계방향 피자 조각, 논문 tab:octacc/text 순서와 동일)
GCL_ORDER = ['sup_t', 'sup', 'sup_n', 'inf_n', 'inf', 'inf_t']
GCL_LABELS = {
    'sup_t': 'ST', 'sup': 'S', 'sup_n': 'SN',
    'inf_n': 'IN', 'inf': 'I', 'inf_t': 'IT',
}

# RNFL 사분면 표시 순서 (S/T/I/N, tab:octacc 순서와 동일)
RNFL_ORDER = ['q_s', 'q_t', 'q_i', 'q_n']
RNFL_LABELS = {'q_s': 'S', 'q_t': 'T', 'q_i': 'I', 'q_n': 'N'}


# 24-2 표준 시야좌표 (row, col), 0-indexed, OD 프레임, 10칸 캔버스(-27~+27, 6도 간격).
# col: 0=-27 1=-21 2=-15 3=-9 4=-3 5=+3 6=+9 7=+15 8=+21 9=+27.
# 중앙 두 행은 비측 -27도 추가점 때문에 9칸(col 0-8), 측두 +27도(col 9)는 24-2에 없음.
# 생리적 암점(col 7=+15도, row 4)은 '검사되는 점'이라 갭이 아님(OD GT 33/172눈 0dB로 검증).
# (구버전은 col 5에 갭을 두어 blind spot 위치를 +3도로 잘못 표시했음 — 그림 전용, 수치 무영향.)
HVF_COORDS = [
    (0,3),(0,4),(0,5),(0,6),
    (1,2),(1,3),(1,4),(1,5),(1,6),(1,7),
    (2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(2,7),(2,8),
    (3,0),(3,1),(3,2),(3,3),(3,4),(3,5),(3,6),(3,7),(3,8),
    (4,0),(4,1),(4,2),(4,3),(4,4),(4,5),(4,6),(4,7),(4,8),
    (5,1),(5,2),(5,3),(5,4),(5,5),(5,6),(5,7),(5,8),
    (6,2),(6,3),(6,4),(6,5),(6,6),(6,7),
    (7,3),(7,4),(7,5),(7,6),
]
assert len(HVF_COORDS) == 54


def _load_correlations():
    path = results_dir() / 'phase5_2_sector_correlations.json'
    return json.load(open(path, encoding='utf-8'))


def _load_vf_accuracy():
    path = results_dir() / 'phase3_extraction_accuracy.json'
    d = json.load(open(path, encoding='utf-8'))
    return d['vf_accuracy']['per_point_normal']


def make_vf_grid_heatmap(out_dir: Path):
    pp = _load_vf_accuracy()
    n_rows, n_cols = 8, 10

    miss_grid  = [[float('nan')] * n_cols for _ in range(n_rows)]
    exact_grid = [[float('nan')] * n_cols for _ in range(n_rows)]
    for idx, (r, c) in enumerate(HVF_COORDS):
        stats = pp[f'p{idx+1:02d}']
        miss_grid[r][c]  = stats['missing_rate_raw_pct']
        exact_grid[r][c] = stats['exact_match_pct']

    import numpy as np
    # 단일 패널: exact-match 를 오류율(100-exact)로 착색. missing 패널은 제거
    # (54칸 중 47칸이 0.0 → 정보량 없음; 본문에 한 줄로 대체). 셀 숫자는 exact-match(%).
    # 순차형 이산 컬러맵(Reds + BoundaryNorm): 정상(99%+)은 거의 흰색이라 가운데 두 행
    # 패턴이 즉시 보이고, 1~2건 차이가 색을 안 바꾸며(n≈320, 1%p≈3건), 흰 배경 위 검은 글씨.
    from matplotlib.colors import BoundaryNorm
    err_grid = [[float('nan')] * n_cols for _ in range(n_rows)]
    for idx, (r, c) in enumerate(HVF_COORDS):
        err_grid[r][c] = 100.0 - exact_grid[r][c]
    err_arr = np.ma.masked_invalid(err_grid)

    bounds = [0, 1, 2, 3, 4, 5, 6]           # 오류율(%) 구간. 실측 최대 5.3%라 상한 6.
    cmap = plt.get_cmap('Reds', len(bounds) - 1).copy()
    cmap.set_bad(color='#00000000')  # 격자 밖 = 완전 투명. 회색 배경은 칸칸이 나뉘면 그 자체가
    # 값처럼 읽혀 24-2 특유의 십자 실루엣(비측 돌출부 포함)이 사각형으로 뭉개짐 — 투명이 정답.
    norm = BoundaryNorm(bounds, cmap.N)

    valid = ~np.isnan(np.array(err_grid, dtype=float))

    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    im = ax.imshow(err_arr, cmap=cmap, norm=norm, aspect='equal')
    # 데이터 셀에만 얇은 테두리(격자 밖에는 그리지 않음) — 옅은 셀들이 뭉쳐 보이지 않고
    # HFA 프린트아웃처럼 칸칸이 구분됨. 굵으면 선이 색보다 강해져 스프레드시트로 보이므로 lw=0.8.
    for r in range(n_rows):
        for c in range(n_cols):
            if valid[r, c]:
                ax.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1,
                                           fill=False, ec='#333333', lw=0.8))
    ax.set_facecolor('white')
    ax.axis('off')  # 바깥 검은 사각 테두리 제거 — 격자가 사각형이 아니므로

    cbar = fig.colorbar(im, ax=ax, ticks=bounds, shrink=0.8,
                        label='Error rate, 100$-$exact-match (%)')
    # 그림 내부 제목 없음 — LaTeX 캡션이 제목 역할(중복 제거).

    # 셀 텍스트 = exact-match(%), 소수 첫째 자리. 오류율 구간에 따라 글자색 대비 확보.
    for idx, (r, c) in enumerate(HVF_COORDS):
        ev = exact_grid[r][c]
        err = 100.0 - ev
        tc = 'white' if err >= 4 else 'black'    # 최상위 구간(4-6%)만 흰 글씨
        ax.text(c, r, f'{ev:.1f}', ha='center', va='center', fontsize=8, color=tc)

    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png'):
        out_path = out_dir / f'vf_grid_heatmap.{ext}'
        fig.savefig(out_path, bbox_inches='tight', dpi=300)
        print(f'저장: {out_path}')
    plt.close(fig)


def make_sf_heatmap(out_dir: Path):
    d = _load_correlations()

    # 왼쪽 패널: RNFL 대응쌍(diagonal)만, (4,1) 열벡터
    rnfl_vals = [d['rnfl_quad'][q]['spearman']['rho'] for q in RNFL_ORDER]

    # 오른쪽 패널: GCIPL 6x6 풀 매트릭스
    gcl_matrix = [[d['gcl_sector'][r][c]['spearman']['rho'] for c in GCL_ORDER]
                  for r in GCL_ORDER]

    vmin, vmax = 0.0, 0.6  # 두 패널 공유 스케일
    cmap = plt.get_cmap('YlOrRd')
    norm = Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(9.5, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 6, 0.25], wspace=0.35)

    # --- 왼쪽: RNFL 대응쌍 (4x1) ---
    ax_l = fig.add_subplot(gs[0, 0])
    col = [[v] for v in rnfl_vals]
    ax_l.imshow(col, cmap=cmap, norm=norm, aspect='auto')
    for i, v in enumerate(rnfl_vals):
        txt_color = 'white' if v > (vmin + vmax) / 2 else 'black'
        ax_l.text(0, i, f'{v:.3f}', ha='center', va='center',
                   color=txt_color, fontsize=11, fontweight='bold')
    ax_l.set_yticks(range(len(RNFL_ORDER)))
    ax_l.set_yticklabels([RNFL_LABELS[q] for q in RNFL_ORDER], fontsize=11)
    ax_l.set_xticks([0])
    ax_l.set_xticklabels(['matched\nVF region'], fontsize=8)
    ax_l.set_title('RNFL quadrant\n(matched-pair only)', fontsize=10)

    # --- 오른쪽: GCIPL 6x6 풀 매트릭스 ---
    ax_r = fig.add_subplot(gs[0, 1])
    ax_r.imshow(gcl_matrix, cmap=cmap, norm=norm, aspect='auto')
    n = len(GCL_ORDER)
    for i in range(n):
        for j in range(n):
            v = gcl_matrix[i][j]
            txt_color = 'white' if v > (vmin + vmax) / 2 else 'black'
            marker = '*' if i == j else ''
            ax_r.text(j, i, f'{v:.3f}{marker}', ha='center', va='center',
                      color=txt_color, fontsize=9)
    ax_r.set_xticks(range(n))
    ax_r.set_xticklabels([GCL_LABELS[s] for s in GCL_ORDER], fontsize=10)
    ax_r.set_yticks(range(n))
    ax_r.set_yticklabels([GCL_LABELS[s] for s in GCL_ORDER], fontsize=10)
    ax_r.set_xlabel('VF region', fontsize=9)
    ax_r.set_ylabel('GCIPL sector', fontsize=9)
    ax_r.set_title('GCIPL sector × VF region (full 6×6 cross-correlation)',
                    fontsize=10)

    # --- 공유 컬러바 ---
    ax_cb = fig.add_subplot(gs[0, 2])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=ax_cb, label='Spearman ρ')

    fig.suptitle('Structure–function correlation: RNFL topographic specificity '
                 'vs. GCIPL global pattern', fontsize=11)

    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png'):
        out_path = out_dir / f'sf_heatmap.{ext}'
        fig.savefig(out_path, bbox_inches='tight', dpi=300)
        print(f'저장: {out_path}')
    plt.close(fig)


def main() -> None:
    set_seed()
    out_dir = _ROOT / 'paper' / 'figures'
    make_vf_grid_heatmap(out_dir)
    make_sf_heatmap(out_dir)
    print(f"[{STAGE}] 완료 - vf_grid_heatmap/sf_heatmap.{{pdf,png}} 생성.")


if __name__ == "__main__":
    main()
