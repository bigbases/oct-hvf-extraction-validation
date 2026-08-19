"""
scripts/make_oct_field_accuracy.py
Fig: OCT 필드별 추출 정확도 (VF 54점 히트맵의 OCT 대응물).

VF 쪽은 vf_grid_heatmap.pdf 가 54점 공간분포를 보여주는데 OCT 는 표뿐이라
modality 간 표현이 비대칭이었다. 같은 색상 규약(정확도 낮을수록 진한 색)으로
OCT 필드별 exact-match 를 원본 차트 배치 위에 얹는다.

패널:
  (a) GCIPL 6섹터 — OD눈 / OS눈 도넛 2개 병치
  (b) RNFL 4사분면
  (c) RNFL 12 clock-hour

좌표계(중요): 이 그림은 **on-screen 프레임**(기기가 출력한 원본 차트 배치)이다.
논문 본문의 구조-기능 분석은 OS를 OD로 미러링한 OD-normalized 값을 쓰지만,
여기서 보여주는 건 '추출이 화면의 어느 위치에서 실패했는가'이므로 미러링하면
그 위치 정보가 사라진다. GCA 6섹터 파이차트는 스키마틱 오버레이라 양안의 섹터
각도가 동일하다(src/hvf/ocr_oct.py::parse_sectors, angles_deg 고정, eye 분기 없음).
따라서 두 도넛은 같은 기하에 값만 다르다.

각도 정본 (전부 src/hvf/ocr_oct.py 에서 그대로 가져옴):
  GCA 6섹터   : sup 90, sup_t 30, inf_t 330, inf 270, inf_n 210, sup_n 150 (각 ±30)
  clock-hour  : h{n} = 90 - n*30 (각 ±15) → h12 위, h03 우측(temporal)
  RNFL 사분면 : S 90, T 0, I 270, N 180 (각 ±45) — clock 과 같은 프레임

수치 출처: results/phase3_extraction_accuracy.json (GT 풀 n=280눈). 새 계산 없음.
요약 평균은 포인트 풀링(가중) 기준 — 원고 Table 5 와 동일 정의
(clock-hour 50.4%, 사분면 94.8%).

출력: paper/manuscript/figures/oct_field_accuracy.{pdf,png}
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['font.family'] = 'sans-serif'
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Wedge
from matplotlib.cm import ScalarMappable

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))

REGISTRY = _ROOT / 'results' / 'phase3_extraction_accuracy.json'
OUT_DIR  = _ROOT / 'paper' / 'manuscript' / 'figures'

# ── 색상 규약: vf_grid_heatmap 과 동일 방향(오차율이 높을수록 진함) ──
# 실측 오차율 범위 0.7–78.3% 이므로 10%p 간격 8구간.
BOUNDS = [0, 10, 20, 30, 40, 50, 60, 70, 80]
CMAP   = plt.get_cmap('Reds', len(BOUNDS) - 1)
NORM   = BoundaryNorm(BOUNDS, CMAP.N)
# Reds 는 오차율 45%대부터 이미 상당히 진해 검은 글씨 대비가 떨어진다 — 45 기준.
WHITE_TEXT_ABOVE = 45.0

R_OUT, R_IN = 1.0, 0.42
LIM = 1.40          # 네 패널 공통 축 범위 — 도넛 지름을 동일하게 유지
# 이 논문의 다른 그림(fig_input_example_full)과 같은 8pt 로 통일. 캔버스가 본문
# 폭과 1:1 이라 지정한 pt 가 지면에서 그대로 8pt 다.
FS_VAL = FS_LAB = FS_PANEL = FS_CBAR = 8.0

# 섹터 중심각 (src/hvf/ocr_oct.py 정본)
GCA_ANGLES  = {'sup': 90, 'sup_t': 30, 'inf_t': 330, 'inf': 270, 'inf_n': 210, 'sup_n': 150}
GCA_ORDER   = ['sup', 'sup_t', 'inf_t', 'inf', 'inf_n', 'sup_n']
QUAD_ANGLES = {'s': 90, 't': 0, 'i': 270, 'n': 180}
QUAD_ORDER  = ['s', 't', 'i', 'n']


def load_accuracy():
    d = json.load(open(REGISTRY, encoding='utf-8'))
    gcl = d['gcl_accuracy']['per_column']
    rnf = d['rnfl_accuracy']['per_column']
    gca = {eye: {s: gcl[f'{eye}_s_{s}']['exact_match_pct'] for s in GCA_ORDER}
           for eye in ('od', 'os')}
    quad = {q: rnf[f'rnfl_q_{q}']['exact_match_pct'] for q in QUAD_ORDER}
    clock = {n: rnf[f'rnfl_h{n:02d}']['exact_match_pct'] for n in range(1, 13)}

    def pooled(cols):
        ex = nv = 0.0
        for c in cols:
            v = rnf[c]
            nv += v['n_valid']
            ex += v['n_valid'] * v['exact_match_pct'] / 100
        return 100 * ex / nv

    means = {
        'quad':  pooled([f'rnfl_q_{q}' for q in QUAD_ORDER]),
        'clock': pooled([f'rnfl_h{n:02d}' for n in range(1, 13)]),
    }
    return gca, quad, clock, means


def draw_donut(ax, values, angles, half_width, title=None):
    """values: {key: exact_match_pct}, angles: {key: 중심각}. 오차율로 착색."""
    for key, exact in values.items():
        err = 100.0 - exact
        c   = angles[key]
        ax.add_patch(Wedge((0, 0), R_OUT, c - half_width, c + half_width,
                           width=R_OUT - R_IN,
                           facecolor=CMAP(NORM(err)), edgecolor='white', linewidth=0.8))
        rad = (R_OUT + R_IN) / 2
        import math
        x = rad * math.cos(math.radians(c))
        y = rad * math.sin(math.radians(c))
        ax.text(x, y, f'{exact:.1f}', ha='center', va='center', fontsize=FS_VAL,
                color=('white' if err >= WHITE_TEXT_ABOVE else 'black'))
    # 네 도넛의 지름을 같게 하려면 축 범위가 모두 같아야 한다(equal aspect).
    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM)
    ax.set_aspect('equal'); ax.axis('off')
    if title:
        ax.set_title(title, fontsize=FS_LAB, pad=1)


def main():
    gca, quad, clock, means = load_accuracy()

    import math

    # 축 박스가 정사각형에 가깝도록 높이를 잡아야 도넛 주변 여백이 안 생긴다.
    fig = plt.figure(figsize=(6.93, 3.72))   # 17.6 cm 폭
    # clock-hour 는 웨지가 12개라 8pt 라벨을 담으려면 지름이 더 커야 한다.
    # 열 너비를 키워 (c) 도넛만 크게 잡는다(다른 패널과 지름이 달라지는 건 감수).
    WR = [1, 1, 1, 1.55]
    gs = fig.add_gridspec(1, 4, left=0.010, right=0.990, top=0.885, bottom=0.250,
                          width_ratios=WR, wspace=0.06)
    ax_od = fig.add_subplot(gs[0, 0])
    ax_os = fig.add_subplot(gs[0, 1])
    ax_q  = fig.add_subplot(gs[0, 2])
    ax_c  = fig.add_subplot(gs[0, 3])

    draw_donut(ax_od, gca['od'], GCA_ANGLES, 30, title='OD')
    draw_donut(ax_os, gca['os'], GCA_ANGLES, 30, title='OS')
    draw_donut(ax_q, quad, QUAD_ANGLES, 45)
    draw_donut(ax_c, {n: clock[n] for n in range(1, 13)},
               {n: 90 - n * 30 for n in range(1, 13)}, 15)

    # 방위 라벨 — (a),(b)만. (c)는 시계 번호가 방위를 대신하고, 문자와 겹친다.
    for ax in (ax_od, ax_os, ax_q):
        for x, y, t in ((0, 1.14, 'S'), (0, -1.14, 'I'), (1.14, 0, 'T'), (-1.14, 0, 'N')):
            ax.text(x, y, t, ha='center', va='center', fontsize=FS_VAL, color='0.35')
    for n in range(1, 13):
        a = math.radians(90 - n * 30)
        ax_c.text(1.24 * math.cos(a), 1.24 * math.sin(a), str(n),
                  ha='center', va='center', fontsize=FS_LAB, color='0.45')

    # 패널 라벨 — 각 그룹 중앙 정렬 (열 중심은 gridspec 기하에서 계산)
    tot = sum(WR) + 3 * 0.06 * (sum(WR) / 4)
    unit = (0.990 - 0.010) / tot
    cen, acc = [], 0.010
    for wr in WR:
        cen.append(acc + unit * wr / 2)
        acc += unit * wr + unit * 0.06 * (sum(WR) / 4)
    for x, text in (((cen[0] + cen[1]) / 2, '(a) GCIPL 6-sector'),
                    (cen[2], '(b) RNFL quadrant'),
                    (cen[3], '(c) RNFL clock-hour')):
        fig.text(x, 0.935, text, fontsize=FS_PANEL, fontweight='bold', ha='center')

    # 요약 수치 — 해당 도넛 바로 아래
    fig.text(cen[2], 0.198, f'pooled mean {means["quad"]:.1f}%',
             fontsize=FS_VAL, ha='center', color='0.30')
    fig.text(cen[3], 0.198, f'pooled mean {means["clock"]:.1f}%',
             fontsize=FS_VAL, ha='center', color='0.30')
    fig.text(cen[3], 0.118, 'clock-hour 3 = temporal',
             fontsize=FS_LAB, ha='center', color='0.45')

    # 공통 컬러바 — (a),(b) 아래 중앙
    cax = fig.add_axes([0.115, 0.088, 0.38, 0.042])
    sm = ScalarMappable(norm=NORM, cmap=CMAP); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation='horizontal', ticks=BOUNDS)
    cb.set_label('Error rate, 100 $-$ exact-match (%)', fontsize=FS_CBAR, labelpad=1.5)
    cb.ax.tick_params(labelsize=FS_CBAR, length=2, pad=1)
    cb.outline.set_linewidth(0.6)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # tight 크롭을 쓰면 캔버스가 잘려 본문 폭과의 1:1 대응이 깨지고, 결국 8pt 가
    # 지면에서 8pt 가 아니게 된다. 캔버스를 그대로 저장한다.
    for ext in ('pdf', 'png'):
        fig.savefig(OUT_DIR / f'oct_field_accuracy.{ext}', dpi=300)
    print(f'저장: {OUT_DIR / "oct_field_accuracy.pdf"}')
    print(f'  GCIPL OD 범위 {min(gca["od"].values()):.1f}-{max(gca["od"].values()):.1f}%  '
          f'OS 범위 {min(gca["os"].values()):.1f}-{max(gca["os"].values()):.1f}%')
    print(f'  RNFL quadrant pooled mean {means["quad"]:.2f}%  '
          f'clock-hour pooled mean {means["clock"]:.2f}%')


if __name__ == '__main__':
    main()
