"""
scripts/make_failure_cases.py
Fig: 추출 실패 사례 3건 — 수정가능성 기반 오류 분류의 시각 자료.

fig_input_example.pdf(Figure 2)는 성공한 크롭만 보여준다. 이 그림은 그 뒤에
(d)~(f)로 이어붙일 수 있는 형식으로 실패 사례를 보여준다(배치 방침 미확정이라
독립 파일로도 저장).

시각 언어(세 패널 공통):
  · 색은 그림 안의 사각형에만 쓴다 — 빨강 = 잘못된 좌표, 초록 = 올바른 좌표.
    텍스트는 전부 검정(흑백 인쇄 대비). 범례는 그림 안 색상 키로 넣는다.
  · 도형은 사각형. 단 (f)의 "읽히지 못한 값"은 좌표가 아니므로 색 없는 점선 표시.
  · 각 패널 3단 — (1) 좌표가 가리킨 위치, (2) 그 좌표에서 읽힌 픽셀,
    (3) extracted / reference 대비.
  · 폰트 크기는 기존 그림(stage_30_figures.py, make_sf_matrix.py)에 맞춘다:
    제목 10, 본문 라벨 8~9, 주석 8. font.family 는 지정하지 않는다(기존과 동일).

좌표는 전부 커밋된 코드/설정에서 가져온다 — 재구성하지 않는다:
  (d) 빨강 = 옛 parse_quadrants 탐색 지점(ocr_rnfl_detail.py, 커밋 b930e92),
      초록 = 현재 src/hvf/ocr_rnfl.py 의 셀 창.
  (e) config ocr.gca.od_sectors + ocr_oct.parse_sectors 의 270도 셀(지금도 실패).
  (f) config ocr.sfa.grid_crop 기준 격자 열(빨강, 실제로 읽은 곳).

개인정보: 리포트 기기 UI 한글 라벨을 마스킹한다(국제 저널 제출본).
차트·격자 영역 자체에는 식별정보가 없다.

입력: data/_failure_case_candidates.csv (환자ID·경로 — gitignore 대상)
출력: paper/manuscript/figures/fig_failure_cases.{pdf,png}
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
import numpy as np
import yaml
from PIL import Image, ImageDraw

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))

CAND = _ROOT / 'data' / '_failure_case_candidates.csv'
OUT_DIR = _ROOT / 'paper' / 'manuscript' / 'figures'
CFG = yaml.safe_load(open(_ROOT / 'config' / 'params.yaml', encoding='utf-8'))['ocr']

C_FAIL, C_FIX, C_NEUTRAL = '#d62728', '#2ca02c', '0.25'
# 1단 이미지의 종횡비 — 열 너비를 여기에 비례시켜야 패널 높이가 서로 맞는다.
PANEL_ASPECT = {'rnfl': 522 / 315, 'gcipl': 216 / 130, 'vf': 186 / 56}
# fig_input_example.pdf 와 같은 캔버스 폭(806pt)·같은 글자 크기(DejaVuSans 9pt)를
# 쓴다. 두 그림을 세로로 이어 붙인 뒤 본문 폭에 맞춰 실으면 같은 배율(0.619)로
# 축소되므로 지면에서 글자 크기가 일치한다.
# 캔버스를 본문 폭(17.6cm)과 1:1 로 둔다. 배율이 1 이므로 지정한 pt 가 그대로
# 지면 pt 가 된다 — 806pt 캔버스를 쓰면 0.62배로 축소돼 9pt 가 5.6pt 로 보였다.
PAGE_W_IN = 17.6 / 2.54
FS_TITLE, FS_LAB, FS_NOTE = 8, 8, 8
LW = 1.4

# (d) RNFL 사분면 패널 — Display item 1안 확정(2026-08-16)으로 제외됐다.
# "값을 못 읽어 결측"이라 결과만으로 이해되고, (e)와 같은 implementation defect
# 범주라 중복이기 때문. 코드는 지우지 않고 이 플래그로만 끈다.
INCLUDE_PANEL_D = False

# ── (d) 근거 ──────────────────────────────────────────────────────
# 결함은 크롭 사각형이 아니라 크롭 안의 셀 좌표였다(커밋 98eb992:
# "parse_quadrants 좌표버그 수정(B-scan 영역을 읽고 있었음)").
# 옛 비율은 프로덕션 크롭 635x201(차트 경계) 기준인데 현재 크롭은 config
# ocr.rnfl.quadrants = 889x311(차트보다 250px 왼쪽에서 시작)이라, 좁은 크롭용
# 비율을 넓은 크롭에 적용하면 네 지점이 왼쪽으로 밀린다.
#   교차검증: 옛 S(0.21x635)=페이지 x723, 새 S(0.430x889)=페이지 x722.
# 실측: 넷 중 셋(T/S/I)은 B-scan 패널 위, 넷째(N)는 차트 배경(숫자 없는 여백)에
# 떨어진다 — 그래서 네 값이 전부 결측이 된다.
OLD_QUAD_FRAC = {'S': (0.21, 0.13), 'T': (0.040, 0.55),
                 'I': (0.21, 0.92), 'N': (0.378, 0.55)}
NEW_QUAD_FRAC = {'S': (0.430, 0.088), 'T': (0.317, 0.413),
                 'I': (0.430, 0.736), 'N': (0.542, 0.413)}
CELL_HW, CELL_HH = 28, 18          # src/hvf/ocr_rnfl.py::parse_quadrants
QUAD_REFERENCE = 'S 98, T 59, I 98, N 52'
QUAD_NOTE = ('Left two crops: adjacent B-scan panel,\n'
             'not analysed in this study. Three points\n'
             'land there, the fourth on chart background.')

# ── (f) 근거 ──────────────────────────────────────────────────────
# 기본 grid_crop 으로 읽으면 3·4행이 한 칸씩 밀린다. 아래 값은 보정 패스를 끄고
# ocr_threshold_img() 를 직접 돌려 얻은 실측치(재현: _validate_gridcrop_fix.py).
VF_READ = ['31', '34', '33', '35']        # 빨간 사각형 4곳에서 읽힌 값
VF_TRUE = ['29', '31', '34', '33']        # 같은 자리의 참값 — 한 칸씩 밀렸다

# 기기 UI 한글 라벨 위치(페이지 좌표, RNFL 리포트) — 확대 시 노출되므로 가린다.
KO_MASK_RNFL = [(875, 1420, 965, 1478), (872, 1550, 968, 1625),
                (335, 1640, 495, 1675), (855, 1695, 970, 1800),
                (320, 1385, 500, 1425)]

# VF 리포트 원본 경로는 환자ID와 로컬 경로를 담고 있어 코드에 두지 않는다.
# 다른 패널과 같이 data/_failure_case_candidates.csv(=gitignore 대상)의 라벨로만
# 참조한다. 이 파일이 없으면 load_candidates() 가 시끄럽게 실패한다.
VF_CASE_LABEL = 'F1'


def load_candidates():
    if not CAND.exists():
        raise FileNotFoundError(f'후보 파일이 없다: {CAND}')
    return {r['label']: r for r in csv.DictReader(open(CAND, encoding='utf-8-sig'))}


def ratio_box(im, ratio):
    x1r, y1r, x2r, y2r = ratio
    return (int(im.width * x1r), int(im.height * y1r),
            int(im.width * x2r), int(im.height * y2r))


def bbox(cx, cy, hw, hh):
    return (cx - hw, cy - hh, cx + hw, cy + hh)


def union_view(boxes, margin):
    xs = [v for b in boxes for v in (b[0], b[2])]
    ys = [v for b in boxes for v in (b[1], b[3])]
    return (int(min(xs) - margin), int(min(ys) - margin),
            int(max(xs) + margin), int(max(ys) + margin))


def _bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def show(ax, im, view, boxes=(), ovals=()):
    """1단: 좌표가 가리킨 위치. boxes=[(page_box,color)], ovals=[page_box] (중립 점선)."""
    x0, y0, x1, y1 = view
    ax.imshow(np.asarray(im.crop(view)))
    for b in ovals:
        ax.add_patch(Ellipse(((b[0] + b[2]) / 2 - x0, (b[1] + b[3]) / 2 - y0),
                             (b[2] - b[0]) * 1.35, (b[3] - b[1]) * 1.35,
                             fill=False, edgecolor=C_NEUTRAL, linewidth=1.0,
                             linestyle=(0, (1.6, 1.6))))
    for b, color in boxes:
        ax.add_patch(Rectangle((b[0] - x0, b[1] - y0), b[2] - b[0], b[3] - b[1],
                               fill=False, edgecolor=color, linewidth=LW))
    ax.set_xlim(0, x1 - x0); ax.set_ylim(y1 - y0, 0)
    ax.set_anchor('N')
    _bare(ax)
    ax.add_patch(Rectangle((0, 0), x1 - x0, y1 - y0, fill=False,
                           edgecolor='0.6', linewidth=0.5))


def show_segments(axes, im, segs, y0, y1, boxes):
    """1단(분할판): x 구간 여러 개를 잘라 나란히 보여준다.

    (d) 는 잘못된 좌표 셋이 B-scan 안에 흩어져 있어 한 장으로 담으면 분석 대상이
    아닌 B-scan 이 지면의 절반을 먹는다. 아무것도 없는 구간을 빼고 이어 붙인다.
    """
    for ax, (sx0, sx1) in zip(axes, segs):
        ax.imshow(np.asarray(im.crop((sx0, y0, sx1, y1))))
        for b, color in boxes:
            if b[0] >= sx0 - 1 and b[2] <= sx1 + 1:
                ax.add_patch(Rectangle((b[0] - sx0, b[1] - y0), b[2] - b[0], b[3] - b[1],
                                       fill=False, edgecolor=color, linewidth=LW))
        ax.set_xlim(0, sx1 - sx0); ax.set_ylim(y1 - y0, 0)
        ax.set_anchor('N')
        _bare(ax)
        ax.add_patch(Rectangle((0, 0), sx1 - sx0, y1 - y0, fill=False,
                               edgecolor='0.6', linewidth=0.5))


def strip(ax, im, boxes, under=None, gap=7):
    """2단: 그 좌표에서 실제로 읽힌 픽셀. under 를 주면 각 조각 아래에 참값을 세로로 붙인다."""
    crops = [im.crop(b) for b in boxes]
    h = max(c.height for c in crops)
    w = sum(c.width for c in crops) + gap * (len(crops) - 1)
    canvas = Image.new('RGB', (w, h), (255, 255, 255))
    x, centers = 0, []
    for c in crops:
        canvas.paste(c, (x, (h - c.height) // 2))
        centers.append(x + c.width / 2)
        x += c.width + gap
    ax.imshow(np.asarray(canvas))
    pad = h * 1.05 if under else 0
    ax.set_xlim(0, w); ax.set_ylim(h + pad, 0)
    ax.set_anchor('N')
    _bare(ax)
    ax.add_patch(Rectangle((0, 0), w, h, fill=False, edgecolor=C_FAIL, linewidth=LW))
    if under:
        for cx, txt in zip(centers, under):
            ax.text(cx, h + pad * 0.62, txt, ha='center', va='center', fontsize=FS_LAB)


# ── 패널 기하 (그리기와 분리 — 결합 그림에서도 같은 좌표를 쓴다) ──────
def gcipl_geometry(cands):
    """반환: (이미지, 1단 view, 실패한 셀 창)"""
    im = Image.open(cands['E6']['image_path']).convert('RGB')
    sec = ratio_box(im, CFG['gca']['od_sectors'])
    sw, sh = sec[2] - sec[0], sec[3] - sec[1]
    win = bbox(sec[0] + sw // 2, sec[1] + int(sh // 2 + min(sw, sh) * 0.35), 55, 20)
    return im, (556, 1120, 772, 1250), win


def vf_geometry(cands=None):
    """반환: (이미지, 1단 view, 실제로 읽은 셀 4개, 밀려나 못 읽은 셀)"""
    from hvf.ocr_vf import _ROW34_OD_SHIFT_PX, get_grid_info
    cands = cands or load_candidates()
    im = Image.open(cands[VF_CASE_LABEL]['image_path']).convert('RGB')
    base = ratio_box(im, CFG['sfa']['grid_crop'])
    row_c, col_c, _, _ = get_grid_info(base[2] - base[0], base[3] - base[1])
    y3 = base[1] + row_c[3]
    read_b = [bbox(base[0] + col_c[i], y3, 15, 17) for i in range(4)]
    lost = bbox(base[0] - _ROW34_OD_SHIFT_PX + col_c[0], y3, 15, 17)
    # 좌측 파란 표식과 우측의 잘린 숫자를 제외한 좁은 범위.
    view = (lost[0] - 5, y3 - 28, read_b[3][2] + 5, y3 + 28)
    return im, view, read_b, lost


def draw_panel(fig, spec, kind, cands=None):
    """spec(SubplotSpec) 안에 1단+2단을 그린다. 반환: (제목, extracted, reference, 비고)"""
    sub = spec.subgridspec(2, 1, height_ratios=[1.0, 0.34], hspace=0.10)
    if kind == 'gcipl':
        im, view, win = gcipl_geometry(cands)
        show(fig.add_subplot(sub[0]), im, view, [(win, C_FAIL)])
        strip(fig.add_subplot(sub[1]), im, [win])
        return 'GCIPL 6-sector, inferior cell', '8', '69', None
    if kind == 'vf':
        im, view, read_b, lost = vf_geometry(cands)
        show(fig.add_subplot(sub[0]), im, view,
             [(b, C_FAIL) for b in read_b], ovals=[lost])
        strip(fig.add_subplot(sub[1]), im, read_b, under=VF_TRUE)
        return ('VF 24-2 grid, row 3',
                ', '.join(VF_READ) + ', …', ', '.join(VF_TRUE) + ', …', None)
    raise ValueError(kind)


def build():
    cands = load_candidates()
    # 패널 순서: INCLUDE_PANEL_D 가 켜져 있으면 RNFL 패널이 맨 앞에 붙고 라벨이
    # 한 칸씩 밀린다. 꺼져 있으면 GCIPL=(d), VF=(e).
    n_panels = 3 if INCLUDE_PANEL_D else 2
    letters = ['(d)', '(e)', '(f)'] if INCLUDE_PANEL_D else ['(d)', '(e)']
    # 열 너비는 각 1단 이미지의 종횡비에 비례시켜야 표시 높이가 같아진다.
    ASPECTS = ([522 / 315] if INCLUDE_PANEL_D else []) + [PANEL_ASPECT['gcipl'], PANEL_ASPECT['vf']]

    fig = plt.figure(figsize=(PAGE_W_IN, 3.5))
    gs = fig.add_gridspec(2, n_panels, height_ratios=[1.0, 0.46],
                          width_ratios=ASPECTS,
                          left=0.100, right=0.900, top=0.900, bottom=0.190,
                          hspace=0.05, wspace=0.07)
    cols = [gs[0, i].get_position(fig).x0 for i in range(n_panels)]

    def title(i, text):
        # 기존 (a)(b)(c) 는 DejaVuSans 9pt 일반 굵기 — 같은 사양으로 맞춘다.
        fig.text(cols[i], 0.935, f'{letters[i]} {text}', fontsize=FS_TITLE, ha='left')

    def caption(i, ex, ref, note=None):
        ax = fig.add_subplot(gs[1, i]); ax.axis('off')
        y = 0.86 if note else 0.72
        for lab, val in (('extracted:', ex), ('reference:', ref)):
            ax.text(0.49, y, lab, fontsize=FS_LAB, fontweight='bold',
                    ha='right', va='center', color='black', transform=ax.transAxes)
            ax.text(0.53, y, val, fontsize=FS_LAB, ha='left', va='center',
                    color='black', transform=ax.transAxes)
            y -= 0.26
        if note:
            ax.text(0.5, y + 0.06, note, fontsize=FS_NOTE - 0.8, ha='center', va='top',
                    color='black', transform=ax.transAxes, linespacing=1.5)

    col = 0   # 다음에 그릴 열
    # ══ RNFL 사분면 셀 좌표가 B-scan 패널을 가리킨 사례 (기본 제외) ══
    if INCLUDE_PANEL_D:
      im_d = Image.open(cands['D1']['image_path']).convert('RGB')
      dd = ImageDraw.Draw(im_d)
      for b in KO_MASK_RNFL:
          dd.rectangle(b, fill=(255, 255, 255))
      cb = ratio_box(im_d, CFG['rnfl']['quadrants'])
      cw, ch = cb[2] - cb[0], cb[3] - cb[1]
      old_b = [bbox(cb[0] + cw * fx, cb[1] + ch * fy, CELL_HW, CELL_HH)
               for fx, fy in OLD_QUAD_FRAC.values()]
      new_b = [bbox(cb[0] + cw * fx, cb[1] + ch * fy, CELL_HW, CELL_HH)
               for fx, fy in NEW_QUAD_FRAC.values()]
      sub = gs[0, col].subgridspec(2, 1, height_ratios=[1.0, 0.34], hspace=0.10)
      # 세 구간으로 잘라 붙인다: 잘못된 좌표가 앉은 B-scan 두 곳 + 차트 본체.
      # 사이의 빈 B-scan 을 빼면 차트가 패널의 주인공이 된다(36% vs 64%).
      dy0, dy1 = 1427, 1742
      segs = [(337, 413), (488, 564), (584, 862)]
      seg_gs = sub[0].subgridspec(1, 3, width_ratios=[t[1] - t[0] for t in segs],
                                  wspace=0.055)
      show_segments([fig.add_subplot(seg_gs[0, i]) for i in range(3)], im_d, segs,
                    dy0, dy1,
                    [(b, C_FAIL) for b in old_b] + [(b, C_FIX) for b in new_b])
      strip(fig.add_subplot(sub[1]), im_d, old_b)
      title(col, 'RNFL quadrant chart')
      caption(col, '—  (all four missing)', QUAD_REFERENCE, QUAD_NOTE)
      col += 1

    # ══ (e) 하부 GCIPL 섹터 창이 숫자 하단만 걸친 사례 ═══════════════
    im_e = Image.open(cands['E6']['image_path']).convert('RGB')
    sec = ratio_box(im_e, CFG['gca']['od_sectors'])
    sw, sh = sec[2] - sec[0], sec[3] - sec[1]
    win = bbox(sec[0] + sw // 2, sec[1] + int(sh // 2 + min(sw, sh) * 0.35), 55, 20)
    sub = gs[0, col].subgridspec(2, 1, height_ratios=[1.0, 0.34], hspace=0.10)
    show(fig.add_subplot(sub[0]), im_e, (556, 1120, 772, 1250), [(win, C_FAIL)])
    strip(fig.add_subplot(sub[1]), im_e, [win])
    title(col, 'GCIPL 6-sector, inferior cell')
    caption(col, '8', '69')
    col += 1

    # ══ (f) 격자 열이 한 칸 밀린 사례 ═══════════════════════════════
    from hvf.ocr_vf import _ROW34_OD_SHIFT_PX, get_grid_info
    im_f = Image.open(cands[VF_CASE_LABEL]['image_path']).convert('RGB')
    base = ratio_box(im_f, CFG['sfa']['grid_crop'])
    row_c, col_c, _, _ = get_grid_info(base[2] - base[0], base[3] - base[1])
    y3 = base[1] + row_c[3]
    read_b = [bbox(base[0] + col_c[i], y3, 15, 17) for i in range(4)]
    # 창 밖으로 밀려나 읽히지 못한 값 — 좌표가 아니므로 색 없는 점선으로만 표시.
    lost = bbox(base[0] - _ROW34_OD_SHIFT_PX + col_c[0], y3, 15, 17)
    sub = gs[0, col].subgridspec(2, 1, height_ratios=[1.0, 0.34], hspace=0.10)
    # 좌측 파란 표식과 우측의 잘린 숫자가 들어오지 않도록 여백을 좁게 잡는다.
    vf_view = (lost[0] - 5, y3 - 28, read_b[3][2] + 5, y3 + 28)
    show(fig.add_subplot(sub[0]), im_f, vf_view,
         [(b, C_FAIL) for b in read_b], ovals=[lost])
    strip(fig.add_subplot(sub[1]), im_f, read_b, under=VF_TRUE)
    title(col, 'VF 24-2 grid, row 3')
    caption(col, ', '.join(VF_READ) + ', …', ', '.join(VF_TRUE) + ', …')

    # ── 색상 키 (텍스트 색 대신 그림 안 사각형으로) ──────────────────
    # 초록(올바른 좌표)은 RNFL 패널에만 쓰이므로, 그 패널을 빼면 키에서도 뺀다.
    items = [(C_FAIL, 'rect', 'superseded coordinates')]
    if INCLUDE_PANEL_D:
        items.append((C_FIX, 'rect', 'corrected coordinates'))
    items.append((C_NEUTRAL, 'oval', 'value not read'))
    kw, kh = 0.014, 0.020
    step = 0.245
    kx = 0.5 - (len(items) * step - (step - 0.16)) / 2
    ky = 0.105
    for color, kind, lab in items:
        if kind == 'rect':
            fig.patches.append(Rectangle((kx, ky), kw, kh, transform=fig.transFigure,
                                         fill=False, edgecolor=color, linewidth=LW))
        else:
            fig.patches.append(Ellipse((kx + kw / 2, ky + kh / 2), kw * 1.2, kh,
                                       transform=fig.transFigure, fill=False,
                                       edgecolor=color, linewidth=1.0,
                                       linestyle=(0, (1.6, 1.6))))
        fig.text(kx + kw + 0.010, ky + kh / 2, lab, fontsize=FS_NOTE, va='center',
                 color='black')
        kx += step

    fig.text(0.5, 0.062,
             'Top: where the coordinates pointed.  Middle: the pixels read there, with '
             f'the reference value beneath each in {letters[-1]}.\n'
             'Coordinates are those used by the released code; report regions shown '
             'contain no identifiers.',
             fontsize=FS_NOTE, color='black', ha='center', va='top', linespacing=1.6)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    # bbox_inches='tight' 를 쓰면 캔버스가 내용에 맞춰 잘려 폭이 806pt 에서 벗어난다.
    # 결합 시 배율이 어긋나 글자 크기가 안 맞으므로 캔버스를 그대로 저장한다.
    for ext in ('pdf', 'png'):
        fig.savefig(OUT_DIR / f'fig_failure_cases.{ext}', dpi=300)
    print(f'저장: {OUT_DIR / "fig_failure_cases.pdf"}')
    print(f'  패널 {n_panels}개 (INCLUDE_PANEL_D={INCLUDE_PANEL_D}), 라벨 {letters}')
    print(f'  GCIPL 창 {win} — extracted 8 vs reference 69')
    print(f'  VF 읽힌 열 {[b[0] for b in read_b]}, 밀려난 값 {lost[0]} (shift {_ROW34_OD_SHIFT_PX}px)')


if __name__ == '__main__':
    build()
