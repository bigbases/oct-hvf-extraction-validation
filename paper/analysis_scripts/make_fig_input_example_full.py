"""
scripts/make_fig_input_example_full.py
Figure 2 결합본: 입력 예시 (a)(b)(c) + 추출 실패 사례 (d)(e).

Display item 1안(2026-08-16 확정)으로 실패 사례를 독립 Figure 로 두지 않고 기존
입력 예시 Figure 에 흡수한다.

왜 PDF 를 이어 붙이지 않고 다시 그리나:
  원본 fig_input_example.pdf 는 806pt 폭 캔버스에 9pt 글자로 만들어져 있어, 본문
  폭(17.6cm = 499pt)에 넣으면 0.62배로 축소돼 글자가 지면에서 5.6pt 밖에 안 된다.
  두 PDF 를 이어 붙이는 방식으로는 이걸 못 고친다 — 블록이 폭을 채우는 한 배율이
  상쇄돼 언제나 5.6pt 가 되기 때문. 그래서 원본 PDF 에 박혀 있는 세 장의 원본 크롭
  이미지를 그대로 꺼내 와서, 캔버스를 본문 폭과 1:1(17.6cm)로 두고 전부 다시 그린다.
  배율이 1 이므로 지정한 8pt 가 지면에서 그대로 8pt 다.

  부수 효과로 Figure 2 가 스크립트 산출물이 된다(원본은 생성 코드가 없는 수작업물
  이었다). 원본 PDF 는 이미지 출처로만 쓰고 수정하지 않는다.

(d)(e) 패널은 make_failure_cases.py 의 draw_panel() 을 그대로 쓴다 — 좌표 정의가
한 곳에만 있어야 두 그림이 어긋나지 않는다.

출력: paper/manuscript/figures/fig_input_example_full.{pdf,png}
"""
import sys
from pathlib import Path

import fitz
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
import numpy as np
from PIL import Image

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'scripts'))
from make_failure_cases import (                                  # noqa: E402
    load_candidates, draw_panel, PANEL_ASPECT,
    C_FAIL, C_NEUTRAL, FS_TITLE, FS_LAB, FS_NOTE, LW, PAGE_W_IN, _bare)

FIGS = _ROOT / 'paper' / 'manuscript' / 'figures'
SRC_PDF = FIGS / 'fig_input_example.pdf'      # (a)(b)(c) 원본 크롭의 출처
OUT = FIGS / 'fig_input_example_full'

# 원본에 박혀 있는 순서대로. 라벨의 픽셀 크기 표기는 뺀다 — 다시 그리면서 배율이
# 1:1 이 아니게 되어 "native" 가 사실이 아니게 되기 때문.
TOP_TITLES = ['VF threshold grid', 'RNFL clock-hour', 'GCIPL 6-sector']

# 세로 배치(인치). 각 항목의 높이를 명시해 아래쪽에 빈 공간이 남지 않게 한다.
H_LABEL = 0.16      # 패널 제목 줄
H_GAP_ROW = 0.30    # (c) 아래 ~ (d) 제목 사이
H_TIER2 = 0.42      # 2단(읽힌 픽셀) 띠
H_TIER_GAP = 0.04   # 1단 ~ 2단 사이
H_CAP_TOP = 0.09    # 2단 ~ extracted 라벨 사이 (붙여 둔다)
H_CAPTION = 0.34    # extracted / reference 두 줄
H_KEY = 0.22        # 색상 키
H_FOOT = 0.50       # 하단 설명 두 줄 — 두 줄이 온전히 들어갈 높이
MARGIN_X = 0.030    # 좌우 여백(figure 비율)


def extract_top_panels():
    """원본 PDF 에서 (a)(b)(c) 크롭 이미지를 꺼낸다(재저장하지 않음)."""
    if not SRC_PDF.exists():
        sys.exit(f'원본 PDF 없음: {SRC_PDF}')
    doc = fitz.open(str(SRC_PDF))
    page = doc[0]
    out = []
    for info in page.get_images(full=True):
        pix = fitz.Pixmap(doc, info[0])
        if pix.n > 4:
            pix = fitz.Pixmap(fitz.csRGB, pix)
        img = Image.frombytes('RGB', (pix.width, pix.height), pix.samples)
        # PDF 이미지 공간은 원점이 좌하단이라, 배치 행렬을 거치지 않고 꺼내면
        # 상하가 뒤집혀 있다. 되돌린다.
        out.append(img.transpose(Image.FLIP_TOP_BOTTOM))
    doc.close()
    if len(out) != 3:
        sys.exit(f'원본에서 이미지 3장을 기대했으나 {len(out)}장 발견')
    return out


def main():
    tops = extract_top_panels()
    cands = load_candidates()

    usable_w = PAGE_W_IN * (1 - 2 * MARGIN_X)
    # ── 상단 (a)(b)(c): 세 크롭을 같은 배율로 — 원본과 동일한 규칙 ──
    gap_a = 0.14
    widths_px = [im.width for im in tops]
    scale_a = (usable_w - gap_a * 2) / sum(widths_px)          # in/px
    h_top = max(im.height for im in tops) * scale_a
    # ── 하단 (d)(e): 열 너비를 1단 종횡비에 비례 → 두 패널 높이가 같아진다 ──
    asp = [PANEL_ASPECT['gcipl'], PANEL_ASPECT['vf']]
    gap_b = 0.22
    k = (usable_w - gap_b) / sum(asp)
    h_tier1 = k

    total = (H_LABEL + h_top + H_GAP_ROW + H_LABEL + h_tier1 + H_TIER_GAP + H_TIER2
             + H_CAP_TOP + H_CAPTION + H_KEY + H_FOOT)
    fig = plt.figure(figsize=(PAGE_W_IN, total))

    def fy(inch_from_top):
        return 1 - inch_from_top / total

    def fx(inch_from_left):
        return inch_from_left / PAGE_W_IN

    left0 = PAGE_W_IN * MARGIN_X

    # ── (a)(b)(c) ────────────────────────────────────────────────
    y = 0.0
    x = left0
    for i, im in enumerate(tops):
        w = im.width * scale_a
        h = im.height * scale_a
        fig.text(fx(x), fy(y + H_LABEL * 0.30), f'({chr(97 + i)}) {TOP_TITLES[i]}',
                 fontsize=FS_TITLE, ha='left', va='center')
        ax = fig.add_axes([fx(x), fy(y + H_LABEL + h), w / PAGE_W_IN, h / total])
        ax.imshow(np.asarray(im))
        _bare(ax)
        ax.add_patch(Rectangle((0, 0), im.width - 1, im.height - 1, fill=False,
                               edgecolor='0.6', linewidth=0.5, transform=ax.transData))
        x += w + gap_a
    y += H_LABEL + h_top + H_GAP_ROW

    # ── (d)(e) ───────────────────────────────────────────────────
    panel_h = h_tier1 + H_TIER_GAP + H_TIER2
    gs = fig.add_gridspec(1, 2, width_ratios=asp, wspace=gap_b / (usable_w / 2),
                          left=fx(left0), right=fx(left0 + usable_w),
                          top=fy(y + H_LABEL), bottom=fy(y + H_LABEL + panel_h))
    col_x = [left0, left0 + asp[0] * k + gap_b]
    for i, kind in enumerate(('gcipl', 'vf')):
        title, ex, ref, _ = draw_panel(fig, gs[0, i], kind, cands)
        fig.text(fx(col_x[i]), fy(y + H_LABEL * 0.30), f'({chr(100 + i)}) {title}',
                 fontsize=FS_TITLE, ha='left', va='center')
        cy = y + H_LABEL + panel_h + H_CAP_TOP + H_CAPTION * 0.26
        cw = asp[i] * k
        for lab, val in (('extracted:', ex), ('reference:', ref)):
            fig.text(fx(col_x[i] + cw * 0.44), fy(cy), lab, fontsize=FS_LAB,
                     fontweight='bold', ha='right', va='center')
            fig.text(fx(col_x[i] + cw * 0.47), fy(cy), val, fontsize=FS_LAB,
                     ha='left', va='center')
            cy += H_CAPTION * 0.42
    y += H_LABEL + panel_h + H_CAP_TOP + H_CAPTION

    # ── 색상 키 ──────────────────────────────────────────────────
    items = [(C_FAIL, 'rect', 'coordinates as used'),
             (C_NEUTRAL, 'oval', 'value not read')]
    kw, kh = 0.16, 0.10          # 인치
    step = 2.0
    kx = left0 + (usable_w - (step * (len(items) - 1) + 1.6)) / 2
    ky = y + H_KEY * 0.28
    for color, kind, lab in items:
        if kind == 'rect':
            fig.patches.append(Rectangle((fx(kx), fy(ky + kh / 2)), kw / PAGE_W_IN,
                                         kh / total, transform=fig.transFigure,
                                         fill=False, edgecolor=color, linewidth=LW))
        else:
            fig.patches.append(Ellipse((fx(kx + kw / 2), fy(ky)), kw * 1.15 / PAGE_W_IN,
                                       kh / total, transform=fig.transFigure, fill=False,
                                       edgecolor=color, linewidth=1.0,
                                       linestyle=(0, (1.6, 1.6))))
        fig.text(fx(kx + kw + 0.07), fy(ky), lab, fontsize=FS_NOTE, va='center')
        kx += step
    y += H_KEY

    fig.text(0.5, fy(y + H_FOOT * 0.42),
             '(d), (e) top: where the coordinates pointed;  middle: the pixels read '
             'there, with the reference value beneath each in (e).\n'
             'Coordinates are those used by the released code; report regions shown '
             'contain no identifiers.',
             fontsize=FS_NOTE, ha='center', va='center', linespacing=1.6)

    FIGS.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{OUT}.{ext}', dpi=300)
    print(f'저장: {OUT}.pdf')
    print(f'  캔버스 {PAGE_W_IN * 2.54:.1f} x {total * 2.54:.1f} cm  '
          f'(= 본문 폭 1:1, 글자 {FS_TITLE}pt 가 지면에서 그대로 {FS_TITLE}pt)')
    print(f'  상단 크롭 배율 {scale_a * 300:.2f}x @300dpi, 높이 {h_top * 2.54:.1f}cm')


if __name__ == '__main__':
    main()
