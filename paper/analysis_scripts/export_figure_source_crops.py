"""
scripts/export_figure_source_crops.py
그림 재료 추출 — 다른 도구로 그림을 다시 만들 때 쓰는 원본 크롭 + 좌표.

fig_input_example_full.pdf / fig_failure_cases.pdf 가 쓰는 원본 크롭을 원해상도
그대로(크기 조정 없이) PNG 로 저장하고, 그 위에 그려야 할 크롭 창 좌표를
**각 PNG 의 픽셀 좌표계**로 함께 기록한다.

마스킹: 기기 UI 한글 라벨 마스크(make_failure_cases.KO_MASK_RNFL)를 적용한 뒤
자른다. 다만 GCIPL/VF 크롭 영역은 마스크 위치와 겹치지 않아 결과가 동일하다
(RNFL 리포트 전체를 보여주던 비활성 패널용 마스크였기 때문). 리포트 헤더의
식별정보도 크롭 범위 밖이라 잘려 나간다.

출력: paper/manuscript/figures/_source_crops/*.png + manifest.md
"""
import json
import sys
from pathlib import Path

from PIL import Image, ImageDraw

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'scripts'))
from make_failure_cases import (                                   # noqa: E402
    load_candidates, gcipl_geometry, vf_geometry, KO_MASK_RNFL, VF_READ, VF_TRUE)

OUT = _ROOT / 'paper' / 'manuscript' / 'figures' / '_source_crops'
SRC_PDF = _ROOT / 'paper' / 'manuscript' / 'figures' / 'fig_input_example.pdf'


def rel(box, view):
    """페이지 좌표 박스 → view 크롭 안의 (x, y, w, h)."""
    return (box[0] - view[0], box[1] - view[1], box[2] - box[0], box[3] - box[1])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []          # (파일명, 크기, 설명, [(라벨, x, y, w, h)])

    # ── (a)(b)(c): 원본 PDF 에 박혀 있는 크롭 3장 ────────────────
    import fitz
    doc = fitz.open(str(SRC_PDF))
    names = ['vf_threshold_grid', 'rnfl_clock_hour', 'gcipl_6sector']
    for i, info in enumerate(doc[0].get_images(full=True)):
        pix = fitz.Pixmap(doc, info[0])
        if pix.n > 4:
            pix = fitz.Pixmap(fitz.csRGB, pix)
        im = Image.frombytes('RGB', (pix.width, pix.height), pix.samples)
        im = im.transpose(Image.FLIP_TOP_BOTTOM)     # PDF 이미지 공간은 원점이 좌하단
        f = f'{i + 1:02d}_{names[i]}.png'
        im.save(OUT / f)
        rows.append((f, im.size, f'Figure 2 패널 ({chr(97 + i)}) — 그리는 창 없음', []))
    doc.close()

    cands = load_candidates()

    # ── (d) GCIPL 하부 섹터 ───────────────────────────────────────
    im_g, view_g, win_g = gcipl_geometry(cands)
    d = ImageDraw.Draw(im_g)
    for b in KO_MASK_RNFL:
        d.rectangle(b, fill=(255, 255, 255))
    chart = im_g.crop(view_g)
    chart.save(OUT / 'gcipl_sector_chart.png')
    rows.append(('gcipl_sector_chart.png', chart.size,
                 '(d) 상단 — GCIPL 6섹터 차트 (하부 섹터 주변)',
                 [('실패한 셀 창 (빨강)',) + rel(win_g, view_g)]))
    winc = im_g.crop(win_g)
    winc.save(OUT / 'gcipl_window_crop.png')
    rows.append(('gcipl_window_crop.png', winc.size,
                 '(d) 중단 — 그 창 안에 실제로 들어온 픽셀 (extracted 8 / reference 69)', []))

    # ── (e) VF 24-2 3행 ──────────────────────────────────────────
    im_v, view_v, read_b, lost = vf_geometry()
    row = im_v.crop(view_v)
    row.save(OUT / 'vf_row3.png')
    marks = [(f'읽은 셀 {i + 1} (빨강, 값 {VF_READ[i]})',) + rel(b, view_v)
             for i, b in enumerate(read_b)]
    marks.append(('못 읽은 값 (중립 점선, 참값 29)',) + rel(lost, view_v))
    rows.append(('vf_row3.png', row.size, '(e) 상단 — 24-2 격자 3행 좌측', marks))

    gap = 7
    crops = [im_v.crop(b) for b in read_b]
    w = sum(c.width for c in crops) + gap * (len(crops) - 1)
    h = max(c.height for c in crops)
    strip = Image.new('RGB', (w, h), (255, 255, 255))
    x, centers = 0, []
    for c in crops:
        strip.paste(c, (x, 0)); centers.append(x + c.width // 2); x += c.width + gap
    strip.save(OUT / 'vf_row3_window.png')
    rows.append(('vf_row3_window.png', strip.size,
                 f'(e) 중단 — 읽은 셀 4개를 {gap}px 간격으로 이어 붙임. '
                 f'각 조각 중심 x = {centers}. 아래에 참값 {VF_TRUE} 를 세로로 붙일 것.',
                 []))

    # ── manifest ─────────────────────────────────────────────────
    md = ['# 그림 재료 — 원본 크롭 + 크롭 창 좌표', '',
          '이미지는 전부 원해상도 그대로(크기 조정 없음). 좌표는 각 PNG 의 픽셀 기준 '
          '(x, y, w, h), 원점은 좌상단.', '']
    for f, size, desc, marks in rows:
        md.append(f'## `{f}`  ({size[0]} x {size[1]} px)')
        md.append(f'{desc}')
        if marks:
            md.append('')
            md.append('| 표시 | x | y | w | h |')
            md.append('|---|---|---|---|---|')
            for lab, x, y, mw, mh in marks:
                md.append(f'| {lab} | {x} | {y} | {mw} | {mh} |')
        md.append('')
    (OUT / 'manifest.md').write_text('\n'.join(md), encoding='utf-8')

    print(f'저장 위치: {OUT}')
    for f, size, _, marks in rows:
        print(f'  {f:28s} {size[0]:4d} x {size[1]:4d} px   창 {len(marks)}개')
    print(f'  manifest.md')


if __name__ == '__main__':
    main()
