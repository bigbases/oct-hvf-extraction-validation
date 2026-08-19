#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""공개(push) 직전 PHI 무결성 검증 — 세 층을 각각 따로 판정한다.

이 스크립트가 존재하는 이유: 2026-08-17~19 감사에서 검사 자체가 다섯 번 틀렸다.
(1) 작업 트리만 읽어 HEAD 오염을 놓쳤고, (2) `\\b1[0-9]{7}\\b` 가 밑줄에 붙은
환자ID(`_19999999_`)를 놓쳤고, (3) Windows 에서 subprocess 가 stdin 개행을
`\\r\\n` 으로 바꿔 `git check-ignore --stdin` 이 전 파일을 오판했다.
그래서 판정 전에 **양성 대조**로 패턴과 gitignore 판정이 실제 작동하는지
증명하고, "깨끗하다"고 쓸 때 어느 층이 깨끗한지 반드시 명시한다.

  L0  파일명    — 이름에 PHI 가 박힌 파일(바이너리도 여기서 걸린다)
  L1  작업 트리 — gitignore 제외분을 뺀, 새 저장소로 복사될 파일 집합
  L1b 컨테이너  — xlsx/docx(zip+XML), pdf(압축 스트림)를 풀어서 검사
  L2  HEAD      — 지금 커밋돼 있는 내용 (`git show HEAD:<path>`)
  L3  이력 전체 — 모든 커밋의 모든 blob (`git rev-list --all --objects`)

L3 는 이 저장소에서 **오염 상태로 남아 있는 것이 정상**이다(과거 커밋을
재작성하지 않기로 했으므로). 공개는 L1 만 새 저장소로 옮기는 방식이며,
따라서 push 게이트는 L1·L2 이고 L3 는 "이 .git 을 올리면 안 된다"의 근거다.

종료 코드: L0·L1·L1b·L2 중 하나라도 적출되면 1, 아니면 0.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# ── 패턴 ────────────────────────────────────────────────────────────────
# 경계에 \b 를 쓰지 않는다 — 밑줄이 단어문자라 `_19999999_` 를 놓친다.
PATTERNS = {
    # 이 코호트의 식별자는 5~8자리이고 1 로 시작하지 않는 것도 있다.
    # 초판이 8자리·1 시작만 보고 있어 6자리 식별자를 공개본에 그대로
    # 내보냈다(예시 값을 여기 적으면 이 파일이 다시 적출 대상이 되므로
    # 적지 않는다). 'patient'/'pid'/'환자' 문맥을 요구해 SHA-256 조각·
    # 바이트수 같은 오탐을 피한다.
    '환자ID(문맥)': re.compile(
        b'(?:patient|pid|' + '환자'.encode('utf-8') + b')'
        rb'[\s:=_#]{0,3}(?<![0-9])\d{5,8}(?![0-9])', re.I),
    '환자ID(8자리)': re.compile(rb'(?<![0-9])1[0-9]{7}(?![0-9])'),
    'DOB/검사일': re.compile(rb'(?<![0-9])(?:19[0-9]{2}|20[0-2][0-9])'
                             rb'(?:0[1-9]|1[0-2])(?:0[1-9]|[12][0-9]|3[01])(?![0-9])'),
    # 하이픈/슬래시 표기(2022-12-19)도 검사일이다. 다만 작업 일자 주석이 많아
    # 'patient'/'환자'/'exam' 문맥이 있을 때만 적출한다.
    '검사일(문맥)': re.compile(
        b'(?:patient|exam|' + '환자'.encode('utf-8') + b'|' + '검사'.encode('utf-8') + b')'
        rb'[^\n]{0,30}?(?:19|20)\d{2}[-/](?:0[1-9]|1[0-2])[-/](?:0[1-9]|[12]\d|3[01])', re.I),
    '개인 경로': re.compile(rb'C:[\\/]{1,2}Users[\\/]{1,2}\w+'),
    '실명(폴더명 규칙)': re.compile(rb'(?<![0-9])1[0-9]{7}__([a-z][a-z_ -]{2,30}?)__', re.I),
}
# 확인된 오탐만 개별 문자열로 면제한다(패턴을 느슨하게 만들지 않는다).
ALLOW = (
    b'5.5.0.20241111',   # Tesseract 버전
    b'PMID 18998889',    # 참고문헌 PMID
    # 이 파일 자신의 양성 대조 프로브(합성값). 실제 환자ID·실명이 아니다.
    # 문자열 단위로만 면제하므로, 여기에 진짜 값이 들어오면 그대로 적출된다.
    b'19999999', b'test_person', b'19010101', b'20200101',
    # 문맥형 패턴용 프로브. 숫자만 면제하면 실제 값까지 가려지므로 문구째 등재.
    b'patient 999999', b'exam 2000-01-01',
    # 소스에는 이스케이프된 형태(역슬래시 2개)로 적혀 있으므로 두 표기를 모두 등재
    rb'C:\Users\nobody', rb'C:\\Users\\nobody',
    # Office 문서 내부값(환자ID 아님). theme1.xml 의 그라디언트 각도(270°를
    # 60000분의 1도 단위로 표기), settings.xml 의 rsid. 맥락 확인 완료.
    b'16200000', b'13174835',
)
# ── 보조 도구 흔적 ──────────────────────────────────────────────────────
# PHI 는 아니지만 공개 저장소에 남기지 않는 것: 집필·코딩 보조 도구의 이름,
# 커밋 트레일러, 그리고 공개본에 없는 내부 문서를 이름으로 지목하는 참조.
# (requirements-train.txt 가 CURSOR_BRIEFING 을 지목하고 있던 것을 이렇게 놓쳤다.)
TOOL_PATTERNS = {
    '보조도구 이름': re.compile(
        rb'claude|anthropic|chatgpt|gpt-4|openai|copilot|cursor|windsurf|'
        rb'codeium|tabnine|github\s*copilot', re.I),
    '커밋 트레일러': re.compile(rb'co-authored-by|generated\s+with\s+\[?', re.I),
    '내부 문서 참조': re.compile(rb'CURSOR_BRIEFING|CLAUDE\.md|NEW_VF_BATCH_HANDOFF'),
}
# ── 산문 검사 ──────────────────────────────────────────────────────────
# 패턴 매칭으로는 "문서에 뭐라고 써 있는지"를 못 잡는다. ENVIRONMENT.md 의
# 내부 서버 사양과 REPRODUCIBILITY_AUDIT.md 의 재현 실패 기록이 이렇게 공개
# 저장소에 그대로 나갔다. 아래는 그 재발을 막기 위한 최소 신호다 — 걸리면
# 자동으로 지우지 말고 사람이 그 문서를 읽어야 한다는 뜻이다.
PROSE_PATTERNS = {
    '미해결·TODO': ('TODO', 'FIXME', 'XXX:', 'STUB', '미해결', '미구현', '확보 필요',
                    '확정 필요', '이관 예정', '진행 중'),
    '재현 실패 서술': ('재현안됨', '부분재현', '재현 실패', '드리프트', '워킹트리',
                     '삭제됨', '불일치 확인', '격차'),
    '내부 인프라': ('conda', 'Ubuntu', '서버', 'GPU 서버', 'C:\\Program Files',
                   'C:/Program Files'),
}
# ── 한글 검사 ──────────────────────────────────────────────────────────
# 범위를 유니코드로 명시한다 — `[가-힣]` 을 셸 grep 에 그대로 쓰면 로케일에
# 따라 em dash(—)·ellipsis(…) 같은 다바이트 문자를 오탐한다(실제로 그렇게
# 잘못 센 적이 있다). 아래 세 블록만 한글이다:
#   U+AC00–U+D7A3 완성형 음절, U+1100–U+11FF 자모, U+3130–U+318F 호환 자모
HANGUL = re.compile('[가-힣ᄀ-ᇿ㄰-㆏]')
#
# **적출 대상은 아래 확장자뿐이다.** 독자가 저장소를 쓰기 위해 반드시 읽는
# 산문·인용 메타데이터가 여기 해당한다. 소스 주석의 한글은 의도적으로 남기기로
# 했고(README 에 그 사실을 밝혀 둠), 레지스트리 JSON 의 한글은 스크립트가 써 넣은
# 봉인 기록이라 손대면 재생성물과 어긋난다.
HANGUL_TARGETS = {'.md', '.txt', '.cff', '.rst'}
# 대상 밖 파일의 한글은 적출하지 않되 **줄 수를 참고로 출력한다** — 예전에
# ".md 한글 0줄" 을 "저장소 전체에 한글 없음" 으로 오해하게 만든 적이 있다.
HANGUL_INFO_ONLY = {'.py', '.json', '.yaml', '.yml', '.toml', '.cfg', '.ini', ''}

# 문서 성격상 정당한 표현은 등재해 면제한다.
PROSE_ALLOW = (
    'STUB. build_ml_final',      # stage_10_dataset.py 가 스텁임을 밝히는 주석
)

# 확인된 오탐만 개별 등재한다.
TOOL_ALLOW = (
    b'Jean-Claude',                 # 참고문헌 저자명 (Mwanza, Jean-Claude)
    b'generated with the docstrip', # LaTeX 패키지 상용구
)

# 파일명 오탐: 환자 검사일이 아니라 **작업 일자**가 이름에 박힌 경우.
# 하나씩 눈으로 확인하고 등재한다 — 패턴을 넓히지 않는다.
FILENAME_ALLOW = (
    'kcc/session_log_20260427.md',
    'results/_backup_phase5_1_275_20260803.json',
    'results/_backup_phase5_3_275_20260803.json',
    'results/_backup_phase5_4_275_20260803.json',
    'results/_backup_sita_276_20260803.json',
    'results/_backup_sita_sensitivity_276_20260803.json',
)
# FILENAME_ALLOW 에 적힌 파일명이 이 소스 안에 리터럴로 들어 있어, 내용 검사에서
# 이 파일 자신이 적출된다. 등재된 파일명 그대로만 면제한다(자기유지형).
ALLOW = ALLOW + tuple(f.encode('utf-8') for f in FILENAME_ALLOW)

TEXT_EXT = {'.py', '.md', '.txt', '.yaml', '.yml', '.json', '.sh', '.tex', '.cfg',
            '.ini', '.toml', '.cff', '.csv', '.bat', '.bib', '.log', ''}
SKIP_DIRS = {'.git', '__pycache__', 'node_modules', '.venv', 'venv', '.pytest_cache'}


def scan(body: bytes, apply_allow: bool = True):
    """면제 문자열 제거 후 패턴별 적출 수를 센다.

    apply_allow=False 는 양성 대조 전용이다 — 프로브 값이 ALLOW 에 등재돼
    있으므로, 면제를 그대로 적용하면 대조가 항상 통과해버려 무의미해진다.
    """
    if apply_allow:
        for a in ALLOW:
            body = body.replace(a, b'')
    out = {}
    for name, rx in PATTERNS.items():
        hits = rx.findall(body)
        if hits:
            out[name] = len(set(hits))
    return out


def git(*args, stdin: bytes = None) -> bytes:
    """git 호출. core.quotepath=false 가 필수다 — 기본값이면 비ASCII 경로를
    8진 이스케이프로 인용해 돌려주므로(`"...\354\265\234..."`), 한글 파일명이
    gitignore 판정에서 통째로 누락된다(`sfa_review_master_최종.xlsx` 를 이렇게
    놓쳤다). ls-files·rev-list 출력에도 같은 문제가 있다."""
    return subprocess.run(['git', '-c', 'core.quotepath=false', *args],
                          cwd=ROOT, input=stdin, capture_output=True).stdout


# ── 0. 양성 대조 ────────────────────────────────────────────────────────
def positive_control() -> bool:
    """패턴이 실제로 작동하는지, gitignore 판정이 맞는지 먼저 증명한다."""
    print('=' * 74)
    print(' 0. 양성 대조 — 검사기 자체가 작동하는지 먼저 증명')
    print('=' * 74)
    ok = True

    # 프로브에는 반드시 가짜 값을 쓴다 — 실제 환자ID·실명을 넣으면 검사기가
    # 스스로 PHI 보유 파일이 된다(초판에서 실제로 이 실수를 했다).
    # 프로브는 패턴 종류마다 하나씩 대응하는 합성 문자열이어야 한다.
    # 문맥형 패턴(5~8자리 ID, 하이픈 날짜)을 빠뜨리면 대조가 통과해도
    # 그 패턴이 죽었는지 알 수 없다.
    probe = ('cirrus_out/by_case/19999999__test_person__19999999_19010101_'
             '20200101_OPT/x.png C:\\Users\\nobody\\OneDrive '
             'patient 999999 exam 2000-01-01').encode('utf-8')
    got = scan(probe, apply_allow=False)
    for name in PATTERNS:
        hit = name in got
        print('  %-22s %s' % (name, '검출 ✓' if hit else '검출 실패 ✗'))
        ok &= hit

    # 밑줄에 붙은 ID 를 놓치지 않는지 (과거 실패 재발 방지)
    underscore_ok = '환자ID(8자리)' in scan(b'_19999999_', apply_allow=False)
    print('  %-22s %s' % ('밑줄 인접 ID', '검출 ✓' if underscore_ok else '검출 실패 ✗'))
    ok &= underscore_ok

    # 음성 대조 — SHA-256 조각·작업일자를 환자ID/검사일로 오인하지 않는지
    fp = scan(b'sha256 aa4a75bba04464120e2b5a9 / 2026-07-23 patch', apply_allow=False)
    no_fp = not fp
    print('  %-22s %s' % ('SHA·작업일자 오탐', '없음 ✓' if no_fp else '오탐 %s ✗' % list(fp)))
    ok &= no_fp

    # 음성 대조 — 면제 문자열이 실제로 면제되는지
    clean = not scan(b'Tesseract 5.5.0.20241111 / PMID 18998889')
    print('  %-20s %s' % ('오탐 면제', '정상 ✓' if clean else '면제 실패 ✗'))
    ok &= clean

    print('\n  => 검사기 %s\n' % ('정상 — 아래 판정을 신뢰할 수 있다' if ok
                                  else '고장 — 아래 판정은 무의미하다'))
    return ok


def ignored_set(rels):
    """gitignore 판정. 입력을 bytes 로 넘긴다 — text=True 면 Windows 에서
    개행이 \\r\\n 으로 바뀌어 git 이 \\r 를 파일명 일부로 받고 전부 오판한다."""
    out = git('check-ignore', '--stdin', stdin='\n'.join(rels).encode('utf-8'))
    got = {x.decode('utf-8') for x in out.split(b'\n') if x.strip()}
    # 교차검증: 표본 몇 건을 단건 호출로 다시 확인
    sample = list(got)[:3] + [r for r in rels if r not in got][:3]
    for s in sample:
        single = subprocess.run(['git', 'check-ignore', '-q', s],
                                cwd=ROOT).returncode == 0
        if single != (s in got):
            print('  ⚠ gitignore 판정 불일치: %s (일괄=%s, 단건=%s)'
                  % (s, s in got, single))
    return got


def report(title, rows, note=''):
    print('=' * 74)
    print(' %s' % title)
    print('=' * 74)
    if note:
        print(' %s' % note)
    if not rows:
        print('  적출 0건 ✓\n')
        return 0
    for rel, hits in sorted(rows):
        print('  %-50s %s' % (rel[:50], ', '.join('%s %d' % kv for kv in hits.items())))
    print('  => 적출 %d개 파일\n' % len(rows))
    return len(rows)


def main():
    if not positive_control():
        print('양성 대조 실패 — 검증을 중단한다.')
        return 2

    files = [p for p in ROOT.rglob('*')
             if p.is_file() and not (SKIP_DIRS & set(p.parts))]
    rels = [p.relative_to(ROOT).as_posix() for p in files]
    ign = ignored_set(rels)
    tracked = {x for x in git('ls-files').decode('utf-8').split('\n') if x.strip()}

    # ── L0 파일명 ───────────────────────────────────────────────────────
    # 내용이 아니라 이름에 PHI 가 박힌 경우. 초판이 놓쳤던 층이다 —
    # `_crop_gcl_19999999.png` 처럼 바이너리라 내용 검사도 안 되는 파일이
    # 이름만으로 환자를 특정한다.
    l0 = [(rel, h) for rel in rels if rel not in ign and rel not in FILENAME_ALLOW
          for h in [scan(rel.encode('utf-8'))] if h]
    n0 = report('L0. 파일명 (gitignore 제외분 제외)', l0,
                '이름만으로 환자를 특정할 수 있는 파일. 바이너리도 여기서 걸린다.')

    # ── L1 작업 트리 (새 저장소로 복사될 집합) ──────────────────────────
    l1, bins = [], []
    for p, rel in zip(files, rels):
        if rel in ign:
            continue
        if p.suffix.lower() not in TEXT_EXT:
            bins.append(rel)
            continue
        h = scan(p.read_bytes())
        if h:
            l1.append((rel, h))
    n1 = report('L1. 작업 트리 내용 (gitignore 제외분 제외)',
                l1, '이 층이 곧 공개 대상이다. 여기가 깨끗해야 push 가능.')

    # ── L1c 보조 도구 흔적 ──────────────────────────────────────────────
    # 파일 내용 + 커밋 메시지 + 커밋 저자/커미터 메타데이터를 함께 본다.
    l1c = []
    for p, r in zip(files, rels):
        # 이 파일 자신은 제외한다 — TOOL_PATTERNS 리터럴이 정의상 여기 들어 있어
        # 항상 자기 자신을 적출한다. 대신 아래 자기검사로 패턴 작동을 증명한다.
        if r in ign or p.suffix.lower() not in TEXT_EXT or r == 'scripts/verify_no_phi.py':
            continue
        body = p.read_bytes()
        for a in TOOL_ALLOW:
            body = body.replace(a, b'')
        h = {n: len(set(rx.findall(body))) for n, rx in TOOL_PATTERNS.items()
             if rx.search(body)}
        if h:
            l1c.append((r, h))
    meta = git('log', '--format=%B%n%an%n%ae%n%cn%n%ce')
    for a in TOOL_ALLOW:
        meta = meta.replace(a, b'')
    mh = {n: len(set(rx.findall(meta))) for n, rx in TOOL_PATTERNS.items()
          if rx.search(meta)}
    if mh:
        l1c.append(('(커밋 메시지·저자 메타데이터)', mh))
    # 자기검사: 이 파일을 제외했으므로, 패턴이 죽지 않았음을 합성 문자열로 증명
    probe_ok = all(rx.search(b'co-authored-by claude cursor CURSOR_BRIEFING')
                   for rx in TOOL_PATTERNS.values())
    n1c = report('L1c. 보조 도구 흔적 (파일 내용 + 커밋 메시지·저자)', l1c,
                 'PHI 는 아니지만 공개 저장소에 남기지 않는다. '
                 '패턴 자기검사: %s' % ('정상 ✓' if probe_ok else '고장 ✗'))

    # ── L1d 산문 검사 (문서·주석의 내부 상태 서술) ──────────────────────
    # 대상은 사람이 읽는 문서다(.md). 코드 주석까지 걸면 잡음이 너무 커진다.
    l1d = []
    for p, r in zip(files, rels):
        if r in ign or p.suffix.lower() != '.md' or r == 'scripts/verify_no_phi.py':
            continue
        txt = p.read_text(encoding='utf-8', errors='ignore')
        for a in PROSE_ALLOW:
            txt = txt.replace(a, '')
        h = {}
        for name, words in PROSE_PATTERNS.items():
            found = sorted({w for w in words if w.lower() in txt.lower()})
            if found:
                h[name] = len(found)
        if h:
            l1d.append((r, h))
    n1d = report('L1d. 산문 검사 (.md 의 내부 상태 서술)', l1d,
                 '자동 판정이 아니라 "이 문서를 사람이 다시 읽으라"는 신호다.')

    # ── L1e 한글 검사 (공개 산문·설정) ──────────────────────────────────
    l1e = []
    for p, r in zip(files, rels):
        if r in ign or p.suffix.lower() not in HANGUL_TARGETS:
            continue
        if r == 'scripts/verify_no_phi.py':
            continue
        lines = p.read_text(encoding='utf-8', errors='ignore').split('\n')
        n = sum(1 for ln in lines if HANGUL.search(ln))
        if n:
            l1e.append((r, {'한글 포함 줄': n}))
    n1e = report('L1e. 한글 검사 — 적출 대상: %s' % ' '.join(sorted(HANGUL_TARGETS)),
                 l1e,
                 '독자가 저장소를 쓰기 위해 읽는 산문·인용 메타데이터만 적출한다.\n'
                 ' 소스 주석과 레지스트리 JSON 은 대상이 아니다(아래 참고 수치).')
    # 대상 밖 한글은 참고로만 보고한다 — 이 층의 0건을 "저장소 전체 0건"으로
    # 오해하지 않도록, 어디에 얼마나 남아 있는지 항상 함께 보여준다.
    info = {}
    for p, r in zip(files, rels):
        if r in ign or p.suffix.lower() in HANGUL_TARGETS:
            continue
        if p.suffix.lower() not in HANGUL_INFO_ONLY:
            continue
        n = sum(1 for ln in p.read_text(encoding='utf-8', errors='ignore').split('\n')
                if HANGUL.search(ln))
        if n:
            key = p.suffix.lower() or '(확장자 없음)'
            f_cnt, l_cnt = info.get(key, (0, 0))
            info[key] = (f_cnt + 1, l_cnt + n)
    if info:
        print(' 참고 — 적출 대상 밖의 한글(의도적으로 남김, README 에 명시):')
        for k in sorted(info, key=lambda x: -info[x][1]):
            print('   %-16s %3d개 파일  %5d줄' % (k, info[k][0], info[k][1]))
        print('   합계 %d줄\n' % sum(v[1] for v in info.values()))

    # ── L2 HEAD ────────────────────────────────────────────────────────
    l2 = []
    for rel in sorted(tracked):
        if Path(rel).suffix.lower() not in TEXT_EXT:
            continue
        h = scan(git('show', 'HEAD:' + rel))
        if h:
            l2.append((rel, h))
    n2 = report('L2. HEAD — 지금 커밋돼 있는 내용', l2,
                '작업 트리가 아니라 커밋된 blob 을 읽는다.')

    # ── L3 이력 전체 ───────────────────────────────────────────────────
    pairs = []
    for line in git('rev-list', '--all', '--objects').decode('utf-8', 'replace').split('\n'):
        if ' ' in line.strip():
            sha, path = line.strip().split(' ', 1)
            pairs.append((sha, path))
    checks = git('cat-file', '--batch-check',
                 stdin='\n'.join(s for s, _ in pairs).encode()).decode().split('\n')
    isblob = {c.split()[0]: (c.split()[1] == 'blob') for c in checks if len(c.split()) >= 2}
    l3, seen = {}, set()
    for sha, path in pairs:
        if not isblob.get(sha) or (sha, path) in seen:
            continue
        seen.add((sha, path))
        if Path(path).suffix.lower() not in TEXT_EXT:
            continue
        h = scan(git('cat-file', 'blob', sha))
        if h:
            cur = l3.setdefault(path, {})
            for k, v in h.items():
                cur[k] = max(cur.get(k, 0), v)
    report('L3. 이력 전체 (%d 커밋)' % len(git('rev-list', '--all').decode().split()),
           list(l3.items()),
           '오염이 남아 있는 것이 예정된 상태다 — 이력은 재작성하지 않기로 했다.\n'
           ' 이 목록이 비어 있지 않은 한 이 .git 은 절대 원격에 올리지 않는다.')

    # ── L1b 바이너리 컨테이너 ───────────────────────────────────────────
    # xlsx/docx 는 zip+XML, pdf 는 압축 스트림이라 원시 정규식으로는 안 잡힌다.
    # 풀어서 검사한다. 이 층을 안 봐서 kcc/sfa_review_A.xlsx(환자ID 29개)를
    # 처음에 놓쳤다. 나머지 바이너리(png/jpg 등)는 여전히 육안 확인 대상이다.
    l1b, opaque = [], []
    for rel in bins:
        p = ROOT / rel
        low = rel.lower()
        try:
            if low.endswith(('.xlsx', '.docx', '.pptx')):
                import zipfile
                z = zipfile.ZipFile(p)
                body = b''.join(z.read(n) for n in z.namelist())
            elif low.endswith('.pdf'):
                import zlib
                raw = p.read_bytes()
                parts = []
                for m in re.finditer(rb'stream\r?\n(.*?)endstream', raw, re.S):
                    try:
                        parts.append(zlib.decompress(m.group(1)))
                    except Exception:
                        parts.append(m.group(1))
                # PDF 좌표 오탐 제거: 소수점 뒤 숫자열은 ID·날짜가 될 수 없다.
                # matplotlib 이 찍는 Td 피연산자(텍스트 위치 오프셋)의 소수부가
                # 8자리라 환자ID 패턴과 그대로 겹친다.
                body = re.sub(rb'\.\d+', b'', b''.join(parts))
            else:
                opaque.append(rel)
                continue
        except Exception:
            opaque.append(rel)
            continue
        h = scan(body)
        if h:
            l1b.append((rel, h))
    n1b = report('L1b. 바이너리 컨테이너 내용 (xlsx/docx/pdf 압축 해제 후)', l1b,
                 'zip+XML 과 PDF 스트림을 풀어서 본다. 원시 바이트로는 안 잡히는 층.')

    # ── 바이너리: 자동 판정 불가, 육안 확인 대상 ────────────────────────
    print('=' * 74)
    print(' 참고. 불투명 바이너리 %d개 — 자동 판정 불가(육안 확인 필요)' % len(opaque))
    print('=' * 74)
    print(' PDF/PNG 안의 문자는 압축돼 있어 정규식으로 못 잡는다.')
    print(' 그림에 환자 식별정보가 인쇄돼 있지 않은지는 사람이 봐야 한다.')
    for b in sorted(opaque)[:20]:
        print('   %s' % b)
    if len(opaque) > 20:
        print('   … 외 %d개' % (len(opaque) - 20))

    print('\n' + '=' * 74)
    remote = git('remote', '-v').decode().strip()
    print(' 원격: %s' % (remote if remote else '없음 (push 이력 없음)'))
    print(' 판정: L0 %d / L1 %d / L1b %d / L1c %d / L1d %d / L1e %d / L2 %d → %s'
          % (n0, n1, n1b, n1c, n1d, n1e, n2, 'PUSH 불가' if (n0 or n1 or n1b or n1c or n1d or n1e or n2) else '전부 깨끗'))
    if l3:
        print(' L3 에 적출이 있다 — 이 .git 은 원격에 올릴 수 없다.')
        print(' 공개는 이력을 옮기지 않는 새 저장소로만 한다.')
    else:
        print(' L3 도 깨끗하다 — 이 저장소는 이력째로 공개해도 된다.')
    print('=' * 74)
    return 1 if (n0 or n1 or n1b or n1c or n1d or n1e or n2) else 0


if __name__ == '__main__':
    sys.exit(main())
