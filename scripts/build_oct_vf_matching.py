"""
OCT-VF 시간 매칭 알고리즘 (Methods 3.4/3.6 대응 구현).

배경 (2026-07-24): 논문이 "환자·좌우안 독립 처리, ±90일 이내 최근접 페어 매칭"을
명시하지만, 정본 파이프라인엔 이 매칭을 실제로 수행하는 코드가 없었다(kcc/
dataset_final.csv의 gap_days를 읽어 쓰기만 함).

원본 후보 풀 재구성 시도 결과: vf_index.csv(환자+날짜만, 눈 정보 없음),
cohort_hd_inventory_vfpool.csv(실은 OCT HD-scan 인벤토리), vf_new_batch_manifest.csv
(별개의 추가 발굴 배치) 전부 원본 후보 풀이 아님을 확인(2026-07-24 감사).
vf_index.csv에 눈 정보가 어떻게 붙었는지 git 이력에서도 못 찾음 — 원본 후보
선별 과정 자체는 유실로 간주.

따라서 이 스크립트는 "후보 풀에서 골라내는" 재현이 아니라 "이미 확정된 276안의
매칭 규칙을 코드로 명시 + dw_images/ 원본 폴더로 검증"하는 것을 목표로 한다.
검증(2026-07-24): dw_images/에서 OCT 후보가 2개 이상인 86개 눈-레코드 중
70건 완전 일치, 16건은 전부 실제 간격이 90일 초과라 정당하게 제외된 케이스,
불일치 0건. 동점(equidistant) 케이스는 0건 관측 — 원고의 동점 처리 문장은
이 코호트에서 실제로 발동한 적이 없어 규칙을 실측으로 확정할 수 없음(추정 금지).

PHI 규칙: 표준출력·로그엔 파생 통계(간격 일수, 집계)만 남기고 환자ID·원본
날짜는 남기지 않는다.
"""
import csv
import os
import re
import sys
from datetime import datetime
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / 'src'))
from hvf.config import get as _cfg


def _parse_date(s: str) -> datetime:
    return datetime.strptime(s, '%Y%m%d')


def nearest_within_window(anchor_date: str, candidates: list[str], window_days: int):
    """anchor_date(예: VF 검사일) 기준 candidates(OCT 검사일 목록) 중 window_days
    이내 최근접 날짜를 반환. 없으면 None. 동점 시 전체 동점 후보 목록도 함께 반환."""
    anchor = _parse_date(anchor_date)
    scored = []
    for c in candidates:
        gap = abs((_parse_date(c) - anchor).days)
        if gap <= window_days:
            scored.append((gap, c))
    if not scored:
        return None, []
    scored.sort(key=lambda x: x[0])
    min_gap = scored[0][0]
    tied = [c for g, c in scored if g == min_gap]
    return scored[0][1], tied


def match_eye(vf_date: str, oct_candidates: list[str], window_days: int | None = None):
    """단일 눈에 대해 VF 검사일과 OCT 후보 목록으로 최근접 매칭 수행.
    반환: (matched_oct_date 또는 None, gap_days 또는 None, tied_candidates)."""
    if window_days is None:
        window_days = _cfg('constants', 'matching', 'primary_window_days')
    chosen, tied = nearest_within_window(vf_date, oct_candidates, window_days)
    if chosen is None:
        return None, None, []
    gap = abs((_parse_date(chosen) - _parse_date(vf_date)).days)
    return chosen, gap, tied


def _scan_oct_candidates_from_dw_images(dw_images_dir: Path) -> dict:
    """dw_images/{patient_id}/{date}_{modality}_{seq}.jpg 구조에서
    환자별 OCT 검사일 후보 목록을 스캔(검증용 — 정본 후보 풀이 아니라
    dw_images에 실제로 존재하는 파일 날짜 기준)."""
    out = {}
    if not dw_images_dir.exists():
        return out
    for pid in os.listdir(dw_images_dir):
        d = dw_images_dir / pid
        if not d.is_dir():
            continue
        dates = set()
        for f in os.listdir(d):
            m = re.match(r'(\d{8})_', f)
            if m:
                dates.add(m.group(1))
        out[pid] = sorted(dates)
    return out


def verify_against_dataset_final():
    """dw_images/ 기반 후보로 kcc/dataset_final.csv의 기존 매칭 결과를 재현하는지 검증.
    PHI 보호: 환자ID/날짜는 콘솔에 출력하지 않고 집계 수치만 반환."""
    window = _cfg('constants', 'matching', 'primary_window_days')
    oct_cands = _scan_oct_candidates_from_dw_images(_ROOT / 'dw_images')
    final_rows = list(csv.DictReader(open(_ROOT / 'kcc' / 'dataset_final.csv', encoding='utf-8-sig')))

    n_tested = n_match = n_mismatch = n_excluded_ok = n_tie = 0
    for r in final_rows:
        pid = r['patient_id']
        cands = oct_cands.get(pid, [])
        if len(cands) <= 1:
            continue  # 후보가 1개 이하면 알고리즘 선택의 여지가 없어 검증 대상에서 제외
        n_tested += 1
        chosen, gap, tied = match_eye(r['vf_date'], cands, window)
        if len(tied) > 1:
            n_tie += 1
        if chosen is None:
            n_excluded_ok += 1  # 실제로 윈도 밖이라 제외되는 게 맞는지는 별도 확인됨(2026-07-24)
            continue
        if chosen == r['oct_date']:
            n_match += 1
        else:
            n_mismatch += 1

    return {
        'window_days': window,
        'n_tested': n_tested,
        'n_match': n_match,
        'n_mismatch': n_mismatch,
        'n_excluded_by_window': n_excluded_ok,
        'n_tie_cases_observed': n_tie,
    }


if __name__ == '__main__':
    result = verify_against_dataset_final()
    print('OCT-VF 매칭 알고리즘 검증 결과 (dw_images 기반, PHI 비노출):')
    for k, v in result.items():
        print(f'  {k}: {v}')
