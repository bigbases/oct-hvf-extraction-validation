"""STAGE 20 — 데이터셋 → 논문 수치 (results registry).

STUB. KCC/저널 수치(카운트, OCR 누락률, 구조-기능 상관/회귀, Bland-Altman)를
계산해 hvf.registry.record_result 로 봉인. 하드코딩 숫자 금지.

이관 시 각 수치는 반드시:
  record_result(key, value, script=__file__, inputs=[...쓴 CSV들...], extra={포함/제외기준, n, 보정변수})
"""
from hvf.config import set_seed

STAGE = "stage_20_analysis"


def main() -> None:
    set_seed()
    print(f"[{STAGE}] STUB — 수치 계산+registry 봉인 이관 대기.")


if __name__ == "__main__":
    main()
