"""STAGE 10 — 추출 CSV → 통합 데이터셋 (ml_final_*).

STUB. build_ml_final / make_flip 로직 이관 예정. 매칭 윈도우·결측 규칙은 config 에서.
"""
from hvf.config import get, set_seed

STAGE = "stage_10_dataset"


def main() -> None:
    set_seed()
    windows = get("constants", "matching", "windows_days")
    print(f"[{STAGE}] STUB — 데이터셋 빌드 이관 대기. windows={windows}일")


if __name__ == "__main__":
    main()
