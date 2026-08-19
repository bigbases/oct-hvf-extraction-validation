"""STAGE 00 — raw(PHI) 이미지 → 추출 CSV (OCR).

STUB. 다음 단계에서 기존 OCR 코드(ocr_threshold/ocr_oct_values/ocr_rnfl_detail 등)를
src/hvf/ 로 이관해 여기서 호출한다. 상수·경로·seed 는 hvf.config 에서만 읽는다.
"""
from hvf.config import set_seed, configure_tesseract

STAGE = "stage_00_extract"


def main() -> None:
    set_seed()
    configure_tesseract()
    print(f"[{STAGE}] STUB — OCR 추출 이관 대기. (raw → data/)")


if __name__ == "__main__":
    main()
