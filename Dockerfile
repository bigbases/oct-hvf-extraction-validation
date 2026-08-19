# 분석·OCR 재현용 이미지 (학습/GPU 아님 — 그건 서버 별도).
# ⚠️ 결정론 주의: apt 의 tesseract 가 로컬(5.5.0.20241111)과 정확히 같지 않으면
#    OCR 포인트 추출이 미세하게 달라질 수 있다. OCR 재현이 목적이면 tesseract 5.5.0 을
#    소스빌드로 맞출 것. 여기서는 distro 제공본을 쓰고 버전을 로그로 남긴다.
FROM python:3.14-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=42 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        tesseract-ocr \
        git \
    && tesseract --version \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt requirements-dev.txt pyproject.toml ./
RUN pip install --no-cache-dir -r requirements-dev.txt

COPY . .
RUN pip install --no-cache-dir -e .

# 컨테이너 안의 tesseract 경로로 덮어쓰기 (config 의 Windows 경로 무시)
ENV HVF_TESSERACT_CMD=/usr/bin/tesseract

CMD ["make", "check"]
