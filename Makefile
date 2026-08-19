# HVF 재현 파이프라인 — 원 데이터 → 최종 표/그림 한 번에.
#   make all      전체 파이프라인 (raw → results/)
#   make check    설정 검증 + 테스트 (코드 없이 골격만으로 동작)
#   make env      의존 설치
# 실행: Git Bash / WSL / Docker / 서버. (Windows 순정 cmd 에는 make 없음)
#
# ⚠️ 아래 stage_* 는 현재 STUB. 기존 코드 이관(다음 단계) 후 실제 동작.
#    DAG(의존 순서)만 먼저 고정해 둔다.

PY := python

.PHONY: all check env test config-check clean \
        extract dataset analysis figures paper

# ── 전체 재현 ───────────────────────────────────────────────
all: figures paper
	@echo "[make all] 완료. 산출물: results/ , paper/"

# raw(PHI) → 추출 CSV
extract:
	$(PY) scripts/stage_00_extract.py

# 추출 CSV → 통합 데이터셋 (ml_final_*)
dataset: extract
	$(PY) scripts/stage_10_dataset.py

# 데이터셋 → 논문 수치(registry JSON)
analysis: dataset
	$(PY) scripts/stage_20_analysis.py

# 데이터셋/registry → 그림
figures: analysis
	$(PY) scripts/stage_30_figures.py

# 그림/수치 → 논문 산출물
paper: figures
	@echo "[stub] paper/ 조립은 이관 후 구현"

# ── 검증 (지금 바로 동작) ────────────────────────────────────
check: config-check test

config-check:
	$(PY) -c "from hvf.config import load_config, seed; c=load_config(); print('config OK · seed=', seed(), '· windows=', c['constants']['matching']['windows_days'])"

test:
	pytest -q

env:
	pip install -r requirements-dev.txt
	pip install -e .

clean:
	rm -rf **/__pycache__ .pytest_cache
