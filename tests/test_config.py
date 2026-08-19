"""config 골격 스모크 테스트."""
from hvf.config import get, load_config, seed, set_seed


def test_config_loads():
    cfg = load_config()
    assert isinstance(cfg, dict)
    assert "constants" in cfg and "paths" in cfg


def test_seed_is_int():
    assert isinstance(seed(), int)
    assert set_seed() == seed()


def test_domain_constants_present():
    assert get("constants", "vf", "max_db") == 40.0
    assert get("constants", "vf", "n_points") == 54
    # 2026-07-23: +-180d 창은 원고에서 삭제됐고 config 에서도 제거했다.
    # 코드가 논문에 없는 창을 약속하면 안 되므로 [90] 이 정답이다.
    assert get("constants", "matching", "windows_days") == [90]
    # 미확정 상수는 명시적으로 null 이어야 함(임의 기본값 금지)
    assert get("constants", "oct", "signal_strength_min") is None
