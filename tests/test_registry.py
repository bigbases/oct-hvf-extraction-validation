"""results registry 라운드트립 + 입력 해시 봉인 테스트."""
import json

from hvf.registry import load_result, record_result, sha256_file


def test_record_and_load(tmp_path):
    # 입력 파일 하나 만들어 해시가 봉인되는지 확인
    inp = tmp_path / "input.csv"
    inp.write_text("a,b\n1,2\n", encoding="utf-8")

    out = record_result(
        "unit_test_metric",
        {"n": 3, "mean": 1.5},
        script="tests/test_registry.py",
        inputs=[inp],
        extra={"note": "smoke"},
    )
    assert out.exists()

    got = load_result("unit_test_metric")
    assert got["value"] == {"n": 3, "mean": 1.5}
    assert got["script"] == "tests/test_registry.py"
    assert got["inputs"][0]["sha256"] == sha256_file(inp)
    assert got["extra"]["note"] == "smoke"
    assert "created_utc" in got

    # JSON 으로 다시 읽히는지
    json.loads(out.read_text(encoding="utf-8"))

    out.unlink()  # 테스트 아티팩트 정리


def test_missing_input_flagged(tmp_path):
    out = record_result(
        "unit_test_missing",
        1,
        script="tests/test_registry.py",
        inputs=[tmp_path / "does_not_exist.csv"],
    )
    got = load_result("unit_test_missing")
    assert got["inputs"][0]["missing"] is True
    assert got["inputs"][0]["sha256"] is None
    out.unlink()
