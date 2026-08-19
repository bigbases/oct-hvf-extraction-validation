# `results/` — the value registry

Every number that appears in the manuscript exists here as JSON. No figure is
transcribed by hand.

- One file per value: `results/<key>.json`
- Written through `hvf.registry.record_result(key, value, script=..., inputs=[...], extra={...})`
- Each record holds the value, the script that computed it, the SHA-256 of every
  input file, the UTC timestamp, the git commit, and the Python version
- Because the input hashes are recorded, a change in the input data is visible
  as a mismatch rather than passing silently

These files hold aggregate values only — no per-eye rows — and are therefore
part of the published repository.
