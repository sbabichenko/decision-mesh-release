# Release validation - 2026-07-29

The archive contents were validated from a clean extracted-style tree with GCC 14.2.0 and CMake in Release mode.

- `scripts/build_release.sh`: build succeeded.
- CTest: `lindsey_guard` passed (1/1).
- `scripts/verify_reference.sh`: direct-BH and adaptive-Lindsey mesh/hierarchy SHA-256 checks both passed.
- Updated paper: 65-page PDF compiled successfully and passed render inspection.

Frozen full-pipeline wall times on the release machine were 7.43 seconds for direct candidate BH and 6.51 seconds for adaptive Lindsey. Stage-adaptation timestamps were approximately 4.9 and 4.5 seconds, respectively. See `PERFORMANCE_UPDATE_2026-07-29.md` for the precise timing scope and configuration.
