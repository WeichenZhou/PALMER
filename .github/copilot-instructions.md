<!-- Copilot / AI agent instructions for the PALMER repo -->
# PALMER — Copilot Instructions (concise)

Goal: Help AI coding agents become productive quickly in this repository.

- **Primary build**: the top-level `Makefile` builds the `PALMER` binary. It expects htslib available via `pkg-config`.
  - Build: `make` (or ensure `pkg-config` finds htslib; see `makefile` for error message).
  - If you need more debugging output, set `CXXFLAGS` before invoking `make`, e.g.
    - `export CXXFLAGS="-g -O0 -fsanitize=address,undefined" && make`

- **Key entry points & structure**
  - `PALMER.cpp` — the primary program built by `make`.
  - `scp/` — a pipeline-style folder containing numbered components (e.g. `1_samtools.cpp`, `3_read_masker.cpp`, ..., `8_calling.cpp`) plus `common.hpp` and `tube.cpp`. These implement the stepwise processing used by the tool.
  - `lib/` — fasta libraries used by the tool (e.g. `AluY.fasta`, `L1.3.fasta`).
  - `index/` and `example/` — index resources and small example inputs you can use for quick manual testing.

- **Coding patterns and conventions to respect**
  - The codebase uses manual dynamic arrays and `new[]`/`delete[]` frequently (see `scp/8_calling.cpp`). Avoid sweeping refactors that replace raw arrays with STL containers in many files at once.
  - Many algorithms rely on global-style 2D arrays (`loc`, `loc_TP`, `loc_tsd`) and integer magic constants (`S`, `L`, `BIN_buff`, `le_3`, `le_5`). Search for definitions in `scp/common.hpp` before changing these constants.
  - I/O is typically done with `std::ofstream` objects named like `file2`, `file4` inside the `scp` logic; these are relied on for downstream steps — preserve their output format when editing logic.
  - Conversion between ints and strings is done with `std::stringstream` repeatedly; string concatenation forms identifiers like `info[j][0] + "." + loc_0 + ...` — be careful when altering these formats, they are used as keys (see `info_tsd` / `loc_tsd`).

- **Important dataflow notes**
  - The `scp/` programs interact by producing and consuming tabular files and identifiers created from concatenated fields (sequence identifiers + location fields). Changing the delimiter/order will break matching logic across steps.
  - `example/` contains small `.bam` and `.bai` files to exercise the pipeline locally — use these when iterating on code changes.

- **Where to look first for common change targets**
  - `scp/common.hpp` — definitions for constants and shared declarations.
  - `scp/8_calling.cpp` — large, complex clustering and TSD-identification logic. Good to read before modifying higher-level behavior.
  - `scp/tube.cpp` — likely orchestration for invoking the steps; useful to understand expected inputs/outputs.

- **Debugging and testing advice**
  - Use `make` for normal builds. For runtime memory issues, prefer AddressSanitizer via `CXXFLAGS` rather than valgrind on macOS.
  - Run small end-to-end checks with the files in `example/` to verify no format changes broke the pipeline.
  - Grep helpers: `grep -R "loc_tsd\|loc_TP\|info_tsd\|file2<<\|file4<<" -n` helps locate critical matching logic and output producers.

- **Safe-editing rules for AI agents**
  - Don't globally replace raw arrays with STL containers in a single change — do it in small, testable PRs.
  - Preserve the string-identifiers used as keys (the concatenation format) unless you update all consumers across `scp/` and other code that reads those files.
  - When changing numeric thresholds (S, L, BIN_buff, le_3, le_5), search where those values are used to ensure all dependent logic is updated.

- **Examples of repo-specific searches**
  - Find TSD logic: `grep -R "loc_tsd\[" scp/ -n`
  - Find output files & columns: `grep -R "file2<<\|file4<<" scp/ -n`
  - Find constants: `grep -R "#define S\|int S\|le_3\|le_5" -n scp/`

If anything here is unclear or you want additional examples (e.g., a small checklist for making safe refactors in `8_calling.cpp`), tell me which area to expand and I will iterate.
