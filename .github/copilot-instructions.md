<!-- Copilot / AI agent instructions for PALMER -->
# PALMER — Quick AI Agent Guide

Purpose: get an AI coding agent productive fast — build, test, and safely modify the PALMER pipeline.

- **Build**: top-level `Makefile` builds the `PALMER` binary. Ensure `htslib` is installable via `pkg-config`.
  - Quick: `make`
  - Debug build (AddressSanitizer):
    - `export CXXFLAGS="-g -O0 -fsanitize=address,undefined" && make`
  - Verify: `pkg-config --cflags --libs htslib` must succeed before `make`.

- **High-level architecture**:
  - `PALMER.cpp` — main program and CLI wrapper.
  - `scp/` — the core pipeline steps implemented as numbered components (`1_samtools.cpp`, `3_read_masker.cpp`, ... `8_calling.cpp`) plus `common.hpp` and `tube.cpp` (orchestration).
  - `lib/` — reference fasta assets (e.g. `AluY.fasta`, `L1.3.fasta`).
  - `index/` and `example/` — pre-built indexes and small example BAMs for manual testing.

- **Dataflow & conventions (critical)**:
  - `scp/` steps exchange tabular files and string identifiers built by concatenation (sequence id + location fields). These identifiers are treated as keys across steps — do NOT change delimiter or ordering without updating all consumers.
  - I/O often uses `std::ofstream` names like `file2` and `file4`; their column order and file names are relied on downstream.
  - Many algorithms use raw `new[]`/`delete[]` arrays and global 2D arrays (`loc`, `loc_TP`, `loc_tsd`). Constants (e.g., `S`, `L`, `BIN_buff`, `le_3`, `le_5`) are defined in `scp/common.hpp` and affect behavior widely.

- **Where to read first**:
  - `scp/common.hpp` — constants and shared structures.
  - `scp/8_calling.cpp` — complex clustering & TSD logic; read before changing calling behavior.
  - `scp/tube.cpp` — orchestration; shows expected inputs/outputs and execution order.

- **Repo-specific grep examples**:
  - TSD logic: `grep -R "loc_tsd\[" scp/ -n`
  - Output producers: `grep -R "file2<<\|file4<<" scp/ -n`
  - Constants: `grep -R "#define S\|int S\|le_3\|le_5" scp/ -n`

- **External dependencies & runtime**:
  - `htslib` (pkg-config), `ncbi-blast+` **2.10.0** recommended (older versions cause bugs).
  - Example builds and runs are in `README.md` and `example/` (use `--intermediate 1` to keep intermediate files for debugging).

- **Editing guidance for agents**:
  - Make minimal, focused changes. Avoid repo-wide replacements (e.g., raw arrays -> STL) in a single PR.
  - Preserve identifier concatenation formats and tabular column orders unless you update ALL producers/consumers.
  - When adjusting numeric thresholds or constants, search usages across `scp/` and run the pipeline on `example/` to validate.

- **Debug & test tips**:
  - Use `export CXXFLAGS="-g -O0 -fsanitize=address,undefined" && make` and run small inputs from `example/`.
  - To examine intermediate outputs, run with `--intermediate 1` and review files produced by `scp/` steps.

If you want, I can add a short checklist for safely refactoring `scp/8_calling.cpp` or create example mutation PRs. Tell me what to expand.
