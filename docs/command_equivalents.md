# External command lines used in helper-based calls

PALMER now prefers the helper functions defined in `scp/common.hpp` (such as `run_process` and `stream_process_output`) instead of raw `system()` calls. The underlying command lines invoked by these helpers correspond to the following shell-equivalent invocations:

- **Samtools region extraction (`scp/1_samtools.cpp`)**
  - Command: `samtools view -q <MAPQ> -F 0x100 -F 0x200 -F 0x400 -T <fasta> <input_bam> <chr>:<start>-<end>`
  - Output handling: streamed directly into `region.sam` with whitespace normalized to underscores.

- **Blastn run (`scp/4_blastn.cpp`)**
  - Command: `blastn -evalue 0.001 -task blastn -gapopen 4 -max_target_seqs 10000000 -query <query_fasta> -subject <WD_dir>/SEQ.masked -outfmt "7 qacc sacc evalue qstart qend sstart send"`
  - Output handling: written to `<WD_dir>/blastn.txt` while suppressing BLAST usage reporting via `BLAST_USAGE_REPORT=0`.

These commands are executed without shell wrapping; arguments are passed directly to the binaries through `execvp` in the helper functions.
