# PALMER

PALMER detects non-reference mobile element insertions (LINE, Alu, SVA, HERVK) and other insertions (for example, NUMTs and HPV) from long-read sequencing and assembled contigs.

## Documentation

Full documentation is in the wiki:

- [PALMER Wiki Home](https://github.com/WeichenZhou/PALMER/wiki)
- [Installation](https://github.com/WeichenZhou/PALMER/wiki/Installation)
- [Quickstart](https://github.com/WeichenZhou/PALMER/wiki/Quickstart)
- [Parameters](https://github.com/WeichenZhou/PALMER/wiki/Parameters)
- [Examples](https://github.com/WeichenZhou/PALMER/wiki/Examples)
- [Output Files](https://github.com/WeichenZhou/PALMER/wiki/Output-Files)
- [Recommendations](https://github.com/WeichenZhou/PALMER/wiki/Recommendations)
- [Changelog](https://github.com/WeichenZhou/PALMER/wiki/Changelog)

## Issues

- [Submit a bug report](https://github.com/WeichenZhou/PALMER/issues/new)
- [Browse all issues](https://github.com/WeichenZhou/PALMER/issues)

## Quick Build

```bash
git clone https://github.com/WeichenZhou/PALMER.git
cd PALMER
make
```

Note: install `htslib` development headers and `ncbi-blast++/2.10.0` before building.

## Citation

For PALMER /* where it begins (Pre-mAsking Long reads for Mobile Element inseRtion) */:
- Weichen Zhou et al., [Identification and characterization of occult human-specific LINE-1 insertions using long-read sequencing technology](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1173/5680708), Nucleic Acids Research, 2019, gkz1173, `https://doi.org/10.1093/nar/gkz1173`

For TEnCATS/NanoPal:
- Torrin L. McDonald and Weichen Zhou et al., [Cas9 targeted enrichment of mobile elements using nanopore sequencing](https://www.nature.com/articles/s41467-021-23918-y), Nature Communications, 2021, `https://doi.org/10.1038/s41467-021-23918-y`

For PALMER2.0 multi-functional usage in WGS, capture, assembled contigs, and single-cell data:
- Weichen Zhou et al., [A personalized multi-platform assessment of somatic mosaicism in the human frontal cortex](https://doi.org/10.1101/2024.12.18.629274), bioRxiv, 2024, `https://doi.org/10.1101/2024.12.18.629274`

For SMaHT benchmarking:
- Wang et al., [Multi-platform framework for mapping somatic retrotransposition in human tissues](https://www.biorxiv.org/content/10.1101/2025.10.07.680917), bioRxiv, 2025, `https://doi.org/10.1101/2025.10.07.680917`

## Contact

- arthurz@umich.edu

## Copyright

Copyright (c) 2018-2026 Weichen Arthur Zhou @ UMich&Fudan&HUST
