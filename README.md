# BAM Seek

BAM Seek is a cross-platform C++/Qt desktop application for targeted inspection of BAM/CRAM evidence. It is being built for reproducible clinical-oriented review workflows, but is not a validated diagnostic device.

## Current capability

This first implementation has a Qt desktop shell and a reusable C++ evidence engine for indexed BAM/CRAM files. It accepts one query per line:

```text
chr7:140453136 A>T
chr1:100000-101000
```

Variant coordinates are one-based. The engine evaluates SNVs, MNVs, and VCF-style small indels against aligned reads; reports ref/alt/other read counts, alt forward/reverse counts, grouped molecule counts, and per-read summaries. Molecule grouping can use raw reads, an auto-detected `MI`, `RX`, or `UB` tag, or a selected tag.

Four-column `CHROM POS REF ALT` and standard VCF rows are also accepted.

An existing BAM/CRAM index is required. CRAM additionally needs its appropriate reference FASTA. The reference field is intended for an hg19 FASTA in the main release.

Region scans require an indexed reference FASTA and currently discover small SNV candidates in windows up to 100 kb. Small indel discovery and a graphical pileup are planned follow-on work; the current per-read text panel supplies inspectable evidence.

## Build

The default build expects this checkout beside `igv-cpp` and requires CMake 3.24+, C++20, Qt 6.5+ Widgets/Concurrent, HTSlib, zlib, and pkg-config.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

If `igv-cpp` is elsewhere, pass `-DIGVCPP_SOURCE_DIR=/path/to/igv-cpp`.

Audit JSON exports record a UTC timestamp, SHA-256 file fingerprints for the alignment/index/reference, filters, version, query metrics, and individual read evidence.
