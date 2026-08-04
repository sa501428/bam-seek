# BAM Seek

BAM Seek is a cross-platform C++/Qt desktop application for targeted inspection of BAM/CRAM evidence. It is being built for reproducible clinical-oriented review workflows, but is not a validated diagnostic device.

## Read-only operation and secure remote BAM access

All computation happens on the local machine. BAM Seek accepts local filesystem paths and HTTPS BAM/CRAM resources with byte-range support, so it can read an indexed BAM from an approved remote server without uploading local files or result data. HTTP and other URL schemes are rejected. The app has no telemetry behavior and never creates indexes or modifies BAM/CRAM/reference files. The sole optional write is an audit JSON file that the user explicitly selects and saves locally.

When using a remote resource, the remote server necessarily receives HTTPS requests for the BAM/index and the byte ranges needed for the selected loci. Use only approved, access-controlled servers for clinical data.

## Current capability

This first implementation has a Qt desktop shell and a reusable C++ evidence engine for indexed BAM/CRAM files. It accepts one query per line:

```text
chr7:140453136 A>T
chr1:100000-101000
```

Variant coordinates are one-based. The engine evaluates SNVs, MNVs, and VCF-style small indels against aligned reads; reports ref/alt/other read counts, strand-specific counts, a two-sided Fisher strand-bias p-value, molecule counts when available, and per-read summaries. Molecule grouping can use an auto-detected `MI`, `RX`, or `UB` tag or a selected two-character BAM tag. Auto-detection requires the tag on at least 90% of callable reads; otherwise molecule counts are explicitly unavailable. `RX`/`UB` grouping incorporates fragment coordinates, and `UB` also incorporates `CB` when present.

Four-column `CHROM POS REF ALT` and standard VCF rows are also accepted.

An existing BAM/CRAM index is required. The app supports an explicit local or HTTPS index location when it is not discoverable beside the alignment. CRAM additionally needs its appropriate reference FASTA. The reference field is intended for an hg19 FASTA in the main release. When a reference is configured, targeted REF alleles are checked against it before evidence is reported.

Region scans require an indexed reference FASTA and currently discover small SNV candidates in windows up to 100 kb. Small indel discovery remains follow-on work; the current per-read text panel and graphical pileup provide inspectable evidence.

Targeted variant rows can be opened in the **Pileup** tab. This local rendering colors forward reads blue, reverse reads red, reference mismatches yellow, and indels green; it can also group and link read pairs by name.

## Build

The default build expects this checkout beside `igv-cpp` and requires CMake 3.24+, a C++20 compiler, Qt 6.5+ (Widgets and Concurrent), HTSlib, zlib, and pkg-config.

### macOS

Install the dependencies with Homebrew:

```sh
brew install cmake htslib pkg-config qt
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$(brew --prefix qt)"
cmake --build build -j
ctest --test-dir build --output-on-failure
open build/bam-seek.app
```

For a distributable app bundle, Qt's deployment helper copies the Qt frameworks into the generated bundle:

```sh
"$(brew --prefix qt)/bin/macdeployqt" build/bam-seek.app
```

### Windows

Use Visual Studio 2022, CMake, vcpkg, and Qt 6 (MSVC 64-bit). Install the native dependencies and configure from a Developer PowerShell:

```powershell
vcpkg install htslib:x64-windows zlib:x64-windows
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
  -DCMAKE_PREFIX_PATH="C:/Qt/6.x.x/msvc2022_64"
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The executable will be under `build/Release`. Run `windeployqt` from the matching Qt installation to produce a folder containing the necessary Qt DLLs before packaging it in an installer.

If `igv-cpp` is elsewhere, pass `-DIGVCPP_SOURCE_DIR=/path/to/igv-cpp`.

Audit JSON exports are written atomically in a background task and record a UTC timestamp, filters, version, errors, query metrics, and individual read evidence. Local alignment/index/reference files receive SHA-256 fingerprints and a stability check; remote resources record a credential- and query-redacted HTTPS URI because their bytes are not downloaded in full for hashing.

For local memory-safety testing, configure with `-DBAM_SEEK_ENABLE_SANITIZERS=ON` on Clang or GCC.
