# BAM Seek

BAM Seek is a focused C++/Qt desktop application for calculating variant allele frequencies across a set of local or remote BAM files. It is intended for reproducible evidence review, but is not a validated diagnostic device.

The workspace supports persistent light and dark appearances through the small theme toggle at the upper right.

## Workflow

1. Add, remove, or clear local and HTTPS BAMs in the dedicated **BAM sources** panel. Check at most one prepared BAM to designate it as the **Current** sample; every unchecked BAM is **Historical**. If none is checked, all BAMs are historical.
2. BAMs are prepared asynchronously when added: BAM Seek opens the alignment, reads its header, and loads its index so later queries start quickly. It does not calculate read or molecule evidence until **Calculate VAFs** is clicked.
3. For local BAMs, BAM Seek automatically looks beside each file for a `.bai` or `.csi` index. For HTTPS BAMs, the reader performs conventional remote index discovery. There is no manual index field and the application does not create or modify indexes.
4. Enter one variant per line in either the **Current variants** or **Historical variants** box using any of these forms:

```text
chr7:140453136 A>T
chr7:g.140453136A>T
chr7 140453136 A T
chr7 140453136 . A T
BRAF c.1799T>A p.V600E chr7:140453136 A>T
IGV 9 139390861 NOTCH1 stop_gained c.7330C>T p.Q2444* G A 3.1
5 176943930 DDX41 ENST00000507955.1 MISSENSE c.17C>T p.P6L G A GTTCGXGTTCC 4.4 19 45
```

Coordinates are one-based. Clinical labels may be included for display, but must be accompanied by a genomic coordinate and REF/ALT allele; BAM Seek does not use a clinical-mapping table or a reference FASTA.

For report-style rows, `IGV` may be at the start of the row, on a separate line, or omitted. An optional transcript column is recognized before the variant type. BAM Seek extracts chromosome, position, gene, type, coding change, protein change, REF, and ALT, and ignores every trailing column—including a previously reported VAF. Insertions, deletions, and duplications must use left-anchored REF/ALT alleles: the position is the shared normal anchor base, and the inserted or deleted sequence follows it.

Repeated variants are evaluated once. BAM Seek considers annotations duplicates when their gene matches and either their coding or protein notation matches; unannotated genomic inputs fall back to chromosome, position, REF, and ALT. A variant entered in both boxes is labeled **Both** in the results and excluded from the comparison narrative.

Each BAM and variant receives a separate result row containing:

- ALT read count and informative read depth
- read VAF: `ALT reads / (REF reads + ALT reads)`
- ALT molecule count and informative molecule depth
- molecule VAF: `ALT molecules / (REF molecules + ALT molecules)`
- OTHER/N reads and ambiguous paired fragments, which are reported but excluded from the VAF denominators
- a copy-ready current-versus-historical comparison narrative

The evidence table labels BAM and variant roles with separate colors, displays compact BAM filenames with the complete source in a tooltip, places a selected current BAM first, and otherwise sorts BAMs and variants alphabetically. The copy-ready comparison narrative evaluates current-only variants against every historical BAM, separates molecularly detected and absent variants, and evaluates historical-only variants against the current BAM when one is selected.

Molecules are paired fragments, not UMI families. Alignments with the same read name are grouped together, and a REF or ALT molecule call requires an unambiguous majority among that fragment's callable alignments. Ties and conflicting calls are counted as ambiguous molecules.

The optional mapQ and baseQ controls filter low-quality alignments and allele observations. Duplicates, secondary alignments, and supplementary alignments are excluded by default.

## Broadcast receiver

The **Broadcast receiver** has its own top-level tab, separate from the BAM and variant workspace. BAM Seek listens on `127.0.0.1:60151`, inspects every received command for the review inbox, and forwards the complete byte stream unchanged to IGV on `127.0.0.1:60152`. IGV's response is relayed unchanged to the original sender. For a `GET /load` request, BAM Seek also extracts only the decoded `file` value and prepares that BAM in the active BAM panel; locus, genome, index, merge, and other request fields do not change the analysis workspace, but are still forwarded to IGV.

Configure IGV's command listener for port `60152` in **View → Preferences → Advanced**, or launch command-line IGV with `-p 60152`. Keep callers and existing IGV links pointed at port `60151`. Both relay endpoints bind only to localhost.

## Read-only operation

VAF calculation is read-only and runs on the local machine. BAM Seek does not upload data, create indexes, modify BAMs, or export audit files. Remote servers necessarily receive HTTPS range requests while BAMs are prepared and queried, so only approved remote sources should be used. HTTP and other URL schemes are rejected. The copy-ready result paragraph is the primary export surface and is placed on the system clipboard only when requested.

## Build

The default build expects this checkout beside `igv-cpp` and requires CMake 3.24+, a C++20 compiler, Qt 6.5+ (Widgets, Concurrent, and Network), HTSlib, zlib, and pkg-config.

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

The executable will be under `build/Release`. Run `windeployqt` from the matching Qt installation before packaging it in an installer.

If `igv-cpp` is elsewhere, pass `-DIGVCPP_SOURCE_DIR=/path/to/igv-cpp`. For local memory-safety testing, configure with `-DBAM_SEEK_ENABLE_SANITIZERS=ON` on Clang or GCC.
