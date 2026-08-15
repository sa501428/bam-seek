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

For report-style rows, `IGV` may be at the start of the row, combined with the chromosome in the first tab-separated field (for example, `IGV X`), on a separate line, or omitted. An optional transcript column is recognized before the variant type, and the protein-change column may be blank for intronic or splice variants. BAM Seek extracts chromosome, position, gene, type, coding change, protein change when present, REF, and ALT, and ignores every trailing column—including a previously reported VAF. Insertions, deletions, and duplications must use left-anchored REF/ALT alleles: the position is the shared normal anchor base, and the inserted or deleted sequence follows it. Equal-length substitutions and delins alleles are processed as listed.

Repeated variants are evaluated once. BAM Seek considers annotations duplicates when their gene matches and either their coding or protein notation matches; unannotated genomic inputs fall back to chromosome, position, REF, and ALT. A variant entered in both boxes is labeled **Both** in the results and excluded from the comparison narrative.

Each BAM and variant receives a separate result row containing:

- ALT read count and informative read depth
- read VAF: `ALT reads / (REF reads + ALT reads)`
- ALT molecule count and informative molecule depth
- molecule VAF: `ALT molecules / (REF molecules + ALT molecules)`
- OTHER/N reads and ambiguous paired fragments, which are reported but excluded from the VAF denominators
- a copy-ready current-versus-historical comparison narrative

When a named gene contains two or more distinct variants, BAM Seek also evaluates every pair automatically in every loaded BAM. The read-based phasing table and copy-ready narrative report the genomic distance, number of informative independent molecules, `AB`, `Ab`, `aB`, and `ab` molecule counts, and one of four results:

- **Cis** when sufficient ALT/ALT (`AB`) molecules link both variants with no more than the configured conflict fraction.
- **Trans** when sufficient reciprocal ALT/REF (`Ab`) and REF/ALT (`aB`) molecules link both variants with no more than the configured conflict fraction.
- **Indeterminate** when molecules span both loci but support is under the selected threshold or the evidence conflicts.
- **Too far apart to phase** when no high-quality molecule directly observes both variant positions. The row tooltip gives the exact reason, so this result also makes insufficient physical linkage explicit when distance alone is not the cause.

Phasing uses the same mapQ, baseQ, alignment-flag, duplicate, and molecule-grouping rules as VAF calculation. A read pair counts as one fragment, and separate mates may jointly observe the two variants only when they form a proper pair. A single read may also observe both. The inferred insert interval between paired reads is never treated as sequenced: each locus must have an actual high-quality REF or ALT observation. MI/RX/UB or a selected molecule tag is honored when configured. Only direct pairwise molecular phase is reported; phase is not propagated indirectly through an intermediate variant.

The evidence table labels BAM and variant roles with separate colors, displays compact BAM filenames with the complete source in a tooltip, places a selected current BAM first, and otherwise sorts BAMs and variants alphabetically. The copy-ready comparison narrative reports consensus reads as `ALT / (ALT + REF)`, evaluates current-only variants against every historical BAM, separates molecularly detected and absent variants, and evaluates historical-only variants against the current BAM when one is selected.

Successful evidence calculations are cached in memory by BAM, normalized genomic allele, and all calculation filters. Changing which BAM is current immediately reuses and re-labels the displayed evidence without querying any BAM again. Re-running after adding a BAM calculates only uncached BAM–variant pairs. Removing a BAM erases its cache entries, clearing the BAM panel erases the whole cache, and exiting BAM Seek destroys the cache; it is never written to disk.

Molecules are paired fragments, not UMI families. Alignments with the same read name are grouped together, and a REF or ALT molecule call requires an unambiguous majority among that fragment's callable alignments. Ties and conflicting calls are counted as ambiguous molecules.

The **Settings** tab contains the mapQ and baseQ controls, alignment-flag inclusion options, molecule grouping mode, minimum phase-supporting molecule count, and maximum conflicting-molecule percentage. The phasing defaults require two supporting molecules (two `AB` molecules for cis or two in each reciprocal direction for trans) and permit at most 10% conflicting alternate-bearing molecules. Duplicates, secondary alignments, and supplementary alignments are excluded by default. Molecules can be paired by read name, grouped by an automatically detected `MI`, `RX`, or `UB` tag, or grouped by a user-selected two-character BAM tag. These calculation settings are saved and apply to the next run. Hovering over a variant or evidence count shows how many overlapping alignments were excluded by alignment/mapQ filters, lacked a callable allele, or fell below the baseQ threshold. For an insertion, baseQ is evaluated across the left anchor and every inserted base; consequently, IGV may display an insertion alignment that is intentionally excluded from the numeric VAF at the selected baseQ threshold.

## Settings and broadcast receiver

The broadcast receiver and its activity log remain available in **Settings**, separate from the BAM and variant workspace. The receiving and IGV forwarding ports are editable and persisted; their defaults are `127.0.0.1:60151` and `127.0.0.1:60152`, respectively. BAM Seek inspects every received command for the activity log and forwards the complete byte stream unchanged to IGV. IGV's response is relayed unchanged to the original sender. For a `GET /load` request, BAM Seek also extracts only the decoded `file` value and prepares that BAM in the active BAM panel; locus, genome, index, merge, and other request fields do not change the analysis workspace, but are still forwarded to IGV.

Configure IGV's command listener for the selected forwarding port (default `60152`) in **View → Preferences → Advanced**, or launch command-line IGV with `-p PORT`. Keep callers and existing IGV links pointed at the selected receiving port (default `60151`). Both relay endpoints bind only to localhost.

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
open "build/BAM Seek.app"
```

Build a self-contained DMG (including Qt and Homebrew runtime libraries) with:

```sh
./scripts/package-macos.sh
```

The DMG and its SHA-256 checksum are written to `build-package-macos/packages`.
Set `QT_ROOT`, `IGVCPP_SOURCE_DIR`, or `BAM_SEEK_BUILD_DIRECTORY` to override
the script defaults. The bundle is ad-hoc signed after its libraries are fixed
up so it can run locally. Distribution outside a local or managed environment
still requires Apple Developer ID signing and notarization.

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

The executable will be under `build/Release`. To produce a self-contained NSIS
`.exe` installer, a WiX `.msi`, or both, install the packaging tool(s):

```powershell
choco install nsis
dotnet tool install --global wix --version 4.0.4
wix extension add --global WixToolset.UI.wixext/4.0.4
```

Then run one of:

```powershell
./scripts/package-windows.ps1                         # both .exe and .msi
./scripts/package-windows.ps1 -InstallerFormat NSIS  # .exe only
./scripts/package-windows.ps1 -InstallerFormat WIX   # .msi only
```

Set `QT_ROOT` to the Qt MSVC directory (for example,
`C:/Qt/6.8.3/msvc2022_64`) and `VCPKG_ROOT` to the vcpkg checkout first. The
installer(s) and SHA-256 checksums are written to
`build-package-windows/packages`. The install step collects HTSlib and its
vcpkg runtime DLLs, while Qt's deployment helper adds the required Qt DLLs,
platform plugin, image handlers, and TLS backends. CMake 3.30 or newer is
required for the scripted WiX 4 build; CMake 3.24–3.29 can instead use an
installed WiX 3 toolset. The generated Windows installers are not Authenticode
signed, so public distribution may trigger a Microsoft Defender SmartScreen
warning until a trusted signing certificate is added.

The same packages can be generated without the helper scripts after configuring
and building:

```sh
cpack --config build/CPackConfig.cmake -C Release -G DragNDrop  # macOS
```

```powershell
cpack --config build/CPackConfig.cmake -C Release -G NSIS  # Windows .exe
cpack --config build/CPackConfig.cmake -C Release -G WIX   # Windows .msi
```

If `igv-cpp` is elsewhere, pass `-DIGVCPP_SOURCE_DIR=/path/to/igv-cpp`. For local memory-safety testing, configure with `-DBAM_SEEK_ENABLE_SANITIZERS=ON` on Clang or GCC.
