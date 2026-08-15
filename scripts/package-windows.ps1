[CmdletBinding()]
param(
    [ValidateSet("Both", "NSIS", "WIX")]
    [string]$InstallerFormat = "Both",
    [string]$QtRoot = $env:QT_ROOT,
    [string]$VcpkgRoot = $env:VCPKG_ROOT,
    [string]$IgvCppSourceDirectory = $env:IGVCPP_SOURCE_DIR,
    [string]$BuildDirectory = ""
)

$ErrorActionPreference = "Stop"
$sourceDirectory = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot ".."))

if ([string]::IsNullOrWhiteSpace($BuildDirectory)) {
    $BuildDirectory = Join-Path $sourceDirectory "build-package-windows"
}
if ([string]::IsNullOrWhiteSpace($IgvCppSourceDirectory)) {
    $IgvCppSourceDirectory = [System.IO.Path]::GetFullPath((Join-Path $sourceDirectory "../igv-cpp"))
}
if ([string]::IsNullOrWhiteSpace($QtRoot)) {
    throw "Set QT_ROOT or pass -QtRoot (for example C:/Qt/6.8.3/msvc2022_64)."
}
if ([string]::IsNullOrWhiteSpace($VcpkgRoot)) {
    throw "Set VCPKG_ROOT or pass -VcpkgRoot."
}

$toolchainFile = Join-Path $VcpkgRoot "scripts/buildsystems/vcpkg.cmake"
if (-not (Test-Path $toolchainFile)) {
    throw "The vcpkg toolchain was not found at $toolchainFile"
}

$windowsGenerators = switch ($InstallerFormat) {
    "NSIS" { "NSIS" }
    "WIX" { "WIX" }
    default { "NSIS;WIX" }
}

& cmake `
    -S $sourceDirectory `
    -B $BuildDirectory `
    -G "Visual Studio 17 2022" `
    -A x64 `
    "-DCMAKE_TOOLCHAIN_FILE=$toolchainFile" `
    "-DCMAKE_PREFIX_PATH=$QtRoot" `
    "-DIGVCPP_SOURCE_DIR=$IgvCppSourceDirectory" `
    "-DBAM_SEEK_WINDOWS_GENERATORS=$windowsGenerators" `
    -DBUILD_TESTING=OFF `
    -DBAM_SEEK_ENABLE_PACKAGING=ON
if ($LASTEXITCODE -ne 0) { throw "CMake configuration failed." }

& cmake --build $BuildDirectory --config Release --parallel
if ($LASTEXITCODE -ne 0) { throw "BAM Seek build failed." }

$generators = if ($InstallerFormat -eq "Both") { @("NSIS", "WIX") } else { @($InstallerFormat) }
foreach ($generator in $generators) {
    & cpack --config (Join-Path $BuildDirectory "CPackConfig.cmake") -C Release -G $generator
    if ($LASTEXITCODE -ne 0) { throw "$generator packaging failed." }
}

Write-Host "Installer(s) and checksums written to $(Join-Path $BuildDirectory 'packages')"
