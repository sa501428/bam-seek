#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_directory="$(cd "${script_directory}/.." && pwd)"
build_directory="${BAM_SEEK_BUILD_DIRECTORY:-${source_directory}/build-package-macos}"
igvcpp_source_directory="${IGVCPP_SOURCE_DIR:-${source_directory}/../igv-cpp}"
qt_prefix="${QT_ROOT:-}"

if [[ -z "${qt_prefix}" ]]; then
    if ! command -v brew >/dev/null 2>&1; then
        echo "Set QT_ROOT to the Qt installation prefix, or install Qt with Homebrew." >&2
        exit 1
    fi
    qt_prefix="$(brew --prefix qt)"
fi

cmake \
    -S "${source_directory}" \
    -B "${build_directory}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH="${qt_prefix}" \
    -DIGVCPP_SOURCE_DIR="${igvcpp_source_directory}" \
    -DBUILD_TESTING=OFF \
    -DBAM_SEEK_ENABLE_PACKAGING=ON

cmake --build "${build_directory}" --config Release --parallel
cpack --config "${build_directory}/CPackConfig.cmake" -C Release -G DragNDrop

echo "DMG and checksum written to ${build_directory}/packages"
