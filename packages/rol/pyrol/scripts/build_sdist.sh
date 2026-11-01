#!/bin/bash
# Build a PyROL source distribution (sdist), as the build-sdist job of
# .github/workflows/python.yml does with this script: run Binder through the
# Trilinos configuration, package ROL, Teuchos and the generated sources with
# CPack, and turn that into an sdist. Works on Linux and macOS.
#
# Usage: pyrol/scripts/build_sdist.sh [options]
#   -w DIR     work directory, outside this checkout; the sdist is written to
#              DIR/dist (default: ./pyrol-sdist)
#   -t DIR     Trilinos source tree to use, outside this checkout (its packages/rol
#              is replaced by this checkout); by default, the latest Trilinos
#              develop is cloned into DIR
#   -b BINDER  Binder executable (default: binder on the PATH)
#   -l PREFIX  on macOS, the LLVM Binder was built with (default: brew --prefix llvm@16)
#   -p PYTHON  Python with pybind11 and build installed (default: python3)
#
# Also needed: CMake, Ninja and BLAS/LAPACK. See pyrol/DEVELOPING.md for Binder.
set -euo pipefail

usage() { sed -n '2,17p' "$0" | sed -E 's/^# ?//'; exit "$1"; }

work=pyrol-sdist
trilinos=""
binder=$(command -v binder || true)
llvm=""
python=python3
while getopts "w:t:b:l:p:h" opt; do
  case $opt in
    w) work=$OPTARG ;;
    t) trilinos=$OPTARG ;;
    b) binder=$OPTARG ;;
    l) llvm=$OPTARG ;;
    p) python=$OPTARG ;;
    h) usage 0 ;;
    *) usage 1 ;;
  esac
done
[ -x "$binder" ] || { echo "Binder not found; pass it with -b (see pyrol/DEVELOPING.md)." >&2; exit 1; }

rol=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$work"
work=$(cd "$work" && pwd)

# This checkout is copied into the Trilinos tree, so neither that tree nor the
# work directory may be inside it.
inside_checkout() { case "$1/" in "$rol"/*) return 0 ;; esac; return 1; }
if inside_checkout "$work"; then
  rmdir "$work" 2>/dev/null || true
  echo "The work directory $work is inside this checkout; pass -w with one outside it." >&2
  exit 1
fi

# Trilinos, with this checkout of ROL as packages/rol.
if [ -z "$trilinos" ]; then
  trilinos=$work/trilinos
  [ -d "$trilinos" ] || git clone --depth 1 https://github.com/trilinos/Trilinos.git "$trilinos"
fi
trilinos=$(cd "$trilinos" && pwd)
if inside_checkout "$trilinos"; then
  echo "The Trilinos tree $trilinos is inside this checkout; pass -t with one outside it." >&2
  exit 1
fi
rsync -a --delete --exclude=.git "$rol/" "$trilinos/packages/rol/"

# The package version is the date and time in Denver.
cp "$rol/pyrol/pyproject.toml" "$trilinos/pyproject.toml"
version="$(TZ=America/Denver date +%Y.%m.%d).dev$(TZ=America/Denver date +%H%M)"
sed -i.bak "s/version = \"VERSION\"/version = \"$version\"/" "$trilinos/pyproject.toml"
rm "$trilinos/pyproject.toml.bak"

# On macOS, Binder needs the SDK and the builtin headers of its LLVM.
binder_flags=()
if [ "$(uname)" = Darwin ]; then
  llvm=${llvm:-$(brew --prefix llvm@16)}
  resource_dir=$(ls -d "$llvm"/lib/clang/* | head -1)
  binder_flags=(-D "PyROL_BINDER_FLAGS=-isysroot;$(xcrun --show-sdk-path);-resource-dir;$resource_dir")
fi

build=$work/build
cmake -G Ninja -S "$trilinos" -B "$build" \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D CPACK_SOURCE_IGNORE_FILES="/packages/rol/example;/packages/rol/test;/packages/rol/tutorial;/packages/teuchos/core/test;/packages/teuchos/core/example;/packages/teuchos/numerics/test;/packages/teuchos/numerics/example;/packages/teuchos/parameterlist/test;/packages/teuchos/parameterlist/example;/packages/teuchos/parser/test;/packages/teuchos/parser/example;/packages/teuchos/remainder/test;/packages/teuchos/remainder/example;/packages/teuchos/comm/test;/packages/teuchos/comm/example" \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D Python3_EXECUTABLE:FILEPATH="$(command -v "$python")" \
  -D ROL_ENABLE_EXAMPLES:BOOL=OFF \
  -D ROL_ENABLE_PYROL:BOOL=ON \
  -D ROL_ENABLE_TESTS:BOOL=OFF \
  -D PYROL_ENABLE_BINDER:BOOL=ON \
  -D PYROL_ENABLE_BINDER_UPDATE:BOOL=ON \
  -D PyROL_BINDER_EXECUTABLE:FILEPATH="$binder" \
  ${binder_flags[@]+"${binder_flags[@]}"} \
  -D Trilinos_CPACK_SOURCE_GENERATOR="TGZ" \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_CPACK_PACKAGING=ON \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
  -D Trilinos_ENABLE_Fortran:BOOL=OFF \
  -D Trilinos_ENABLE_ROL:BOOL=ON \
  -D Trilinos_ENABLE_TESTS:BOOL=OFF
cmake --build "$build" --target package_source

# Unpack the reduced source tree and turn it into an sdist.
rm -rf "${work:?}/source"
mkdir "$work/source"
tar -xzf "$(find "$build" -maxdepth 1 -name '*-Source.tar.gz' | head -1)" -C "$work/source"
source_tree=$(find "$work/source" -mindepth 1 -maxdepth 1 -type d | head -1)
rm -f "$source_tree/packages/framework/asan_assets/dummy_dlclose.so"
"$python" -m build --sdist --outdir "$work/dist" "$source_tree"
echo "sdist: $(ls "$work"/dist/*.tar.gz)"
