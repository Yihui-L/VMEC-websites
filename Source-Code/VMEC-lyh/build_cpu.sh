#!/usr/bin/env bash
set -euo pipefail
source_dir=$(dirname "$(realpath -e -- "$0")")
work=${1:?Usage: build_cpu.sh OUTSIDE_SOURCE_BUILD_DIRECTORY [debug|release]}
mode=${2:-release}
case "$mode" in debug|release) ;; *) exit 2;; esac
# Check physical paths before creating or copying anything.
work=$(realpath -m -- "$work")
case "$work/" in "$source_dir/"*) echo "Build directory must be outside the source tree" >&2; exit 2;; esac
case "$source_dir/" in "${work%/}/"*) echo "Build directory must not contain the source tree" >&2; exit 2;; esac
mkdir -p -- "$work"
cp -a "$source_dir/." "$work/"
mkdir -p "$work/bin" "$work/LIBSTELL/Debug" "$work/LIBSTELL/Release" "$work/VMEC2000/Debug" "$work/VMEC2000/Release"
export STELLOPT_PATH="$work" MACHINE=vmec_lyh MYHOME="$work/bin"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
jobs=${BUILD_JOBS:-4}
make -C "$work/LIBSTELL" "$mode" -j"$jobs"
case "$mode" in debug) obj=Debug;; release) obj=Release;; esac
while IFS= read -r src; do
    target="$work/LIBSTELL/$obj/$(basename "$src")"
    test -e "$target" || ln -s "$src" "$target"
done < <(find "$work/LIBSTELL/Sources" -type f \( -name '*.f' -o -name '*.f90' \))
archive="$work/LIBSTELL/$obj/libstell.a"
make -C "$work/VMEC2000" "$mode" -j"$jobs" LIB="$archive" \
    LIB_LINK="$archive -lopenblas -lscalapack-openmpi $(nf-config --flibs)"
mpifort -fopenmp -I"$work/LIBSTELL/$obj" -o "$work/bin/vmec_cpu_input" \
    "$work/VMEC2000/cpu_input.f90" "$archive"
cp "$work/VMEC2000/run_vmec.py" "$work/bin/run_vmec.py"
chmod +x "$work/bin/run_vmec.py"
sha256sum "$work/bin/xvmec2000"
