#!/usr/bin/env bash
#
# Install script for HuRNA_test external dependencies:
#   - LAMMPS (29Aug2024)
#   - TMscore
#   - CONTRAfold
#
# Usage:
#   ./install.sh              # build all
#   ./install.sh lammps       # build only LAMMPS
#   ./install.sh TMscore      # build only TMscore
#   ./install.sh contrafold   # build only CONTRAfold
#
# Flags (env vars):
#   JOBS=N        parallel jobs (default: nproc)
#   CLEAN=1       wipe existing build dirs before building
#   GPU=1         enable LAMMPS GPU build (cuda, sm_80)

# Keep the caller's shell intact when this script is sourced.
(
set -Eeuo pipefail

# ---------- pretty output ----------
if [[ -t 1 ]]; then
    C_RED=$'\033[0;31m'
    C_GREEN=$'\033[0;32m'
    C_YELLOW=$'\033[0;33m'
    C_BLUE=$'\033[0;34m'
    C_BOLD=$'\033[1m'
    C_RESET=$'\033[0m'
else
    C_RED=""; C_GREEN=""; C_YELLOW=""; C_BLUE=""; C_BOLD=""; C_RESET=""
fi

log_info()  { printf "%s[INFO]%s  %s\n"  "$C_BLUE"   "$C_RESET" "$*"; }
log_ok()    { printf "%s[ OK ]%s  %s\n"  "$C_GREEN"  "$C_RESET" "$*"; }
log_warn()  { printf "%s[WARN]%s  %s\n"  "$C_YELLOW" "$C_RESET" "$*" >&2; }
log_error() { printf "%s[FAIL]%s  %s\n"  "$C_RED"    "$C_RESET" "$*" >&2; }
log_step()  { printf "\n%s==>%s %s%s%s\n" "$C_BLUE" "$C_RESET" "$C_BOLD" "$*" "$C_RESET"; }

# ---------- error trap ----------
on_error() {
    local exit_code=$?
    local line=${1:-?}
    log_error "Installation failed at line ${line} (exit ${exit_code})."
    log_error "Current step: ${CURRENT_STEP:-unknown}"
    log_error "Working dir : $(pwd)"
    log_error "See the output above for the underlying build error."
    exit "$exit_code"
}
trap 'on_error $LINENO' ERR

CURRENT_STEP=""

# ---------- config ----------
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
CLEAN="${CLEAN:-0}"
GPU="${GPU:-0}"

log_info "Project root : ${ROOT_DIR}"
log_info "Parallel jobs: ${JOBS}"
log_info "Clean build  : ${CLEAN}"
log_info "GPU build    : ${GPU}"

# ---------- prerequisite check ----------
check_tools() {
    log_step "Checking required tools"
    local missing=()
    for tool in cmake make gcc g++ mpirun; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing+=("$tool")
        else
            log_ok "Found ${tool}: $(command -v "$tool")"
        fi
    done
    if (( ${#missing[@]} > 0 )); then
        log_error "Missing required tools: ${missing[*]}"
        log_error "Install the missing tools or add them to PATH, then retry."
        exit 1
    fi

    log_info "Checking g++ for C++14 support..."
    if ! g++ -std=c++14 -x c++ -c -o /dev/null - <<'CPP'
#if !defined(__cplusplus) || __cplusplus < 201402L
#error "C++14 is required"
#endif
#include <memory>
int main() {
    auto value = std::make_unique<int>(14);
    auto identity = [](auto input) { return input; };
    return identity(*value);
}
CPP
    then
        log_error "g++ failed to compile a C++14 test program."
        log_error "Install a working g++ compiler with C++14 support and ensure it is on PATH."
        exit 1
    fi
    log_ok "g++ supports C++14"
}

# ---------- builders ----------
install_lammps() {
    CURRENT_STEP="LAMMPS"
    log_step "Building LAMMPS (29Aug2024)"

    local src="${ROOT_DIR}/external/lammps-29Aug2024"
    local build="${src}/build"

    if [[ ! -d "$src" ]]; then
        log_error "LAMMPS source dir not found: ${src}"
        return 1
    fi

    if [[ "$CLEAN" == "1" && -d "$build" ]]; then
        log_warn "CLEAN=1: removing existing build dir ${build}"
        rm -rf "$build"
    fi

    mkdir -p "$build"
    cd "$build"

    local cmake_args=(
        -C ../cmake/presets/rna.cmake
        -D BUILD_MPI=yes
    )
    if [[ "$GPU" == "1" ]]; then
        log_info "Configuring LAMMPS with GPU (cuda, sm_80)"
        cmake_args+=(-D PKG_GPU=on -D GPU_API=cuda -D GPU_ARCH=sm_80 -D LAMMPS_MACHINE=gpu)
    else
        log_info "Configuring LAMMPS for CPU"
        # cmake_args+=(-D LAMMPS_MACHINE=cpu)
    fi

    log_info "Running cmake configure..."
    cmake "${cmake_args[@]}" ../cmake

    log_info "Compiling LAMMPS with ${JOBS} jobs..."
    cmake --build . -j "${JOBS}"

    log_ok "LAMMPS built: ${build}"
    cd "$ROOT_DIR"
}

install_cifparse() {
    CURRENT_STEP="cifparse-obj"
    log_step "Building cifparse-obj (v7.105)"

    local dir="${ROOT_DIR}/external/cifparse-obj-v7.105-prod-src"
    if [[ ! -d "$dir" ]]; then
        log_error "cifparse source dir not found: ${dir}"
        return 1
    fi

    cd "$dir"
    if [[ "$CLEAN" == "1" ]]; then
        log_warn "CLEAN=1: running make clean"
        make clean || log_warn "make clean returned non-zero (continuing)"
    fi

    log_info "Compiling cifparse with ${JOBS} jobs..."
    make -j "${JOBS}"

    log_ok "cifparse built: ${dir}"
    cd "$ROOT_DIR"
}

install_TMscore() {
    CURRENT_STEP="TMscore"
    log_step "Building TMscore"

    local dir="${ROOT_DIR}/external/TMscore"
    if [[ ! -d "$dir" ]]; then
        log_error "TMscore source dir not found: ${dir}"
        return 1
    fi

    cd "$dir"
    if [[ "$CLEAN" == "1" ]]; then
        log_warn "CLEAN=1: running make clear"
        make clear || log_warn "make clear returned non-zero (continuing)"
    fi

    log_info "Compiling TMscore with ${JOBS} jobs..."
    make -j "${JOBS}"

    log_ok "TMscore built: ${dir}"
    cd "$ROOT_DIR"
}

install_contrafold() {
    CURRENT_STEP="CONTRAfold"
    log_step "Building CONTRAfold"

    local dir="${ROOT_DIR}/external/contrafold/src"
    if [[ ! -d "$dir" ]]; then
        log_error "CONTRAfold source dir not found: ${dir}"
        return 1
    fi

    cd "$dir"
    if [[ "$CLEAN" == "1" ]]; then
        log_warn "CLEAN=1: running make clean"
        make clean || log_warn "make clean returned non-zero (continuing)"
    fi

    log_info "Compiling CONTRAfold with ${JOBS} jobs..."
    make -j "${JOBS}"

    log_ok "CONTRAfold built: ${dir}"
    cd "$ROOT_DIR"
}

install_cobmms() {
    # parallel building may cause bugs
    CURRENT_STEP="Cobmms"
    log_step "Building Cobmms"

    local dir="${ROOT_DIR}/src"

    cd "$dir"
    if [[ "$CLEAN" == "1" ]]; then
        log_warn "CLEAN=1: running make clean"
        make clean || log_warn "make clean returned non-zero (continuing)"
    fi

    log_info "Compiling Cobmms ..."
    make all

    log_ok "Cobmms built: ${dir}"
    cd "$ROOT_DIR"
}

# ---------- dispatch ----------
START_TS=$(date +%s)

# ---------- extract external dependencies ----------
log_step "Extracting external.zip"
unzip -o "${ROOT_DIR}/external.zip" -d "${ROOT_DIR}"

check_tools

targets=("${@:-all}")
declare -a built=()
declare -a failed=()

run_target() {
    local name="$1"
    if "$name"; then
        built+=("${CURRENT_STEP}")
    else
        failed+=("${CURRENT_STEP}")
        log_error "${CURRENT_STEP} failed — continuing with remaining targets."
    fi
}

# Disable ERR trap for selective failures so we can show a summary.
set +e
trap - ERR

for t in "${targets[@]}"; do
    case "$t" in
        all)
            run_target install_lammps
            run_target install_TMscore
            run_target install_contrafold
            run_target install_cobmms
            ;;
        lammps)     run_target install_lammps ;;
        TMscore)   run_target install_TMscore ;;
        contrafold) run_target install_contrafold ;;
        cobmms)     run_target install_cobmms;;
        -h|--help)
            sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            log_error "Unknown target: $t"
            log_error "Valid targets: all | lammps | TMscore | contrafold"
            exit 2
            ;;
    esac
done

# ---------- extract data ----------
log_step "Extracting data.zip"
merged_dir=$(mktemp -d "${ROOT_DIR}/.data-merge.XXXXXX")

if [[ -z "${merged_dir}" ]]; then
    failed+=("data.zip extraction")
    log_error "Could not create a temporary directory for data.zip."
else
    if ! mkdir -p "${ROOT_DIR}/data" ||
       ! zip -s 0 "${ROOT_DIR}/data.zip" --out "${merged_dir}/data.zip" ||
       ! unzip -o "${merged_dir}/data.zip" -d "${ROOT_DIR}"; then
        failed+=("data.zip extraction")
        log_error "data.zip merge or extraction failed."
    fi

    rm -f "${merged_dir}/data.zip"
    rmdir "${merged_dir}"
fi

# ---------- summary ----------
END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

echo
log_step "Installation summary"
printf "  Elapsed: %ds\n" "$ELAPSED"
if (( ${#built[@]} > 0 )); then
    printf "  %sBuilt:%s   %s\n" "$C_GREEN" "$C_RESET" "${built[*]}"
fi
if (( ${#failed[@]} > 0 )); then
    printf "  %sFailed:%s  %s\n" "$C_RED" "$C_RESET" "${failed[*]}"
    exit 1
fi

log_ok "All targets built successfully."
)
