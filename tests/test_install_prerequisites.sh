#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
temp_dir="$(mktemp -d)"
trap 'rm -rf -- "$temp_dir"' EXIT

project="${temp_dir}/project"
tools="${temp_dir}/tools"
build_log="${temp_dir}/build.log"
sentinel="${project}/external/lammps-29Aug2024/build/keep"
real_cxx="$(command -v g++)" || {
    printf 'FAIL: g++ is required to run these tests\n' >&2
    exit 1
}

mkdir -p "$(dirname "$sentinel")" "${project}/external/contrafold/src" "$tools"
cp "${root}/install.sh" "${project}/install.sh"
printf 'keep\n' > "$sentinel"

for tool in dirname date mkdir as sed; do
    ln -s "$(command -v "$tool")" "${tools}/${tool}"
done
ln -s /bin/true "${tools}/gcc"
ln -s "$real_cxx" "${tools}/g++"
cat > "${tools}/build-command" <<'SH'
#!/bin/sh
printf '%s\n' "${0##*/}" >> "$BUILD_LOG"
SH
chmod +x "${tools}/build-command"
ln -s "${tools}/build-command" "${tools}/cmake"
ln -s "${tools}/build-command" "${tools}/make"

fail() {
    printf 'FAIL: %s\n%s\n' "$1" "${output:-}" >&2
    exit 1
}

run_installer() {
    local target="$1" clean="${2:-0}"
    if output="$(PATH="$tools" JOBS=1 CLEAN="$clean" GPU=0 \
        BUILD_LOG="$build_log" REAL_CXX="$real_cxx" \
        /bin/bash "${project}/install.sh" "$target" 2>&1)"; then
        status=0
    else
        status=$?
    fi
}

run_sourced_installer() {
    local target="$1" clean="${2:-0}" errexit="${3:-0}"
    if output="$(PATH="$tools" JOBS=1 CLEAN="$clean" GPU=0 \
        BUILD_LOG="$build_log" REAL_CXX="$real_cxx" PROJECT_DIR="$project" CALLER_ERREXIT="$errexit" \
        /bin/bash --noprofile --norc -c '
            if [[ "$CALLER_ERREXIT" == 1 ]]; then
                set -e
            else
                set +e
            fi
            trap : ERR
            log_info() { :; }
            caller_dir=$PWD
            caller_flags=$-
            caller_pipefail=$(set +o pipefail)
            caller_trap=$(trap -p ERR)
            caller_log_info=$(declare -f log_info)

            if [[ "$CALLER_ERREXIT" == 1 ]]; then
                install_status=0
                source "$PROJECT_DIR/install.sh" "$1" || install_status=$?
            else
                source "$PROJECT_DIR/install.sh" "$1"
                install_status=$?
            fi

            if [[ "$PWD" != "$caller_dir" || "$-" != "$caller_flags" ||
                  "$(set +o pipefail)" != "$caller_pipefail" ||
                  "$(trap -p ERR)" != "$caller_trap" ||
                  "$(declare -f log_info)" != "$caller_log_info" ]]; then
                printf "Calling shell state changed.\n" >&2
                exit 1
            fi
            printf "CALLER_SURVIVED=%d\n" "$install_status"
        ' _ "$target" 2>&1)"; then
        status=0
    else
        status=$?
    fi
}

assert_status() {
    [[ "$status" -eq "$1" ]] || fail "Expected exit $1, got $status"
}

assert_output() {
    [[ "$output" == *"$1"* ]] || fail "Expected output to contain: $1"
}

assert_no_build() {
    [[ ! -e "$build_log" ]] || fail "A build command ran before prerequisite checks passed"
    [[ -f "$sentinel" ]] || fail "A build directory was cleaned before prerequisite checks passed"
}

run_installer all 1
assert_status 1
assert_output "Missing required tools: mpirun"
assert_output "Install the missing tools or add them to PATH"
assert_no_build

run_sourced_installer all 1
assert_status 0
assert_output "Missing required tools: mpirun"
assert_output "CALLER_SURVIVED=1"
assert_no_build

run_sourced_installer all 1 1
assert_status 0
assert_output "Missing required tools: mpirun"
assert_output "CALLER_SURVIVED=1"
assert_no_build

rm "${tools}/g++"
run_installer contrafold
assert_status 1
assert_output "Missing required tools: g++ mpirun"
assert_no_build

ln -s /bin/true "${tools}/mpirun"
run_installer contrafold
assert_status 1
assert_output "Missing required tools: g++"
assert_no_build

cat > "${tools}/g++" <<'SH'
#!/bin/sh
printf '%s\n' 'simulated compiler failure' >&2
exit 1
SH
chmod +x "${tools}/g++"
run_installer lammps 1
assert_status 1
assert_output "simulated compiler failure"
assert_output "g++ failed to compile a C++14 test program."
assert_output "C++14 support"
assert_no_build

run_sourced_installer lammps 1
assert_status 0
assert_output "simulated compiler failure"
assert_output "g++ failed to compile a C++14 test program."
assert_output "CALLER_SURVIVED=1"
assert_no_build

cat > "${tools}/g++" <<'SH'
#!/bin/sh
exec "$REAL_CXX" "$@" -std=c++11
SH
run_installer contrafold
assert_status 1
assert_output "C++14 is required"
assert_output "g++ failed to compile a C++14 test program."
assert_no_build

rm "${tools}/g++"
ln -s "$real_cxx" "${tools}/g++"

run_sourced_installer --help
assert_status 0
assert_output "Usage:"
assert_output "CALLER_SURVIVED=0"
assert_no_build

run_sourced_installer invalid-target
assert_status 0
assert_output "Unknown target: invalid-target"
assert_output "CALLER_SURVIVED=2"
assert_no_build

if tty_command="$(command -v script)"; then
    if output="$(PATH="$tools" SHELL=/bin/sh JOBS=1 CLEAN=0 GPU=0 \
        BUILD_LOG="$build_log" REAL_CXX="$real_cxx" PROJECT_DIR="$project" \
        "$tty_command" -q -e -c '/bin/bash "$PROJECT_DIR/install.sh" --help' /dev/null 2>&1)"; then
        status=0
    else
        status=$?
    fi
    assert_status 0
    assert_output "g++ supports C++14"
    assert_output "Usage:"
    assert_no_build
else
    printf 'SKIP: terminal test requires the script command.\n'
fi

run_installer lammps
assert_status 0
assert_output "Found mpirun:"
assert_output "g++ supports C++14"
assert_output "LAMMPS built:"
[[ "$(<"$build_log")" == $'cmake\ncmake' ]] || fail "Expected LAMMPS configure and build commands"

run_installer contrafold
assert_status 0
assert_output "g++ supports C++14"
assert_output "CONTRAfold built:"
[[ "$(<"$build_log")" == $'cmake\ncmake\nmake' ]] || fail "Expected CONTRAfold make command"

run_sourced_installer contrafold
assert_status 0
assert_output "CONTRAfold built:"
assert_output "CALLER_SURVIVED=0"
[[ "$(<"$build_log")" == $'cmake\ncmake\nmake\nmake' ]] || fail "Expected sourced CONTRAfold make command"

printf 'All installer prerequisite and shell-isolation tests passed.\n'
