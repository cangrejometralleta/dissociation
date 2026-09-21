#!/bin/sh
# Builds the Package, then Runs the Demo.
set -eu

cd -- "$(dirname -- "$0")"

PYTHON="${PYTHON:-python3}"
VENV="${VENV:-.venv}"

verify_interpreter() {
	command -v "$PYTHON" >/dev/null 2>&1 || {
		printf '\342\235\214 wanted a Python interpreter named %s.\n' "$PYTHON" >&2
		printf '   PYTHON=/path/to/python ./test.sh\n' >&2
		exit 1
	}
}

build_package() {
	./build.sh
}

venv_python() {
	[ -x "$VENV/Scripts/python.exe" ] &&
		printf '%s/Scripts/python.exe' "$VENV" ||
		printf '%s/bin/python' "$VENV"
}

run_demo() {
	"$(venv_python)" dissociation.py
}

verify_interpreter
build_package
run_demo
