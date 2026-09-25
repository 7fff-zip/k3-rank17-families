#!/bin/sh
# Runs the checks and compares their output with expected/ (timings stripped).
# Needs gp on PATH, e.g.  nix develop -c ./check.sh
set -eu
cd "$(dirname "$0")"
mkdir -p out
status=0
run() {  # run NAME FILE...: read the files in one fresh gp session
  name=$1; shift
  gp -q --default parisize=256M --default parisizemax=2G "$@" < /dev/null \
    | sed -E 's/ *\([0-9.]+ s\)$//' > "out/$name.out"
  if diff -u "expected/$name.out" "out/$name.out" > "out/$name.diff"; then
    echo "$name: output matches expected/$name.out"
  else
    echo "$name: output DIFFERS from expected/$name.out (see out/$name.diff)"; status=1
  fi
}
run D510 D510_K3_family.gp check.gp
run D546 D546_K3_family.gp check.gp
run D546_to302 D546_K3_family.gp D546_to302.gp
run D546_to845 D546_K3_family.gp D546_to845.gp
exit $status
