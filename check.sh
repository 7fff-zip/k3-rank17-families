#!/bin/sh
# Runs check.gp on both families and compares with expected/ (timings stripped).
# Needs gp on PATH, e.g.  nix develop -c ./check.sh
set -eu
cd "$(dirname "$0")"
mkdir -p out
status=0
for d in D510 D546; do
  gp -q --default parisize=256M --default parisizemax=2G "${d}_K3_family.gp" check.gp \
    | sed -E 's/ *\([0-9.]+ s\)$//' > "out/$d.out"
  if diff -u "expected/$d.out" "out/$d.out" > "out/$d.diff"; then
    echo "$d: output matches expected/$d.out"
  else
    echo "$d: output DIFFERS from expected/$d.out (see out/$d.diff)"; status=1
  fi
done
exit $status
