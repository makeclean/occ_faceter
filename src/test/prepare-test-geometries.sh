#!/bin/bash

set -euo pipefail

oc_workflow="$OVERLAP_CHECKER/tests/demo_workflow.sh"

mkdir -p workdir
cd workdir

for f in ../geometries/*.stp; do
  bash "$oc_workflow" "$f"
done

cp *-merged.brep ../geometries
