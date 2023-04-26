#!/bin/bash

set -euo pipefail

oc_workflow="../overlap_checker_workflow.sh"

mkdir -p workdir
cd workdir

for f in ../geometries/*.stp; do
  bash "$oc_workflow" "$f"
done

cp *-merged.brep ../geometries
