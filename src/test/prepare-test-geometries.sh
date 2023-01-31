#!/bin/bash

set -euo pipefail

oc_workflow="$OVERLAP_CHECKER/tests/demo_workflow.sh"

mkdir -p workdir
cd workdir

for f in ../geometries/*.stp; do
  bash "$oc_workflow" "$f"
done

cp *-merged.brep ../geometries

exit

STP_FILE=toroidal_field_coils_toroidal_field_coils_case_2.stp
ROOT=`basename -s .stp "$STP_FILE"`

occ_faceter -t 1e-3 "$ROOT-merged.brep" -f ../materials.txt --add_mat_ids

make_watertight dagmc_not_watertight.h5m > "$ROOT.out"

mbconvert watertight_dagmc.h5m "$ROOT.stl"
