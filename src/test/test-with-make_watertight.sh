#!/bin/bash

set -euo pipefail
cd workdir

for f in ../geometries/*-merged.brep; do
  base=`basename "$f"`
  dagfile="$base.h5m"

  echo "processing $base"
  
  occ_faceter -t 1e-3 "$f" -o "$dagfile" -f ../test_materials.txt --add_mat_ids

  make_watertight "$dagfile" > "$base.out"
done
