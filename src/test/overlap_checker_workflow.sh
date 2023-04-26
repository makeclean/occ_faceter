#!/bin/bash

set -euo pipefail

if [ ! -f "${1-}" ]; then
    echo "Please pass the name of a STEP file first" 1>&2
    exit 1
fi

PATH=.:$PATH

source="$1"
base="${source##*/}"
base="${base%.*}"
brep="$base.brep"
geometry="$base-geometry.csv"
overlaps="$base-overlaps.csv"
common="$base-common.brep"
imprinted="$base-imprinted.brep"
merged="$base-merged.brep"

echo "linearising solids into $brep" 1>&2

step_to_brep "$source" "$brep" > "$geometry"

echo "checking for intersecting solids" 1>&2

overlap_checker -j1 "$brep" > "$overlaps"

if grep -q overlap "$overlaps"; then
    echo "writing overlapping solds into $common" 1>&2
    grep overlap "$overlaps" | overlap_collecter "$brep" "$common"
fi

echo "removing overlaps and writing to $imprinted" 1>&2
imprint_solids "$brep" "$imprinted" < "$overlaps"

echo "merging faces, edges and verticies and writing to $merged" 1>&2
merge_solids "$imprinted" "$merged"
