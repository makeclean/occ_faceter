#!/bin/bash

set -euo pipefail

./faceter_test

python dagmc_canon.py cubit_output.h5m > cubit_output.txt
python dagmc_canon.py test_output.h5m > test_output.txt

exec diff -u cubit_output.txt test_output.txt
