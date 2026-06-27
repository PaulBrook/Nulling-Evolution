#!/bin/bash

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

python -u "${BASE_DIR}/python_code/nf_calculator.py" \
    -pulsar J1559-5545 \
    -fb 114 \
    -lb 139 \
    -outdir "${BASE_DIR}/pulsar_results/" \
    -datadir "${BASE_DIR}/pulsar_data/"
