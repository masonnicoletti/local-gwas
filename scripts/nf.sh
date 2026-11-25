#!/usr/bin/env bash

set -euo pipefail
mkdir -p logs
nextflow -log logs/nextflow.log run "$@"
