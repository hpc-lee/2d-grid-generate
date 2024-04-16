#!/bin/bash

set -e

printf "\nUse 1 CPUs on following nodes:\n"

printf "\nStart grid generate ...\n";
time /data/apps/openmpi/4.1.5-cuda-aware/bin/mpiexec -np 1 /data/lihl/code/2d-grid-generate/elliptic/example/../main_grid_2d /data/lihl/code/2d-grid-generate/elliptic/example/../project/test.json 100 2 2>&1 |tee log
if [ 0 -ne 0 ]; then
    printf "\ngrid generate fail! stop!\n"
    exit 1
fi

