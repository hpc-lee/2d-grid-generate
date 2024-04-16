#!/bin/bash

set -e

printf "\nStart grid generate ...\n";
time /data/lihl/code/2d-grid-generate/new-parabolic/example/../main_grid_2d /data/lihl/code/2d-grid-generate/new-parabolic/example/../project4/test.json 100 2>&1 |tee log
if [ 0 -ne 0 ]; then
    printf "\ngrid generate fail! stop!\n"
    exit 1
fi

