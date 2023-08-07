#!/bin/bash

set -e

printf "\nStart grid generate ...\n";
time /data3/lihl/code-lihl/2d-grid-generation/example/../main_grid_2d /data3/lihl/code-lihl/2d-grid-generation/example/../project/test.json 100 2>&1 |tee log
if [ 0 -ne 0 ]; then
    printf "\ngrid generate fail! stop!\n"
    exit 1
fi

