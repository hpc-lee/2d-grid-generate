#!/bin/bash

set -e

printf "\nStart grid post process ...\n";
time /data/lihl/code/2d-grid-generate/grid-post-process/example/../grid_post_proc_2d /data/lihl/code/2d-grid-generate/grid-post-process/example/../project/test.json 100 2>&1 |tee log
if [ 0 -ne 0 ]; then
    printf "\ngrid generate fail! stop!\n"
    exit 1
fi

