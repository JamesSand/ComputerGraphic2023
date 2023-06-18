#!/usr/bin/env bash
set -e

ROUNDS=5000
PHOTONS=1000000
CKPT_INTERVAL=3
METHOD=sppm
# If project not ready, generate cmake file.
if [[ ! -d build ]]; then
    mkdir -p build
    cd build
    cmake ..
    cd ..
fi

# Build project.
cd build
make -j
cd ..

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.

# testcases/scene18_dof_0.txt

# TEST_CASE_NAME=scene06_bunny_1k_vn
# TEST_CASE_NAME=scene19_sibenik
# TEST_CASE_NAME=scene20_diamond_r
TEST_CASE_NAME=scene18_dof

mkdir -p sppm_output
mkdir -p sppm_output/$TEST_CASE_NAME
time bin/PA1 testcases/$TEST_CASE_NAME.txt sppm_output/$TEST_CASE_NAME $ROUNDS $PHOTONS $CKPT_INTERVAL





# mkdir -p sppm_output/scene19_sibenik
# time bin/PA1 testcases/scene19_sibenik.txt sppm_output/scene19_sibenik $METHOD $ROUNDS $PHOTONS $CKPT_INTERVAL
# mkdir -p sppm_output/scene21_livingroom
# time bin/PA1 testcases/scene21_livingroom.txt sppm_output/scene21_livingroom $METHOD $ROUNDS $PHOTONS $CKPT_INTERVAL
# mkdir -p sppm_output/scene19_sibenik_lucy
# time bin/PA1 testcases/scene19_sibenik.txt sppm_output/scene19_sibenik_lucy $METHOD $ROUNDS $PHOTONS $CKPT_INTERVAL
# mkdir -p sppm_output/scene19_sibenik_lucy_halton
# time bin/PA1 testcases/scene19_sibenik.txt sppm_output/scene19_sibenik_lucy_halton $METHOD $ROUNDS $PHOTONS $CKPT_INTERVAL