#!/bin/bash
PREFIX=$1
SEED=$2
rm -rf "$PREFIX"
python generate_experiment_files.py "$PREFIX" $2
pushd "$PREFIX" 
./run_simulation.sh 
popd 
python collect_run_data.py "$PREFIX"
