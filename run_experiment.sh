#!/bin/bash
PREFIX="template-test1"
rm -rf "$PREFIX"
python generate_experiment_files.py "$PREFIX" 
pushd "$PREFIX" 
./run_simulation.sh 
popd 
python collect_run_data.py "$PREFIX"
