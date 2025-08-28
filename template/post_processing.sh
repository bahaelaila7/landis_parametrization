#!/bin/bash
./apply_crs.sh initial-communities.tif output-community-0.tif
timestep=5
for i in $(seq 1 10); do ./apply_crs.sh initial-communities.tif biomass-succession\\biomass-anpp-$(expr $i \* $timestep).tif; done;
