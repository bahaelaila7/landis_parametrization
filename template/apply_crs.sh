#!/bin/bash

ref="$1"     # Reference raster (e.g., initial-communities.tif)
target="$2"  # Target raster to fix (e.g., biomass-spp-Map-001.tif)

# Extract CRS WKT from reference
crs=$(gdalsrsinfo -o wkt "$ref")

# Extract bounding box and resolution
ulx=$(gdalinfo "$ref" | grep "Upper Left" | grep -oE '[-0-9.]+[,)]' | head -n1 | tr -d ',)')
uly=$(gdalinfo "$ref" | grep "Upper Left" | grep -oE '[-0-9.]+[,)]' | tail -n1 | tr -d ',)')
lrx=$(gdalinfo "$ref" | grep "Lower Right" | grep -oE '[-0-9.]+[,)]' | head -n1 | tr -d ',)')
lry=$(gdalinfo "$ref" | grep "Lower Right" | grep -oE '[-0-9.]+[,)]' | tail -n1 | tr -d ',)')
pxsize_x=$(gdalinfo "$ref" | grep "Pixel Size" | grep -oE '[-0-9.]+' | head -n1)
pxsize_y=$(gdalinfo "$ref" | grep "Pixel Size" | grep -oE '[-0-9.]+' | tail -n1)
echo $ulx $uly
echo $lrx $lry
echo $pxsize_x $pxsize_y

# Apply spatial metadata to target
gdal_edit.py -a_nodata 0 -a_srs "$crs" -a_ullr "$ulx" "$uly" "$lrx" "$lry" "$target"
