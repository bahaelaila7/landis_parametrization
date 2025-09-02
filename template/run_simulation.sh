#!/bin/bash
localdir=$(pwd)
#/usr/bin/apptainer exec -C -B "$localdir" ~/data/.apptainer-root/landis-ii-v8-release.sif /bin/bash -c "cd \"$localdir\" && dotnet /bin/LANDIS_Linux/build/Release/Landis.Console.dll scenario.txt" && ./post_processing.sh
#apptainer exec -C -B "$localdir" ~/data/.apptainer-root/landis-ii-v8-release.sif /bin/bash -c "cd \"$localdir\" && dotnet \$LANDIS_CONSOLE scenario.txt" && ./post_processing.sh
apptainer exec -C -B "$localdir" ~/data/.apptainer-root/landis-ii-v8-release.sif /bin/bash -c "cd \"$localdir\" && dotnet \$LANDIS_CONSOLE scenario.txt > exp_logs"
#/usr/bin/apptainer exec -C -B "$localdir" ~/data/.apptainer-root/landis-ii-v8-release.sif /bin/bash 

