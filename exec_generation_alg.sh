#!/bin/bash
#set -x

epi=$1
vd=$2
ve=$3
numfib=$4
fibbase=$5
output_file=$6
dx=$7
dy=$8
dz=$9

dir="MRI-2DBivMesh"

if [ $(basename "$PWD") = "$dir" ]; then

    if [ $# -ne 9 ]; then
        echo "Usage: $0 epi vd ve numfib fibbase output_file dx dy dz"
        exit 1
    fi

    python3 generate_alg.py -epi "$epi" -vd "$vd" -ve "$ve" -numfib "$numfib" -fibbase "$fibbase" -o "$output_file"
    python_exit_code=$?
    if [ $python_exit_code -ne 0 ]; then
        echo "Error: Failed to execute generate_alg.py (exit code: $python_exit_code)"
        exit 1
    fi

    cd hexa-mesh-from-VTK/
    ./bin/HexaMeshFromVTK -i "../$output_file.vtu" --dx "$dx" --dy "$dy" --dz "$dz" -r 1000 -o "../$output_file.alg" -c ../config_file.ini --2d
    
else
    echo "You are not in the desired directory."
    echo "Please go to the ${dir} directory and run the script again."    
    exit 1
fi

exit 0
