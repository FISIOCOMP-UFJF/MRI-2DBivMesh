#!/bin/bash
#set -x

t=$1

if [ $t -eq 0 ]; then
    epi=$2
    vd=$3
    ve=$4
    numfib=$5
    fibbase=$6
    output_file=$7
    dx=$8
    dy=$9
    dz=${10}
    m=0
    slice=0
elif [ $t -eq 1 ]; then
    output_file=$2
    m=$3
    slice=$4
    dx=$5
    dy=$6
    dz=$7
    epi=0
    vd=0
    ve=0
    numfib=0
    fibbase=0
else
    echo "Error: Invalid type segmentation."
    exit 1
fi

dir="MRI-2DBivMesh"

if [ $(basename "$PWD") != "$dir" ]; then
    echo "You are not in the desired directory."
    echo "Please go to the ${dir} directory and run the script again."
    exit 1
else
    if [ $t -eq 0 ]; then
        if [ $# -ge 8 ] && [ $# -le 11 ]; then
            if [ $# -ge 8 ] && [ $# -le 9 ]; then
                if [ $# -eq 9 ] && [ $numfib -eq 0 ]; then
                    dz=$dy
                    dy=$dx
                    dx=$output_file
                    output_file=$fibbase
                    fibbase=-1
                else
                    echo "The number of fibers needs to be zero to omit fibbase"
                    exit 1
                fi
            fi
        fi

        python3 generate_alg.py -epi "$epi" -vd "$vd" -ve "$ve" -numfib "$numfib" -fibbase "$fibbase" -o "$output_file" -t 0 -m 0
        python_exit_code=$?
        if [ $python_exit_code -ne 0 ]; then
            echo "Error: Failed to execute generate_alg.py (exit code: $python_exit_code)"
            exit 1
        fi

        echo "Warning: the script considers that the mesh is in millimeters. For other units, it is necessary to change the unit_factor (-r) in exec_generation.sh."
        cd hexa-mesh-from-VTK/
        ./bin/HexaMeshFromVTK -i "../$output_file.vtu" --dx "$dx" --dy "$dy" --dz "$dz" -r 1000 -o "../$output_file.alg" -c ../config_file.ini --2d

    else
        if [ $# -ne 7 ]; then
            echo "Error: Invalid number of parameters."
            exit 1
        fi
        python3 generate_alg.py -epi 0 -vd 0 -ve 0 -numfib 0 -fibbase 0 -o "$output_file" -t 1 -m "$m" -s "$slice"
        python_exit_code=$?
        if [ $python_exit_code -ne 0 ]; then
            echo "Error: Failed to execute generate_alg.py (exit code: $python_exit_code)"
            exit 1
        fi
        echo "Warning: the script considers that the mesh is in millimeters. For other units, it is necessary to change the unit_factor (-r) in exec_generation.sh."
        cd hexa-mesh-from-VTK/
        ./bin/HexaMeshFromVTK -i "../$output_file.vtu" --dx "$dx" --dy "$dy" --dz "$dz" -r 1000 -o "../$output_file.alg" -c ../config_file.ini --2d
    fi
fi