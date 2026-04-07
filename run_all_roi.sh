#!/bin/bash
# Gera todos os pacientes do pendrive com --scar roi e resolução 100µm (dx=dy=dz=0.1)
set -e

cd /home/filipe/MRI-2DBivMesh

PATIENTS=(Patient_1 Patient_2 Patient_3 Patient_5 Patient_6 Patient_7)

for P in "${PATIENTS[@]}"; do
    echo ""
    echo "=========================================="
    echo "  Gerando $P  (ROI, dx=dy=dz=0.1)"
    echo "=========================================="
    python3 generate_alg.py mat \
        -m segmentation/${P}.mat \
        -o ${P} \
        --scar roi \
        --dx 0.1 --dy 0.1 --dz 0.1
    echo "  ✓ $P concluído → outputs/${P}/"
done

echo ""
echo "Todos os pacientes concluídos."
