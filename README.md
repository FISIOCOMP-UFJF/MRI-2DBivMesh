# MRI-2DBivMesh

Generates 2D biventricular `.alg` meshes from cardiac MRI segmentations for electrophysiology simulators.

---

## Pre-Requisites

All dependencies are included in the provided Conda environment — no manual installation needed.

- FEniCS 2019.1.0
- Gmsh
- meshio / meshio-gmsh
- scikit-image
- h5py, scipy, matplotlib
- CMake / VTK (libvtk9-dev)
- [hexa-mesh-from-VTK](https://github.com/FilipeNamorato/hexa-mesh-from-VTK_vtk9.git) — cloned automatically during configuration

---

## Installation

```sh
git clone https://github.com/FISIOCOMP-UFJF/MRI-2DBivMesh.git
cd MRI-2DBivMesh
conda env create -f mri2d.yml
conda activate mri2d
bash config.sh
```

---

## Usage

All generation is done through `generate_alg.py`, which has two subcommands.

### From a Segment `.mat` file

```sh
# All valid slices (auto-detect fibrosis)
python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2

# Single slice
python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2 --slice 4

# Custom resolution (default: 0.1 mm)
python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2 --dx 0.2 --dy 0.2 --dz 0.2
```

### From plain-text contour files

```sh
# With fibrosis
python3 generate_alg.py txt -epi epi.txt -vd vd.txt -ve ve.txt -o MyMesh --numfib 3 --fibbase f

# Without fibrosis
python3 generate_alg.py txt -epi epi.txt -vd vd.txt -ve ve.txt -o MyMesh
```

---

## Fibrosis marking (`--scar`)

| Mode | Behaviour |
|---|---|
| `auto` *(default)* | Uses `GreyZone.map` from the `.mat` if available; falls back to `Roi` otherwise |
| `greyzone` | Forces `GreyZone.map` — tag 1 = border zone, tag 2 = core fibrosis |
| `roi` | Forces `Roi`-based marking — tag 1 = fibrosis |
| `none` | No fibrosis marking |

```sh
python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2 --scar greyzone
python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2 --scar roi
python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2 --scar none
```

Material tags written to the `.alg`:

| Value | Tissue |
|---|---|
| `0` | Healthy myocardium |
| `1` | Greyzone / border zone fibrosis |
| `2` | Core fibrosis *(only in `greyzone` mode)* |

---

## Output structure

```
outputs/
  <run_name>/
    <run_name>_slice_4.alg   ← simulation input
    <run_name>_slice_4.vtu   ← visualization (ParaView)
    <run_name>_slice_5.alg
    <run_name>_slice_5.vtu
    ...
    tmp/                     ← intermediate files (.msh, .xml, .xdmf, ...)
```

---

## Help

```sh
python3 generate_alg.py --help
python3 generate_alg.py mat --help
python3 generate_alg.py txt --help
```

---

## How to cite

If you use this pipeline, please cite:

NAMORATO, F. L.; SOARES, T. J.; VIOL, R. D.; SCHMAL, T. R.; DOS SANTOS, R. W.; CAMPOS, J. O. *Geração semi-automática de geometrias computacionais biventriculares personalizadas por paciente.* In: XXIX Congresso Brasileiro de Engenharia Biomédica (CBEB 2024), Ribeirão Preto, Brasil, 2024. p. 621–627.
https://sbeb.org.br/wp-content/uploads/2024/10/SBEB-CBEB2024-ProceedingsFinal.pdf

---

Also related:

PEREIRA, J. P. B. et al. *Pipeline para avaliação do risco arrítmico com modelos computacionais personalizados baseados em ressonância magnética cardíaca e eletrocardiograma.* XXVI Encontro Nacional de Modelagem Computacional, Nova Friburgo, 2023.
https://www.even3.com.br/anais/xxvi-encontro-nacional-de-modelagem-computacional-xiv-encontro-de-ciencia-e-tecnologia-dos-materiais-338941/696802-pipeline-para-avaliacao-do-risco-arritmico-com-modelos-computacionais-personalizados-baseados-em-ressonancia-magn/
