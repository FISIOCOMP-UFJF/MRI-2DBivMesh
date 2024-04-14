# MRI-2DBivMesh
Generates 2D biventricular mesh from MRI to electrophysiology simulators.

# Pre-Requisites

  - FEniCS 2019.1.0
  - Gmsh
  - meshio
  - h5py 
  - Scipy
  - CMake
  - VTK
  - [hexa-mesh-from-VTK](https://github.com/rsachetto/hexa-mesh-from-VTK.git): This repository is necessary for the generation of hexahedral meshes from VTK files. It will be cloned during the Configuration.
  
# Configuration
  ```sh
    bash config.sh
  ```

# Description parameters
  - t: Input type: txt or matlab

  - epi: epicardium segmentation

  - vd: right ventricle endocardium segmentation

  - ve: left ventricle endocardium segmentation

  - numfib: number of fibrosis files

  - fibbase: prefix of the filenames with the fibrosis segmentation
  
  - output_file_name: output file name

  - m: .mat file

  - slice: desired slice

  - dx, dy, and dz: refer to the discretization for the .vtu. Conventionally, we use the value of 0.2.
# Running


```sh
conda activate fenicsproject
```
To generate .alg do using .txt:
```sh
bash exec_generation_alg.sh 0 epi vd ve numfib fibbase output_file_name dx dy dz
```
To generate .alg do using .mat:
```sh
bash exec_generation_alg.sh 1 output_file_name patient_mat slice dx dy dz
```

# Running example
```sh
bash exec_generation_alg.sh 0 ./segmentation/epi9.txt ./segmentation/endoVD9.txt ./segmentation/endoVE9.txt 3 ./segmentation/fibr9_ output_file 0.2 0.2 0.2
```
For segmentation without fibrosis, set numFib to zero and optionally skip the fibrosis segmentation:
```sh
bash exec_generation_alg.sh 0 ./segmentation/epi9.txt ./segmentation/endoVD9.txt ./segmentation/endoVE9.txt 0 output_file 0.2 0.2 0.2
```
For using .mat:
```sh
bash exec_generation_alg.sh 1 outputfile ./segmentation/Patient_3.mat 5 0.2 0.2 0.2
```

# How to cite:

PEREIRA, J. P. B. ; SOARES, T. J. ; WERNECK, Y. B. ; ALMEIDA, D. K. ; SANTOS, Y. R. A. ; SANTOS, F. J. M. ; FRANCO, T. D. ; OLIVEIRA, R. S. ; SCHMAL, T. R. ; SOUZA, T. G. S. E. ; ROCHA, B. M. ; CAMPOS, J. O. ; DOS SANTOS, R. W. . PIPELINE PARA AVALIAÇÃO DO RISCO ARRÍTMICO COM MODELOS COMPUTACIONAIS PERSONALIZADOS BASEADOS EM RESSONÂNCIA MAGNÉTICA CARDÍACA E ELETROCARDIOGRAMA. In: XXVI Encontro Nacional de Modelagem Computacional e XIV Encontro de Ciência e Tecnologia dos Materiais, 2023, Nova Friburgo. Anais do XXVI Encontro Nacional de Modelagem Computacional e XIV Encontro de Ciência e Tecnologia dos Materiais, 2023.
https://www.even3.com.br/anais/xxvi-encontro-nacional-de-modelagem-computacional-xiv-encontro-de-ciencia-e-tecnologia-dos-materiais-338941/696802-pipeline-para-avaliacao-do-risco-arritmico-com-modelos-computacionais-personalizados-baseados-em-ressonancia-magn/