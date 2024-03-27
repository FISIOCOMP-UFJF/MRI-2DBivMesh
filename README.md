# MRI-2DBivMesh
Generates 2D biventricular mesh from MRI to electrophysiology simulators.

# Pre-Requisites

  - FEniCS 2019.1.0
  - Gmsh
  - meshio
  - CMake
# Configuration
  ```sh
    bash config.sh
  ```
# Running

To generate .alg do:
```sh
conda activate fenicsproject
```
```sh
bash exec_generation_alg.sh epi vd ve numfib fibbase output_file_name dx dy dz
```

# Running example
```sh
bash exec_generation_alg.sh ./segmentation/epi9.txt ./segmentation/endoVD9.txt ./segmentation/endoVE9.txt 3 ./segmentation/fibr9_ output_file 0.2 0.2 0.2
```
# Dependencies

This project depends on the following repositories:

- [hexa-mesh-from-VTK](https://github.com/rsachetto/hexa-mesh-from-VTK.git): This repository is necessary for the generation of hexahedral meshes from VTK files. It will be cloned during the Configuration.

# How to cite:
----

PEREIRA, J. P. B. ; SOARES, T. J. ; WERNECK, Y. B. ; ALMEIDA, D. K. ; SANTOS, Y. R. A. ; SANTOS, F. J. M. ; FRANCO, T. D. ; OLIVEIRA, R. S. ; SCHMAL, T. R. ; SOUZA, T. G. S. E. ; ROCHA, B. M. ; CAMPOS, J. O. ; DOS SANTOS, R. W. . PIPELINE PARA AVALIAÇÃO DO RISCO ARRÍTMICO COM MODELOS COMPUTACIONAIS PERSONALIZADOS BASEADOS EM RESSONÂNCIA MAGNÉTICA CARDÍACA E ELETROCARDIOGRAMA. In: XXVI Encontro Nacional de Modelagem Computacional e XIV Encontro de Ciência e Tecnologia dos Materiais, 2023, Nova Friburgo. Anais do XXVI Encontro Nacional de Modelagem Computacional e XIV Encontro de Ciência e Tecnologia dos Materiais, 2023.
https://www.even3.com.br/anais/xxvi-encontro-nacional-de-modelagem-computacional-xiv-encontro-de-ciencia-e-tecnologia-dos-materiais-338941/696802-pipeline-para-avaliacao-do-risco-arritmico-com-modelos-computacionais-personalizados-baseados-em-ressonancia-magn/