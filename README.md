# MRI-2DBivMesh
Generates 2D biventricular mesh from MRI to electrophysiology simulators.

# Pre-Requisites

  - FEniCS 2019.1.0
  - Gmsh

# Running examples
----
To generate .alg do:

```sh
$ python generate_alg.py -epi ./segmentation/epi9.txt -vd ./segmentation/endoVD9.txt -ve ./segmentation/endoVE9.txt -numfib 3 -fibbase ./segmentation/fibr9_ -o demo_biv_mesh
```
# How to cite:
----

PEREIRA, J. P. B. ; SOARES, T. J. ; WERNECK, Y. B. ; ALMEIDA, D. K. ; SANTOS, Y. R. A. ; SANTOS, F. J. M. ; FRANCO, T. D. ; OLIVEIRA, R. S. ; SCHMAL, T. R. ; SOUZA, T. G. S. E. ; ROCHA, B. M. ; CAMPOS, J. O. ; DOS SANTOS, R. W. . PIPELINE PARA AVALIAÇÃO DO RISCO ARRÍTMICO COM MODELOS COMPUTACIONAIS PERSONALIZADOS BASEADOS EM RESSONÂNCIA MAGNÉTICA CARDÍACA E ELETROCARDIOGRAMA. In: XXVI Encontro Nacional de Modelagem Computacional e XIV Encontro de Ciência e Tecnologia dos Materiais, 2023, Nova Friburgo. Anais do XXVI Encontro Nacional de Modelagem Computacional e XIV Encontro de Ciência e Tecnologia dos Materiais, 2023.
https://www.even3.com.br/anais/xxvi-encontro-nacional-de-modelagem-computacional-xiv-encontro-de-ciencia-e-tecnologia-dos-materiais-338941/696802-pipeline-para-avaliacao-do-risco-arritmico-com-modelos-computacionais-personalizados-baseados-em-ressonancia-magn/