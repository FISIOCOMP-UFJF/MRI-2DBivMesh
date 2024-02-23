
import argparse
import generate_fiber2D_biv
import numpy as np
import gmsh
from generate_mesh import generate_mesh_from_points
from generate_fiber2D_biv import generate_fiber2D
parser = argparse.ArgumentParser()
parser.add_argument('-epi', type=str, default='epi.txt', help='File with segmentation epicardium points')
parser.add_argument('-vd', type=str, default='vd.txt', help='File with segmentation right ventricle points')
parser.add_argument('-ve', type=str, default='ve.txt', help='File with segmentation left ventricle points')

parser.add_argument('-numfib', type=int, default=0, help='Number of fibroses')

parser.add_argument('-fibbase', type=str, default='f', help='Base name for fibrosis files. Files must start with 0 id.')

parser.add_argument('-o', type=str, default='output', help='Output file name')
args = parser.parse_args()

epi_points = np.loadtxt(args.epi)
vd_points = np.loadtxt(args.vd)
ve_points = np.loadtxt(args.ve)

num_fib = args.numfib
outpuf_file = args.o
fib_points = [] #adiciona todos os pontos de fibroses em um vetor
for i in range(num_fib):
    fibfile = args.fibbase + str(i) + ".txt"
    try: #verificação de erro ao carregar os arquivos
        fib_points.append(np.loadtxt(fibfile))
    except Exception as e:
        print(f"Erro ao carregar pontos do arquivo {fibfile}: {e}")
       
generate_mesh_from_points(epi_points,vd_points,ve_points,fib_points, outpuf_file)
generate_fiber2D(outpuf_file)

