
import argparse
from generate_mesh import generate_mesh_from_points
from generate_fiber2D_biv import generate_fiber2D
from generate_mesh import generate_mesh_from_matlab

parser = argparse.ArgumentParser()
parser.add_argument('-epi', type=str, default='epi.txt', help='File with segmentation epicardium points')
parser.add_argument('-vd', type=str, default='vd.txt', help='File with segmentation right ventricle points')
parser.add_argument('-ve', type=str, default='ve.txt', help='File with segmentation left ventricle points')
parser.add_argument('-numfib', type=int, default=0, help='Number of fibroses')
parser.add_argument('-fibbase', type=str, default='f', help='Base name for fibrosis files. Files must start with 0 id.')
parser.add_argument('-o', type=str, default='output', help='Output file name')
parser.add_argument('-t', type=int, default=0, help='Choose input type: txt or matlab')
parser.add_argument('-m', type=str, default='0', help='matlab file name')
parser.add_argument('-s', type=int, default=5, help='Slice number')
args = parser.parse_args()

if args.t == 0:
    if args.fibbase == "-1":
        generate_mesh_from_points(args.epi,args.vd,args.ve,0, 0, args.o)
    else:
        generate_mesh_from_points(args.epi,args.vd,args.ve,args.fibbase, args.numfib, args.o)
    generate_fiber2D(args.o, args.numfib)
else:

    numfib = generate_mesh_from_matlab(args.m, args.o, args.s)
    generate_fiber2D(args.o, numfib)
