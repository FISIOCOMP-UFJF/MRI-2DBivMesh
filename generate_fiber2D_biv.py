# # Generating fibers for patient specific geometries
#
# In this demo we will show how to generate fiber orientations from a patient specific geometry. We will use a mesh of an LV that is constructed using gmsh (https://gmsh.info), see https://github.com/finsberg/ldrb/blob/main/demos/mesh.msh
#
# It is important that the mesh contains physical surfaces of the endocardium (lv and rv if present), the base and the epicardium. You can find an example of how to generate such a geometry using the python API for gmsh here: https://github.com/finsberg/pulse/blob/0d7b5995f62f41df4eec9f5df761fa03da725f69/pulse/geometries.py#L160
#
# First we import the necessary packages. Note that we also import `meshio` which is used for converted from `.msh` (gmsh) to `.xdmf` (FEnICS).
import dolfin as df
from math import pi, cos, sin
import quaternion
import numpy as np
import argparse
from dolfin_utils.meshconvert import meshconvert

def solve_laplace(mesh, boundary_markers, boundary_values):
    V = df.FunctionSpace(mesh, 'P', 1)

    u_rv, u_lv, u_epi = boundary_values

    bc1 = df.DirichletBC(V, u_rv, boundary_markers, 40) 
    bc2 = df.DirichletBC(V, u_lv, boundary_markers, 20)
    bc3 = df.DirichletBC(V, u_epi, boundary_markers, 30)

    bcs=[bc1, bc2 ,bc3]

    ds = df.Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    dx = df.Measure('dx', domain=mesh)

    # Define variational problem
    u = df.TrialFunction(V)
    v = df.TestFunction(V)
    f = df.Constant(0.0)   
    a = df.dot(df.grad(u), df.grad(v))*dx  
    L = f*v*dx

    # Compute solution
    u = df.Function(V)
    df.solve(a == L, u, bcs, solver_parameters=dict(linear_solver='gmres', preconditioner='hypre_amg')) 

    return u


parser = argparse.ArgumentParser() 
parser.add_argument('-meshname', type=str, default='patient', help='Gmsh file name without extension (.msh)')
argumentos = parser.parse_args()

#convert mesh to fenics format
meshname = argumentos.meshname
ifilename = meshname + '.msh'
ofilename = meshname + '.xml'
iformat = 'gmsh'
meshconvert.convert2xml(ifilename, ofilename, iformat=iformat)


# Create mesh and define function space
mesh = df.Mesh(meshname + '.xml')
materials = df.MeshFunction("size_t", mesh, meshname + '_physical_region.xml')
boundary_markers = df.MeshFunction("size_t", mesh, meshname + '_facet_region.xml')

V = df.FunctionSpace(mesh, 'Lagrange', 1)
Vg = df.VectorFunctionSpace(mesh, 'Lagrange', 1, 3)

# Solve Laplace problems with different boundary conditions
# u=1 on epicardium
phi_epi = solve_laplace(mesh, boundary_markers, [0, 0, 1])
# u=1 on LV endocardium
phi_lv = solve_laplace(mesh, boundary_markers, [0, 1, 0])
# u=1 on RV endocardium
phi_rv = solve_laplace(mesh, boundary_markers, [1, 0, 0])

# Compute field with Laplace solutions
u = -(phi_epi + 2*phi_rv*phi_lv/(phi_rv+phi_lv) ) + 1

gU = df.grad(u)
gradU = df.as_vector([gU[0], gU[1], 0])

# Define rotation angle varying from -60 to +60
theta = df.project((-pi/3 + u*(2*pi/3)), V)

# Define Rodrigues rotation matrix
W = df.as_matrix([[0, -gradU[2], gradU[1]], [gradU[2], 0, -gradU[0]], [-gradU[1], gradU[0], 0]])
I = df.as_matrix([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
#R = I + df.sin(theta)*W + 2*(df.sin(theta/2)**2)*W*W

# Define vector orthogonal to gradient
sl = df.as_vector([-gradU[1], gradU[0], 0])
sl = sl/df.sqrt(df.dot(sl, sl))

# Rotates gradient by theta using Rodrigues rotation matrix
sn = df.cos(theta)*sl + df.sin(theta)*W*sl
sn = sn/df.sqrt(df.dot(sn, sn))

# Compute vectors for each element
#sigma_l = df.project(sl, Vg, solver_type="gmres", preconditioner_type="hypre_amg")
sigma_n = df.project(sn, Vg, solver_type="gmres", preconditioner_type="hypre_amg")
#grad_u = df.project(gradU, Vg, solver_type="gmres", preconditioner_type="hypre_amg")
u = df.project(u, V, solver_type="gmres", preconditioner_type="hypre_amg")


#phi_epi.rename("phi_epi", "phi_epi")
#phi_lv.rename("phi_lv","phi_lv")
#phi_rv.rename("phi_rv","phi_rv")
sigma_n.rename("sigma_n","fiber")
#sigma_l.rename("sigma_l","sigma_l")
u.rename("u","u")
#grad_u.rename("grad_u", "grad_u")

V0 = df.FunctionSpace(mesh, 'DG', 0)
mat  = df.Function(V0)
mat_values = np.array([0, 1])  # valor de condutividade para cada regiao

help = np.asarray(materials.array(), dtype=np.int32)
mat.vector()[:] = np.choose(help, mat_values) #atribuir o valor de condutividade para cada elemento 
mat.rename("material","material")

with df.XDMFFile(mesh.mpi_comm(), meshname + ".xdmf") as xdmf:
    xdmf.parameters.update(
    {
        "functions_share_mesh": True,
        "rewrite_function_mesh": False
    })
    xdmf.write(mesh)
#    xdmf.write(phi_epi, 0)
#    xdmf.write(phi_lv, 0)
#    xdmf.write(phi_rv, 0)
    xdmf.write(u, 0)
    xdmf.write(sigma_n, 0)
#    xdmf.write(sigma_l, 0)
#    xdmf.write(grad_u, 0)
    xdmf.write(mat, 0)
