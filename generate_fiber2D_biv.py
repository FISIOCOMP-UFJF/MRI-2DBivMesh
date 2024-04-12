# # Generating fibers for patient specific 2D geometries

import dolfin as df
from math import pi, cos, sin
import meshio
import numpy as np
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


def generate_fiber2D(mesh_name, numfib):
    
    #convert mesh to fenics format
    meshname = mesh_name
    ifilename = meshname + '.msh'
    ofilename = meshname + '.xml'
    iformat = 'gmsh'
    meshconvert.convert2xml(ifilename, ofilename, iformat=iformat)


    # Create mesh and define function space
    mesh = df.Mesh(meshname + '.xml')
    if numfib==0:
        materials = df.MeshFunction("size_t", mesh, 2)    
        materials.set_all(0)
    else:
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
    sigma_n = df.project(sn, Vg, solver_type="gmres", preconditioner_type="hypre_amg")
    u = df.project(u, V, solver_type="gmres", preconditioner_type="hypre_amg")

    sigma_n.rename("sigma_n","fiber")
    u.rename("u","u")

    # Assigns label (healthy or fibrosis) for each element
    V0 = df.FunctionSpace(mesh, 'DG', 0)
    mat  = df.Function(V0)
    mat_values = np.array([0, 1])  # conductivity value for each region

    if numfib > 0:
        help = np.asarray(materials.array(), dtype=np.int32)
    else:
        help = np.zeros(mesh.num_cells(), dtype=np.int32)

    mat.vector()[:] = np.choose(help, mat_values) 
    mat.rename("material","material")
    # Write mesh, fiber and material label in XDMF file
    with df.XDMFFile(mesh.mpi_comm(), meshname + ".xdmf") as xdmf:
        xdmf.parameters.update(
        {
            "functions_share_mesh": True,
            "rewrite_function_mesh": False
        })
        xdmf.write(mesh)
        xdmf.write(u, 0)
        xdmf.write(sigma_n, 0)
        xdmf.write(mat, 0)
    
    convert_xdmf_to_vtu(meshname)

def convert_xdmf_to_vtu(meshname):
    filename = meshname+".xdmf"

    t, point_data, cell_data = None, None, None

    with meshio.xdmf.TimeSeriesReader(filename) as reader:
        points, cells = reader.read_points_cells()
        t, point_data, cell_data = reader.read_data(0)

    mesh = meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data,)
    mesh.write(meshname+".vtu", file_format="vtu",  )
        