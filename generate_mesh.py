import gmsh
import numpy as np
import argparse

parser = argparse.ArgumentParser() #arguments for each segmentation

def generate_mesh_from_points(epi, vd, ve, fibbase, numfib, outputfile):
    lc = 1

    #declaration of the argument related to the name of fibrosis as F.
    epi_points = np.loadtxt(epi)
    vd_points = np.loadtxt(vd)
    ve_points = np.loadtxt(ve)
    
    if numfib > 0:
        fib_points = [] #add all the fibrosis points into a vector.
        for i in range(numfib):
            fibfile = fibbase + str(i) + ".txt"
            try: #error checking when loading files fib_points.
                fib_points.append(np.loadtxt(fibfile))
            except Exception as e:
                print(f"Error loading points from the file {fibfile}: {e}")
    
    gmsh.initialize()

    #checking outer files
    epi = [] 
    for pt in epi_points:
        epi.append(gmsh.model.geo.addPoint(pt[0], pt[1], pt[2], lc))

    vd = []
    for pt in vd_points:
        vd.append(gmsh.model.geo.addPoint(pt[0], pt[1], pt[2], lc))

    ve = []
    for pt in ve_points:
        ve.append(gmsh.model.geo.addPoint(pt[0], pt[1], pt[2], lc))

    sp_epi = gmsh.model.geo.addSpline(epi)
    sp_epi2 = gmsh.model.geo.addSpline([epi[-1], epi[0]])
    sp_vd = gmsh.model.geo.addSpline(vd)
    sp_vd2 = gmsh.model.geo.addSpline([vd[-1], vd[0]])
    sp_ve = gmsh.model.geo.addSpline(ve)
    sp_ve2 = gmsh.model.geo.addSpline([ve[-1], ve[0]])

    cl_epi = gmsh.model.geo.addCurveLoop([sp_epi, sp_epi2])
    gmsh.model.addPhysicalGroup(1, [sp_epi, sp_epi2], tag = 30)

    cl_vd = gmsh.model.geo.addCurveLoop([sp_vd, sp_vd2])
    gmsh.model.addPhysicalGroup(1, [sp_vd, sp_vd2], tag = 40)

    cl_ve = gmsh.model.geo.addCurveLoop([sp_ve, sp_ve2])
    gmsh.model.addPhysicalGroup(1, [sp_ve, sp_ve2], tag = 20)

    cl_list = [cl_epi, cl_vd, cl_ve]


    splines = []  #list of points representing the curve.
    fib_splines = []
    cl_fibs = []
    
    if numfib > 0:
        for fib_pts in fib_points:  #coordinates of the points
            fpt = []
            for pt in fib_pts:  #add the points to gmsh
                fpt.append(gmsh.model.geo.addPoint(pt[0], pt[1], pt[2], lc))
        
            fib_sp = gmsh.model.geo.addSpline(fpt)
            fib_sp2 = gmsh.model.geo.addSpline([fpt[-1], fpt[0]])
            cl = gmsh.model.geo.addCurveLoop([fib_sp, fib_sp2])
            cl_fibs.append(cl)
            cl_list.append(cl)

    gmsh.model.geo.synchronize()

    #create surface
    surface = gmsh.model.geo.addPlaneSurface(cl_list)
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(2, [surface], tag = 0)

    fib_surf = []
    for cl in cl_fibs:
        fib_surf.append(gmsh.model.geo.addPlaneSurface([cl]))

    if numfib > 0:
        gmsh.model.addPhysicalGroup(2, fib_surf, tag = 1)    

    gmsh.model.geo.synchronize()

    gmsh.write(outputfile + '.geo_unrolled')

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2) #save in ASCII 2 format
    gmsh.write(outputfile + ".msh")
    gmsh.clear()
    gmsh.finalize()
    
    return outputfile+".msh"