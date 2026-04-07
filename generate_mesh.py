import os
import gmsh
import numpy as np
import scipy.io
import meshio
from matplotlib.path import Path as MplPath
from skimage.measure import find_contours


def _mark_fibrosis_post(msh_path, polygons_by_tag):
    """
    Marks triangles by point-in-polygon and rewrites the .msh.

    polygons_by_tag: dict mapping tag (int) -> list of Nx2 polygon arrays.
      Tag 1 = greyzone / ROI fibrosis
      Tag 2 = core fibrosis (overrides tag 1)
    Tags are applied in ascending order so higher tags win.
    """
    mesh = meshio.read(msh_path)

    tri_idx = next((i for i, cb in enumerate(mesh.cells) if cb.type == "triangle"), None)
    if tri_idx is None:
        return

    tri = mesh.cells[tri_idx].data
    centroids = mesh.points[tri].mean(axis=1)[:, :2]
    material = np.zeros(len(tri), dtype=np.int32)

    for tag in sorted(polygons_by_tag.keys()):
        for poly in polygons_by_tag[tag]:
            inside = MplPath(poly[:, :2]).contains_points(centroids)
            material[inside] = tag

    new_cell_data = {
        key: list(arrays) for key, arrays in mesh.cell_data.items()
    }
    if new_cell_data:
        for key in new_cell_data:
            new_cell_data[key][tri_idx] = material
    else:
        new_cell_data = {"gmsh:physical": [None] * len(mesh.cells)}
        new_cell_data["gmsh:physical"][tri_idx] = material

    meshio.write(
        msh_path,
        meshio.Mesh(points=mesh.points, cells=mesh.cells, cell_data=new_cell_data),
        file_format="gmsh22",
        binary=False,
    )


def _read_greyzone_polygons(data_simple, slice_idx):
    """
    Extracts contour polygons from Scar.GreyZone.map.

    data_simple: .mat loaded with simplify_cells=True
    Returns (gz_polygons, core_polygons) or ([], []) if data not available.
    Coordinate convention: find_contours dim0 = X, dim1 = Y (matches mesh coords).
    """
    ss = data_simple['setstruct']

    try:
        gz_map = np.array(ss['Scar']['GreyZone']['map'])
    except (KeyError, TypeError, ValueError):
        return [], []

    if gz_map.ndim != 3 or slice_idx >= gz_map.shape[2]:
        return [], []

    plane = gz_map[:, :, slice_idx]

    def contours_for_value(val):
        mask = (plane == val).astype(float)
        polys = []
        for c in find_contours(mask, level=0.5):
            polys.append(np.column_stack([c[:, 0], c[:, 1]]))
        return polys

    return contours_for_value(1), contours_for_value(2)


def _load_roi_polygons(data, slice_idx):
    """Returns (polygons_by_tag, numfib) from setstruct.Roi for the given slice."""
    dataScar    = data['setstruct']['Roi'][0][0][0]
    numfibTotal = data['setstruct']['RoiN'][0][0][0][0]
    roi_polys   = []
    for i in range(numfibTotal):
        if dataScar[i][3][0] - 1 == slice_idx:  # ROI Z is 1-indexed in Matlab
            scar = np.zeros((40, 3))
            scar[:, 0] = dataScar[i][0][0:80:2, 0]
            scar[:, 1] = dataScar[i][1][0:80:2, 0]
            roi_polys.append(scar)
    polygons_by_tag = {1: roi_polys} if roi_polys else {}
    return polygons_by_tag, len(roi_polys)


def generate_mesh_from_matlab(patient, outputfile, slice_idx,
                               tmp_dir="./outputs_other", scar_mode="auto", lc=1):
    """
    scar_mode: 'auto' | 'greyzone' | 'roi' | 'none'
      auto     – use GreyZone.map if available, else ROI
      greyzone – force GreyZone (tag1=border zone, tag2=core)
      roi      – force ROI-based marking (tag1=fibrosis)
      none     – skip fibrosis entirely
    """
    os.makedirs(tmp_dir, exist_ok=True)
    data        = scipy.io.loadmat(patient)
    data_simple = scipy.io.loadmat(patient, simplify_cells=True)

    endoX = data['setstruct']['EndoX'][0][0]
    xlv = endoX[:, 0, slice_idx]
    endoY = data['setstruct']['EndoY'][0][0]
    ylv = endoY[:, 0, slice_idx]

    RVEndoX = data['setstruct']['RVEndoX'][0][0]
    xrv = RVEndoX[:, 0, slice_idx]
    RVEndoY = data['setstruct']['RVEndoY'][0][0]
    yrv = RVEndoY[:, 0, slice_idx]

    RVEpiX = data['setstruct']['RVEpiX'][0][0]
    xepi = RVEpiX[:, 0, slice_idx]
    RVEpiY = data['setstruct']['RVEpiY'][0][0]
    yepi = RVEpiY[:, 0, slice_idx]

    n = len(xlv)
    sz = len(xlv[0:n:2])

    epi_points = np.zeros((sz, 3))
    vd_points  = np.zeros((sz, 3))
    ve_points  = np.zeros((sz, 3))

    ve_points[:, 0] = xlv[0:n:2]
    ve_points[:, 1] = ylv[0:n:2]
    vd_points[:, 0] = xrv[0:n:2]
    vd_points[:, 1] = yrv[0:n:2]
    epi_points[:, 0] = xepi[0:n:2]
    epi_points[:, 1] = yepi[0:n:2]

    # --- Resolve fibrosis polygons according to scar_mode ---
    polygons_by_tag = {}
    numfib = 0

    if scar_mode == "none":
        print("  [fibrose] skipped (--scar none)")

    elif scar_mode in ("greyzone", "auto"):
        gz_polys, core_polys = _read_greyzone_polygons(data_simple, slice_idx)
        has_gz = bool(gz_polys or core_polys)

        if has_gz:
            if gz_polys:
                polygons_by_tag[1] = gz_polys
            if core_polys:
                polygons_by_tag[2] = core_polys
            numfib = len(gz_polys) + len(core_polys)
            print(f"  [fibrose] GreyZone: {len(gz_polys)} border-zone, {len(core_polys)} core")
        elif scar_mode == "greyzone":
            print("  [fibrose] WARNING: --scar greyzone requested but no GreyZone.map found")
        else:
            # auto fallback to ROI
            polygons_by_tag, numfib = _load_roi_polygons(data, slice_idx)
            print(f"  [fibrose] ROI fallback: {numfib} region(s)")

    elif scar_mode == "roi":
        polygons_by_tag, numfib = _load_roi_polygons(data, slice_idx)
        print(f"  [fibrose] ROI: {numfib} region(s)")

    # --- Gmsh: mesh cardiac wall only (no fibrosis holes) ---
    gmsh.initialize()
    lc = 0.5
    epi_pts = [gmsh.model.geo.addPoint(p[0], p[1], p[2], lc) for p in epi_points]
    vd_pts  = [gmsh.model.geo.addPoint(p[0], p[1], p[2], lc) for p in vd_points]
    ve_pts  = [gmsh.model.geo.addPoint(p[0], p[1], p[2], lc) for p in ve_points]

    sp_epi  = gmsh.model.geo.addSpline(epi_pts)
    sp_epi2 = gmsh.model.geo.addSpline([epi_pts[-1], epi_pts[0]])
    sp_vd   = gmsh.model.geo.addSpline(vd_pts)
    sp_vd2  = gmsh.model.geo.addSpline([vd_pts[-1], vd_pts[0]])
    sp_ve   = gmsh.model.geo.addSpline(ve_pts)
    sp_ve2  = gmsh.model.geo.addSpline([ve_pts[-1], ve_pts[0]])

    cl_epi = gmsh.model.geo.addCurveLoop([sp_epi, sp_epi2])
    gmsh.model.addPhysicalGroup(1, [sp_epi, sp_epi2], tag=30)
    cl_vd  = gmsh.model.geo.addCurveLoop([sp_vd, sp_vd2])
    gmsh.model.addPhysicalGroup(1, [sp_vd, sp_vd2], tag=40)
    cl_ve  = gmsh.model.geo.addCurveLoop([sp_ve, sp_ve2])
    gmsh.model.addPhysicalGroup(1, [sp_ve, sp_ve2], tag=20)

    gmsh.model.geo.synchronize()
    surface = gmsh.model.geo.addPlaneSurface([cl_epi, cl_vd, cl_ve])
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(2, [surface], tag=0)
    gmsh.model.geo.synchronize()

    gmsh.write(os.path.join(tmp_dir, outputfile + ".geo_unrolled"))
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    msh_path = os.path.join(tmp_dir, outputfile + ".msh")
    gmsh.write(msh_path)
    gmsh.clear()
    gmsh.finalize()

    # --- Post-process: mark fibrosis by point-in-polygon ---
    if polygons_by_tag:
        _mark_fibrosis_post(msh_path, polygons_by_tag)

    return numfib


def generate_mesh_from_points(epi, vd, ve, fibbase, numfib, outputfile,
                               tmp_dir="./outputs_other"):
    os.makedirs(tmp_dir, exist_ok=True)
    lc = 1

    epi_points = np.loadtxt(epi)
    vd_points  = np.loadtxt(vd)
    ve_points  = np.loadtxt(ve)

    if numfib > 0:
        fib_points = []
        for i in range(numfib):
            fibfile = fibbase + str(i) + ".txt"
            try:
                fib_points.append(np.loadtxt(fibfile))
            except Exception as e:
                print(f"Error loading points from the file {fibfile}: {e}")

    gmsh.initialize()

    epi_pts = [gmsh.model.geo.addPoint(p[0], p[1], p[2], lc) for p in epi_points]
    vd_pts  = [gmsh.model.geo.addPoint(p[0], p[1], p[2], lc) for p in vd_points]
    ve_pts  = [gmsh.model.geo.addPoint(p[0], p[1], p[2], lc) for p in ve_points]

    sp_epi  = gmsh.model.geo.addSpline(epi_pts)
    sp_epi2 = gmsh.model.geo.addSpline([epi_pts[-1], epi_pts[0]])
    sp_vd   = gmsh.model.geo.addSpline(vd_pts)
    sp_vd2  = gmsh.model.geo.addSpline([vd_pts[-1], vd_pts[0]])
    sp_ve   = gmsh.model.geo.addSpline(ve_pts)
    sp_ve2  = gmsh.model.geo.addSpline([ve_pts[-1], ve_pts[0]])

    cl_epi = gmsh.model.geo.addCurveLoop([sp_epi, sp_epi2])
    gmsh.model.addPhysicalGroup(1, [sp_epi, sp_epi2], tag=30)
    cl_vd  = gmsh.model.geo.addCurveLoop([sp_vd, sp_vd2])
    gmsh.model.addPhysicalGroup(1, [sp_vd, sp_vd2], tag=40)
    cl_ve  = gmsh.model.geo.addCurveLoop([sp_ve, sp_ve2])
    gmsh.model.addPhysicalGroup(1, [sp_ve, sp_ve2], tag=20)

    cl_list = [cl_epi, cl_vd, cl_ve]

    if numfib > 0:
        for fib_pts in fib_points:
            fpt = [gmsh.model.geo.addPoint(p[0], p[1], p[2], lc) for p in fib_pts]
            fib_sp  = gmsh.model.geo.addSpline(fpt)
            fib_sp2 = gmsh.model.geo.addSpline([fpt[-1], fpt[0]])
            cl = gmsh.model.geo.addCurveLoop([fib_sp, fib_sp2])
            cl_list.append(cl)

    gmsh.model.geo.synchronize()
    surface = gmsh.model.geo.addPlaneSurface(cl_list)
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(2, [surface], tag=0)

    if numfib > 0:
        fib_surf = []
        for cl in cl_list[3:]:
            fib_surf.append(gmsh.model.geo.addPlaneSurface([cl]))
        gmsh.model.addPhysicalGroup(2, fib_surf, tag=1)

    gmsh.model.geo.synchronize()
    gmsh.write(os.path.join(tmp_dir, outputfile + ".geo_unrolled"))
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.write(os.path.join(tmp_dir, outputfile + ".msh"))
    gmsh.clear()
    gmsh.finalize()

    return outputfile + ".msh"


def verify_zSize(patient):
    data = scipy.io.loadmat(patient)
    return data['setstruct']['ZSize'][0][0][0][0]


def verify_slices(patient):
    data = scipy.io.loadmat(patient)
    zSize = data['setstruct']['ZSize'][0][0][0][0]
    slices = []
    for i in range(zSize):
        endoX  = data['setstruct']['EndoX'][0][0]
        RVEndoX = data['setstruct']['RVEndoX'][0][0]
        RVEpiX  = data['setstruct']['RVEpiX'][0][0]
        xlv  = endoX[:, 0, i]
        xrv  = RVEndoX[:, 0, i]
        xepi = RVEpiX[:, 0, i]
        n = len(xlv)
        ve_pts = np.zeros((len(xlv[0:n:2]), 3))
        vd_pts = np.zeros_like(ve_pts)
        ep_pts = np.zeros_like(ve_pts)
        ve_pts[:, 0] = xlv[0:n:2]
        vd_pts[:, 0] = xrv[0:n:2]
        ep_pts[:, 0] = xepi[0:n:2]
        if not np.isnan(ve_pts[0, 0]) and not np.isnan(vd_pts[0, 0]) and not np.isnan(ep_pts[0, 0]):
            slices.append(i)
    return slices
