"""
generate_alg.py  –  Entry point for 2D biventricular mesh + ALG generation.

Subcommands
-----------
  mat   Read contours from a Segment .mat file (one slice or all valid slices)
  txt   Read contours from plain-text files

Fibrosis modes  (--scar)
------------------------
  auto      Use GreyZone.map if present in the .mat, otherwise fall back to ROI  [default]
  greyzone  Force GreyZone.map  (tag 1 = border zone, tag 2 = core)
  roi       Force ROI-based marking  (tag 1 = fibrosis)
  none      No fibrosis marking

Examples
--------
  # All slices, auto fibrosis detection
  python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2

  # Single slice, force GreyZone
  python3 generate_alg.py mat -m segmentation/Patient_2.mat -o Patient_2 --slice 4 --scar greyzone

  # Text files, no fibrosis
  python3 generate_alg.py txt -epi epi.txt -vd vd.txt -ve ve.txt -o MyMesh --scar none

  # Custom resolution
  python3 generate_alg.py mat -m segmentation/Patient_3.mat -o Patient_3 --dx 0.2 --dy 0.2 --dz 0.2
"""

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from generate_mesh import generate_mesh_from_points, generate_mesh_from_matlab, verify_slices
from generate_fiber2D_biv import generate_fiber2D

BASE_DIR = Path(__file__).resolve().parent
HEX_BIN  = BASE_DIR / "hexa-mesh-from-VTK_vtk9" / "bin" / "HexaMeshFromVTK"
CFG_FILE = BASE_DIR / "config_file.ini"

SCAR_CHOICES = ("auto", "greyzone", "roi", "none")


# ── Shared resolution arguments ────────────────────────────────────────────────
def add_resolution_args(p):
    p.add_argument("--dx", type=float, default=0.1)
    p.add_argument("--dy", type=float, default=0.1)
    p.add_argument("--dz", type=float, default=0.1)
    p.add_argument("--lc", type=float, default=1.0,
                   help="Gmsh characteristic length (triangle size in .msh). Default: 1.0")


def add_scar_arg(p):
    p.add_argument(
        "--scar",
        choices=SCAR_CHOICES,
        default="auto",
        help=(
            "Fibrosis marking mode: "
            "auto=GreyZone if available else ROI | "
            "greyzone=force GreyZone (tag1=border, tag2=core) | "
            "roi=ROI only (tag1=fibrosis) | "
            "none=skip fibrosis"
        ),
    )


# ── HexaMeshFromVTK call ───────────────────────────────────────────────────────
def run_hexa(mesh_name, out_dir, tmp_dir, dx, dy, dz):
    vtu = Path(tmp_dir) / f"{mesh_name}.vtu"
    alg = Path(out_dir) / f"{mesh_name}.alg"
    cmd = [
        str(HEX_BIN), "-i", str(vtu),
        "--dx", str(dx), "--dy", str(dy), "--dz", str(dz),
        "-r", "1000", "-o", str(alg), "-c", str(CFG_FILE), "--2d",
    ]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"HexaMeshFromVTK failed (exit {result.returncode}) for {mesh_name}")

    # Copy .vtu alongside the .alg for easy visualization (always overwrite)
    if vtu.exists():
        shutil.copy(vtu, Path(out_dir) / f"{mesh_name}.vtu")


# ── Subcommand: mat ────────────────────────────────────────────────────────────
def cmd_mat(args):
    out_dir = BASE_DIR / "outputs" / args.o
    tmp_dir = out_dir / "tmp"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    slices = [args.slice] if args.slice is not None else verify_slices(args.m)
    print(f"Slices to process: {slices}")

    for sl in slices:
        name = f"{args.o}_slice_{sl}"
        print(f"\n{'='*55}\nSlice {sl}  →  {name}\n{'='*55}")

        numfib = generate_mesh_from_matlab(
            args.m, name, sl,
            tmp_dir=str(tmp_dir),
            scar_mode=args.scar,
            lc=args.lc,
        )
        generate_fiber2D(name, numfib, tmp_dir=str(tmp_dir))
        run_hexa(name, out_dir, tmp_dir, args.dx, args.dy, args.dz)

    print(f"\nDone.  Results in: outputs/{args.o}/")


# ── Subcommand: txt ────────────────────────────────────────────────────────────
def cmd_txt(args):
    out_dir = BASE_DIR / "outputs" / args.o
    tmp_dir = out_dir / "tmp"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    name = args.o
    generate_mesh_from_points(
        args.epi, args.vd, args.ve,
        args.fibbase, args.numfib,
        name, tmp_dir=str(tmp_dir),
    )
    generate_fiber2D(name, args.numfib, tmp_dir=str(tmp_dir))
    run_hexa(name, out_dir, tmp_dir, args.dx, args.dy, args.dz)

    print(f"\nDone.  Results in: outputs/{args.o}/")


# ── CLI ────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    prog="generate_alg.py",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
sub = parser.add_subparsers(dest="cmd", required=True)

# -- mat subcommand --
p_mat = sub.add_parser("mat", help="Generate from Segment .mat file")
p_mat.add_argument("-m",  required=True,  help="Path to .mat file")
p_mat.add_argument("-o",  required=True,  help="Run name (output folder)")
p_mat.add_argument("--slice", type=int, default=None,
                   help="Single 0-indexed slice (omit = all valid slices)")
add_scar_arg(p_mat)
add_resolution_args(p_mat)
p_mat.set_defaults(func=cmd_mat)

# -- txt subcommand --
p_txt = sub.add_parser("txt", help="Generate from plain-text contour files")
p_txt.add_argument("-epi",     required=True)
p_txt.add_argument("-vd",      required=True)
p_txt.add_argument("-ve",      required=True)
p_txt.add_argument("-o",       required=True, help="Run name (output folder)")
p_txt.add_argument("--numfib", type=int, default=0)
p_txt.add_argument("--fibbase",type=str, default="f")
add_scar_arg(p_txt)
add_resolution_args(p_txt)
p_txt.set_defaults(func=cmd_txt)

if __name__ == "__main__":
    args = parser.parse_args()
    args.func(args)
