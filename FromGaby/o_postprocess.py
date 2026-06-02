import os
import sys
import glob
import datetime
import subprocess
import shutil
import traceback
from pathlib import Path

def fail(message, exc=None, code=1):
    """Print user-friendly error details and exit with a non-zero code."""
    print(f"ERROR: {message}", file=sys.stderr)
    if exc is not None:
        print(f"ERROR DETAIL: {type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc()
    sys.exit(code)


if __name__ == '__main__':
    # Import croco_sigma_to_z from the same directory as this script
    sys.path.insert(0, str(Path(__file__).parent))
    try:
        from croco_sigma_to_z import croco_sigma_to_z
    except Exception as exc:
        fail(
            "Failed to import croco_sigma_to_z. "
            "Check Python dependencies in the croco_forecast environment.",
            exc,
            code=2,
        )

    script_dir = Path(__file__).parent.resolve()
    if len(sys.argv) > 1:
        scratch_dir = Path(sys.argv[1]).resolve()
    else:
        today = datetime.date.today().strftime('%d-%m-%Y')
        scratch_dir = script_dir.parent / today / 'SCRATCH'

    # Original hardcoded paths
    # gpath = '/media/claudiaa/Elements/eDRIVE/02-02-2026/SCRATCH/'
    # path = '/DATA/CROCO/FromDMZ/27-04-2026/SCRATCH'
    # grid_file = os.path.join(gpath, 'croco_grd.nc')

    path = str(scratch_dir)
    grid_file = str(scratch_dir / 'croco_grd.nc')

    print(f"SCRATCH dir: {path}")
    print(f"Grid file:   {grid_file}")

    if not scratch_dir.exists():
        fail(f"SCRATCH directory does not exist: {path}", code=3)
    if not os.path.exists(grid_file):
        fail(f"Grid file not found: {grid_file}", code=3)

    input_files = sorted(glob.glob(os.path.join(path, 'croco_his.*.nc')))
    if not input_files:
        input_files = sorted(glob.glob(os.path.join(path, 'croco_his.nc')))

    if not input_files:
        fail(f"No croco_his*.nc files found in: {path}", code=4)

    for input_file in input_files:
        suffix = os.path.basename(input_file).replace('croco_his', '').replace('.nc', '').strip('.')
        if suffix:
            output_file = os.path.join(path, f'd1_temp_salt_uv_z.{suffix}.nc')
        else:
            output_file = os.path.join(path, 'd1_temp_salt_uv_z.nc')

        print(f"Processing {os.path.basename(input_file)} -> {os.path.basename(output_file)}")
        try:
            croco_sigma_to_z(input_file, output_file, grid_file)
        except Exception as exc:
            fail(
                f"sigma-to-z conversion failed for input file: {input_file}",
                exc,
                code=5,
            )

    output_files = sorted(glob.glob(os.path.join(path, 'd1_temp_salt_uv_z.*.nc')))
    if not output_files:
        # single file with no suffix
        output_files = sorted(glob.glob(os.path.join(path, 'd1_temp_salt_uv_z.nc')))
    concat_output = os.path.join(path, 'd1_temp_salt_uv_z_all.nc')

    if not output_files:
        fail(
            f"No sigma-to-z output files were created in: {path}",
            code=6,
        )

    if len(output_files) <= 1:
        print("Only one output file found; skipping concatenation.")
        try:
            shutil.copy2(output_files[0], concat_output)
            print(f"Created {os.path.basename(concat_output)} from single output file.")
        except Exception as exc:
            fail(
                f"Failed to create {concat_output} from single output file.",
                exc,
                code=7,
            )
    else:
        print(f"Concatenating {len(output_files)} files -> {os.path.basename(concat_output)}")
        # Make t a record dimension in each file so ncrcat can concatenate
        rec_files = []
        try:
            for f in output_files:
                rec = f.replace('.nc', '_rec.nc')
                subprocess.run(['ncks', '-O', '--mk_rec_dmn', 't', f, rec], check=True)
                rec_files.append(rec)
            cmd = ['ncrcat'] + rec_files + [concat_output]
            subprocess.run(cmd, check=True)
        except FileNotFoundError as exc:
            fail(
                "Required NCO command not found (ncks or ncrcat). "
                "Install NCO in the container image.",
                exc,
                code=8,
            )
        except subprocess.CalledProcessError as exc:
            fail(
                "NCO command failed during concatenation.",
                exc,
                code=8,
            )
        except Exception as exc:
            fail(
                "Unexpected error during concatenation.",
                exc,
                code=8,
            )
        finally:
            for f in rec_files:
                if os.path.exists(f):
                    os.remove(f)
