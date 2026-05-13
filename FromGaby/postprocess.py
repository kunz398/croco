import os
import sys
import glob
import datetime
import subprocess
from pathlib import Path

# Import croco_sigma_to_z from the same directory as this script
sys.path.insert(0, str(Path(__file__).parent))
from croco_sigma_to_z import croco_sigma_to_z


if __name__ == '__main__':
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

    input_files = sorted(glob.glob(os.path.join(path, 'croco_his.*.nc')))
    if not input_files:
        input_files = sorted(glob.glob(os.path.join(path, 'croco_his.nc')))

    if not input_files:
        print("No croco_his.*.nc files found. Exiting.")
        sys.exit(1)

    for input_file in input_files:
        suffix = os.path.basename(input_file).replace('croco_his', '').replace('.nc', '').strip('.')
        if suffix:
            output_file = os.path.join(path, f'd1_temp_salt_uv_z.{suffix}.nc')
        else:
            output_file = os.path.join(path, 'd1_temp_salt_uv_z.nc')

        print(f"Processing {os.path.basename(input_file)} -> {os.path.basename(output_file)}")
        ds_out = croco_sigma_to_z(input_file, output_file, grid_file)

    output_files = sorted(glob.glob(os.path.join(path, 'd1_temp_salt_uv_z.*.nc')))
    if not output_files:
        # single file with no suffix
        output_files = sorted(glob.glob(os.path.join(path, 'd1_temp_salt_uv_z.nc')))
    concat_output = os.path.join(path, 'd1_temp_salt_uv_z_all.nc')

    print(f"Concatenating {len(output_files)} files -> {os.path.basename(concat_output)}")
    if len(output_files) == 1:
        import shutil
        shutil.copy2(output_files[0], concat_output)
        print(f"Single file: copied to {os.path.basename(concat_output)}")
    else:
        # Make t a record dimension in each file so ncrcat can concatenate
        rec_files = []
        for f in output_files:
            rec = f.replace('.nc', '_rec.nc')
            subprocess.run(['ncks', '-O', '--mk_rec_dmn', 't', f, rec], check=True)
            rec_files.append(rec)
        cmd = ['ncrcat'] + rec_files + [concat_output]
        subprocess.run(cmd, check=True)
        for f in rec_files:
            os.remove(f)
