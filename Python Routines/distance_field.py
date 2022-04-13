import argparse
import numpy as np
import os
import SimpleITK as sitk
from utilities import read_mat_file
from utilities import bd_extraction
from utilities import dist_field_comp

# positional arguments for the command line
parser = argparse.ArgumentParser()
parser.add_argument('data_file', type=str, help='data file name')
parser.add_argument('factor', type=int, help='downscaling factor')
filtername = 'Adaptive'
args = parser.parse_args()

# get the inputs
main_file_name, factor = args.data_file, args.factor

# name of the file without pathname and extension
base_name = os.path.basename(main_file_name)
base_name = os.path.splitext(base_name)[0]

# full file name without extension
matfilename = os.path.splitext(main_file_name)[0]
# read the mat file and get the data -- optionally downsample
arr = read_mat_file(matfilename, factor)

# adjust data range (contrast adjustment - imadjust)
low, upp, typ = np.quantile(arr, 0.01), np.quantile(arr, 0.99), arr.dtype
arr = (2**16 - 1) * (arr - low) / (upp - low)
arr = (np.clip(arr, 0, 2**16-1)).astype(typ)

# write the downsampled file for visualization
sitk.WriteImage(sitk.GetImageFromArray(arr), matfilename + '.mhd')

# get the binary volume
ls = bd_extraction(arr, slicewise=True, filterName=filtername)
# get the distance field
dist_field = dist_field_comp(ls)

print('Writing File')
output_path_name = '../ChamferDistance/'
if not os.path.exists(output_path_name):
    os.makedirs(output_path_name)
sitk.WriteImage(sitk.GetImageFromArray(dist_field),
                output_path_name+'chamf_distance_'+base_name+'.mhd')

print('File Written')
