import argparse
import h5py
import os
import numpy as np
import SimpleITK as sitk
from scipy import io
from skimage.measure import block_reduce

parser = argparse.ArgumentParser(description='Convert .mat into .raw')
parser.add_argument('filename', type=str,
                    help='Path to mat file with volm array')
parser.add_argument('--factor', type=int,
                    help='Downscaling factor', nargs='?', const=1, default=1)

args = parser.parse_args()
filepath, filename = os.path.split(args.filename)
filename = os.path.join(filepath, os.path.splitext(filename)[0])
factor = args.factor

try:
    # if saved wih -v7.3 flag
    f = h5py.File(f'{filename}.mat', 'r')
except OSError:
    # else
    f = io.loadmat(f'{filename}.mat')
    sci = True
except:
    print("Something went wrong")
    raise Exception('Not able to read the file.')

data = f.get('volm')
arr = np.array(data)
if sci:
    arr = np.transpose(arr, axes=[2, 1, 0])
print(f'Original: {arr.shape}')

if (factor > 1):
    arr = block_reduce(arr, (factor, factor, factor),
                       func=np.mean).astype(np.uint16)
    print(f'Downsampled: {arr.shape}')
    filename += f'_downsampled_{factor}'

img = sitk.GetImageFromArray(arr)
writer = sitk.ImageFileWriter()
writer.SetFileName(f'{filename}.mhd')
writer.Execute(img)
