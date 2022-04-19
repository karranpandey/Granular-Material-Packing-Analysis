import SimpleITK as sitk
from scipy.io import savemat
import argparse
import os


def from_raw_2_np(filename):
    '''
    The filename should be .mhd associated with the .raw file
    you want to convert into the numpy array.
    '''
    reader = sitk.ImageFileReader()
    reader.SetImageIO("MetaImageIO")
    reader.SetFileName(filename)

    image = reader.Execute()
    return sitk.GetArrayFromImage(image)


parser = argparse.ArgumentParser(description='Convert .raw into .mat')
parser.add_argument('file_path', metavar='filename', type=str, nargs=1,
                    help='Path to the .mhd file')

args = parser.parse_args()
filepath = args.file_path[0]
filename = os.path.basename(filepath)
# filename = os.path.splitext(filename)[0]


arr = from_raw_2_np(f'{filename}')
my_var_dict = dict(volm=arr)
savemat(f'{filename}.mat', my_var_dict)
