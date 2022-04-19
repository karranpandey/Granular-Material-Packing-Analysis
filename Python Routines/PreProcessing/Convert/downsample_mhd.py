from skimage.measure import block_reduce
import SimpleITK as sitk
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('data_file', type=str, help='data file name')
parser.add_argument('factor', type=int, help='downsampling factor')
args = parser.parse_args()
data_file_name = args.data_file
factor = args.factor

fac = (factor, factor, factor)

reader = sitk.ImageFileReader()
reader.SetImageIO("MetaImageIO")
reader.SetFileName(data_file_name)
img_arr = reader.Execute()

img_arr = sitk.GetArrayFromImage(img_arr)
print(np.shape(img_arr))
downsampled_img = block_reduce(img_arr, fac, func=np.mean).astype(np.uint16)
print(np.shape(downsampled_img))

downsampled_img = sitk.GetImageFromArray(downsampled_img)

output_path_name = './'
if not os.path.exists(output_path_name):
    os.makedirs(output_path_name)

base_name = os.path.basename(data_file_name)
base_name = os.path.splitext(base_name)[0]

raw_writer = sitk.ImageFileWriter()
raw_writer.SetFileName(output_path_name+'downsampled_'+base_name+'.mhd')
raw_writer.Execute(downsampled_img)
