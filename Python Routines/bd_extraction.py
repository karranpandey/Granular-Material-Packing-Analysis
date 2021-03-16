import skimage.io as io
from skimage.segmentation import morphological_chan_vese
import SimpleITK as sitk
import itk
import numpy as np
import argparse
import os

def bd_extraction(image):

    print('Commencing Boundary Extraction')

    otsuFilter=sitk.OtsuThresholdImageFilter()
    otsuFilter.Execute(image)
    otsu_thresh=otsuFilter.GetThreshold()

    print('Otsu Threshold Computed')
    print(otsu_thresh)

    dist_image=image - otsu_thresh
    dist_image=sitk.GetArrayFromImage(dist_image)

    image=sitk.GetArrayFromImage(image)

    total_num_points=(image>0).sum()
    inside_num_points=(image>otsu_thresh).sum()
    max_dens_val=np.max(image)

    lambda_1=total_num_points*max_dens_val/(inside_num_points*(max_dens_val-otsu_thresh))
    lambda_2=total_num_points*max_dens_val/((total_num_points-inside_num_points)*(otsu_thresh))

    print('Commencing Chan and  Vese with Lambda 1 and Lambda 2 as follows:')
    print(lambda_1)
    print(lambda_2)
    ls = morphological_chan_vese(image, 20, init_level_set=dist_image, lambda1=abs(lambda_1),lambda2=abs(lambda_2),smoothing=0)
    print('Boundary Extracted.')
    return ls

def dist_field_comp(ls):

    print('Commencing Distance Field Computation')

    itk_image=itk.GetImageFromArray(np.array(ls).astype(np.float64))

    antialiasfilter = itk.AntiAliasBinaryImageFilter.New(itk_image)
    antialiasfilter.SetInput(itk_image)
    antialiasfilter.Update()
    antialias_image=antialiasfilter.GetOutput()

    isoContourFilter=itk.IsoContourDistanceImageFilter.New(antialias_image)
    isoContourFilter.SetLevelSetValue(0.5)
    isoContourFilter.SetFarValue(100) #make sure this is greater than your max distance value for the chamfer filter
    isoContourFilter.SetInput(antialias_image)
    isoContourFilter.Update()
    isoContour_image=isoContourFilter.GetOutput()

    chamferFilter=itk.FastChamferDistanceImageFilter.New(isoContour_image)
    chamferFilter.SetMaximumDistance(50.0)
    chamferFilter.Update()
    chamferFilter.SetInput(isoContour_image)
    chamferFilter.Update()
    chamf_image=chamferFilter.GetOutput()
    chamf_arr=itk.GetArrayFromImage(chamf_image)
    dist_image=sitk.GetImageFromArray(chamf_arr)

    print('Distance Field Computed')

    return dist_image


parser=argparse.ArgumentParser()
parser.add_argument('data_file',type=str,help='data file name')
#parser.add_argument('output_path',type=str,help='output file path')
args=parser.parse_args()

main_file_name = args.data_file

base_name=os.path.basename(main_file_name)
base_name=os.path.splitext(base_name)[0]

reader=sitk.ImageFileReader()
reader.SetImageIO("MetaImageIO")
reader.SetFileName(main_file_name)
image=reader.Execute()

ls=bd_extraction(image)
dist_image=dist_field_comp(ls)
#dist_image = sitk.GetImageFromArray(ls)

print('Writing File')

raw_writer = sitk.ImageFileWriter()
output_path_name='../Outputs/'
raw_writer.SetFileName(output_path_name+'chamf_distance_'+base_name+'.mhd')
raw_writer.Execute(dist_image)

print('File Written')
