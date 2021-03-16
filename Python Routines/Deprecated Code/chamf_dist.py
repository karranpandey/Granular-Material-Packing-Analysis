import itk
import SimpleITK as sitk
import numpy as np
import argparse

parser=argparse.ArgumentParser()

parser.add_argument('data_file',type=str,help='data file name')
args=parser.parse_args()

main_file_name = args.data_file

#main_file_name="chan_vese_standard_Steel_Deposition_362_352_502.mhd"

raw_writer = sitk.ImageFileWriter()

reader=sitk.ImageFileReader()
reader.SetImageIO("MetaImageIO")
reader.SetFileName(main_file_name)

image=reader.Execute()

arr_image=sitk.GetArrayFromImage(image)

itk_image=itk.GetImageFromArray(np.array(arr_image).astype(np.float32))

antialiasfilter = itk.AntiAliasBinaryImageFilter.New(itk_image)
antialiasfilter.SetInput(itk_image)
antialiasfilter.Update()
antialias_image=antialiasfilter.GetOutput()

isoContourFilter=itk.IsoContourDistanceImageFilter.New(antialias_image)

isoContourFilter.SetLevelSetValue(0.1)

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

raw_writer = sitk.ImageFileWriter()
raw_writer.SetFileName('chamf_distance_'+main_file_name)
raw_writer.Execute(dist_image)
