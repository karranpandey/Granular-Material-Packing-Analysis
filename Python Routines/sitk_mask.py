import SimpleITK as sitk
import random
import numpy as np

reader = sitk.ImageFileReader()
reader.SetImageIO("MetaImageIO")
reader.SetFileName("../Outputs/angular_downsampled.mhd")
img = reader.Execute()
img = sitk.GetArrayFromImage(img)
dims = img.shape
max_val = np.max(img)
for x in range(dims[0]):
    for y in range(dims[1]):
            for z in range(dims[2]):
                    if(((y-79)**2 + (z-79)**2)**0.5>76):
                            img[x,y,z]=img[x,y,z]/max_val*30000
img = sitk.GetImageFromArray(img)
writer = sitk.ImageFileWriter()
writer.SetFileName('../Outputs/masked_angular_downsampled.mhd')
writer.Execute(img)

