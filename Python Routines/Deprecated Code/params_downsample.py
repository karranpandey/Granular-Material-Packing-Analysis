import dxchange
import numpy as np
from skimage.measure import block_reduce
import SimpleITK as sitk
import os

#slice by image height

file_path='Depositional_Steel_Air.txm'

file_name=os.path.splitext(file_path)[0]

txm_img_params=dxchange.read_txrm(file_path,[1,1,1])

max_dim_x=txm_img_params[1]['number_of_images'] #1004 (X)

max_dim_y=txm_img_params[1]['image_height'] #1024 (Y)

max_dim_z=txm_img_params[1]['image_width'] #1004 (Z)

#max_dim_x=900

factor=4

step=max_dim_y/factor
curr_step=0
slice_arr=[]

for i in range(factor):
	next_step=curr_step+step
	if(next_step>max_dim_y):
		next_step=max_dim_y
	slice_arr.append([(0,max_dim_x,1),(curr_step,next_step,1),(0,max_dim_z,1)])
	curr_step=next_step

print(slice_arr)

full_ds_img=[]

down_fac=(2,2,2)

for i in range(len(slice_arr)):
	txm_img=dxchange.read_txrm(file_path,slice_arr[i])
	txm_img=txm_img[0]
	#txm_img=txm_img[50:max_dim_x]
	print('Slice '+str(i)+' read')
	print(np.shape(txm_img))
	print(np.max(txm_img))
	print(np.min(txm_img))
	downsampled_img=block_reduce(txm_img,down_fac,func=np.mean)
	print('Slice '+str(i)+' Downsampled')
	print(np.shape(downsampled_img))
	print(np.max(downsampled_img))
	print(np.min(downsampled_img))
	if(i==0):
		full_ds_img=downsampled_img
	else:
		full_ds_img=np.concatenate((full_ds_img,downsampled_img),axis=1)
	print(np.shape(full_ds_img))

	
output_file_name=file_name+'_downsampled_mean_'+str(np.shape(full_ds_img))+'.mhd'

output_file_name=output_file_name.replace('(','')
output_file_name=output_file_name.replace(')','')
output_file_name=output_file_name.replace(',','_')
output_file_name=output_file_name.replace(' ','')
rawWriter=sitk.ImageFileWriter()
rawWriter.SetFileName(output_file_name)
rawWriter.Execute(sitk.GetImageFromArray(full_ds_img))




