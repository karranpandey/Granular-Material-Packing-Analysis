import pyms3d
import itertools
import numpy as np
from collections import defaultdict
from python_algorithms.basic.union_find import UF
import math
import pickle
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('data_file',type=str,help='data file name')
parser.add_argument('dim',nargs=3,type=int,help='dimensions')
parser.add_argument('output_path',type=str,help='output file path')

args=parser.parse_args()
msc_file_name=args.data_file
dim = list(args.dim)
output_path_name=args.output_path

msc=pyms3d.mscomplex()
msc.load(msc_file_name)
msc.collect_geom(dim=3,dir=0)
cps_2sad=msc.cps(2)
x_max=dim[0]-1
y_max=dim[1]-1
z_max=dim[2]-1
dual_pts=msc.dual_points()

particle_segmentation=[]

surv_sads=np.loadtxt(output_path_name+'surviving_sads')

def get_cube_max_val(dual_pt): #actually min val
	inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
	max_val=100
	#flag=0
	for inc in inc_arr:
		curr_pt=np.add(dual_pt,inc)
		if(curr_pt[0]<x_max and curr_pt[1]<y_max and curr_pt[2]<z_max):
			#flag=1
			curr_val=msc.vert_func(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
			if(curr_val<max_val):
				max_val=curr_val
	#print(flag)
	return(max_val)

def get_des_vol(max_id):
	des_geom=msc.des_geom(max_id)
	count=0
	for cube_id in des_geom:
		if(get_cube_max_val(dual_pts[cube_id])>0):
			count+=1
			particle_segmentation.append(np.append(dual_pts[cube_id],(max_id)))
	return count

def dia_from_vol(vol):
	dia=(vol*3.0)/(4.0*math.pi)
	dia=dia**(1/3.0)
	return 2.0*dia

cps_max=msc.cps(3)

#point ids
#C. Number
#Vol
#Assoc Radius

def get_cnum_max(max_id):
	count=0
	for sad_id,k in msc.des(max_id):
		if(sad_id in surv_sads):
			count+=1
	return count

max_data_array=[]
max_data_dict=defaultdict(list)

for i in range(len(cps_max)):
	max_id=cps_max[i]
	curr_max_data=[]
	co_ords=[x/2.0 for x in msc.cp_cellid(max_id)]
	vol=get_des_vol(max_id)
	dia=dia_from_vol(vol)
	c_num=get_cnum_max(max_id)
	curr_max_data.extend(co_ords)
	curr_max_data.append(vol)
	curr_max_data.append(dia)
	curr_max_data.append(c_num)
	max_data_dict[max_id]=curr_max_data
	curr_max_data.append(max_id)
	max_data_array.append(curr_max_data)

dict_file=open(output_path_name+"max_data_dict.pkl","wb")
pickle.dump(max_data_dict,dict_file)
np.savetxt(output_path_name+'max_data_array',max_data_array)
np.savetxt(output_path_name+'segmentation',particle_segmentation)
