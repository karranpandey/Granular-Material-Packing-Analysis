#load msc etc make sure template same as cluster saddles
#

import pyms3d
import itertools
import numpy as np
from collections import defaultdict
from python_algorithms.basic.union_find import UF
import math
import pickle
import argparse
import SimpleITK as sitk


def check_bd_cell_id(cell_id,mfold_1,mfold_2):
    allowed_thresh=0 #max(0,min_sad_val-0.01) # TODO: Should be 0 ideally to make the most sense
    nb_list=[]
    inc_arr=[1,x_max,y_max*x_max]
    if(cell_id not in mfold_1):
        return False
    if(get_cube_max_val(dual_pts[cell_id])<allowed_thresh):
        return False
    for inc in inc_arr:
        if(cell_id+inc in mfold_2):
            if(get_cube_max_val(dual_pts[cell_id+inc])>allowed_thresh):
                return True
        if(cell_id-inc in mfold_2):
            if(get_cube_max_val(dual_pts[cell_id-inc])>allowed_thresh):
                return True
    return False

def get_nbs_cell_id(cell_id):
    nbd_list=[]
    nbd_arr=list(itertools.product([1,-1,0],repeat=3))
    nbd_arr.remove((0,0,0))
    for arr in nbd_arr:
        inc=(arr[0]+x_max*arr[1]+y_max*x_max*arr[2])
        if((cell_id+inc)>0 and (cell_id+inc)<x_max*y_max*z_max):
            nbd_list.append(cell_id+inc)
    return nbd_list



def get_saddle_cofac_cell_id(s):
    cofac_arr=[]
    int_dual_pt=[x/2 for x in msc.cp_cellid(s)]
    cell_id=int_dual_pt[0]+(x_max)*int_dual_pt[1]+(y_max)*(x_max)*int_dual_pt[2]
    dual_pt=[x/2.0 for x in msc.cp_cellid(s)]
    for i in range(3):
        if(dual_pt[i]%1==0):
            if(i==0):
                cofac_arr.append(cell_id+1)
                cofac_arr.append(cell_id-1)
            if(i==1):
                cofac_arr.append(cell_id+x_max)
                cofac_arr.append(cell_id-x_max)
            if(i==2):
                cofac_arr.append(cell_id+y_max*x_max)
                cofac_arr.append(cell_id-y_max*x_max)
    return cofac_arr


#check for redundant ones

def bfs_intersection_ids(m1,m2,sad):
    visited=[]
    uf_struct=UF(0)
    init_list=[]
    init_list.extend(get_saddle_cofac_cell_id(sad))
    queue=init_list
    bd_int=[]
    explored=set()
    mfold_1=msc.des_geom(m1)
    mfold_2=msc.des_geom(m2)
    while queue:
        node=queue.pop(0)
        if node not in explored:
            explored.add(node)
            neighbours=get_nbs_cell_id(node)
            for neighbour in neighbours:
                if(check_bd_cell_id(neighbour,mfold_1,mfold_2)):
                    queue.append(neighbour)
    return explored


def get_cube_max_val(dual_pt): #actually min val
	inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
	max_val=-1
	#flag=0
	for inc in inc_arr:
		curr_pt=np.add(dual_pt,inc)			#flag=1
		curr_val=msc.vert_func(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
		if(curr_val>max_val):
			max_val=curr_val
	#print(flag)
	return(max_val)



parser=argparse.ArgumentParser()
parser.add_argument('data_file',type=str,help='data file name')
parser.add_argument('dim',nargs=3,type=int,help='dimensions')
parser.add_argument('output_path',type=str,help='output file path')

args=parser.parse_args()
msc_file_name=args.data_file
output_path_name=args.output_path

msc_file_name=args.data_file
dim = list(args.dim)
output_path_name=args.output_path

msc=pyms3d.mscomplex()
msc.load(msc_file_name)
x_max=dim[0]-1
y_max=dim[1]-1
z_max=dim[2]-1
msc.collect_geom(dim=3,dir=0)
cps_2sad=msc.cps(2)
max_conn_dict=defaultdict(list)
dual_pts=msc.dual_points()

surv_sads=np.loadtxt(output_path_name+'surviving_sads')

sad_arr=[]
i=0.0
for sad in surv_sads:
    print((i/len(surv_sads))*100)
    i+=1
    sad=int(sad)
    flag=0
    max_val=0
    sad_dist_val=msc.cp_func(int(sad))
    curr_max_list=[]
    for m,k in msc.asc(int(sad)):
        m=int(m)
        curr_max_list.append(m)
    intersection_cells=bfs_intersection_ids(curr_max_list[0],curr_max_list[1],sad)
    for cell_id in intersection_cells:
        curr_pt=dual_pts[cell_id]
        curr_list=[]
        curr_list.extend(curr_pt)
        curr_list.append(sad)
        sad_arr.append(curr_list)

np.savetxt(output_path_name+'intersection_surfaces_2',sad_arr)
#print(max_removal_dict)


#new_dict_file=open(output_path_name+"merged_max_data_dict.pkl","wb")
#pickle.dump(max_data_dict,new_dict_file)
#np.savetxt(output_path_name+'merged_segmentation',merged_seg)
