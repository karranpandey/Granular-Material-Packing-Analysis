import pyms3d
import itertools
import numpy as np
from collections import defaultdict
from python_algorithms.basic.union_find import UF
import math
import pickle
import argparse
import SimpleITK as sitk

def dia_from_vol(vol):
    dia=(vol*3.0)/(4.0*math.pi)
    dia=dia**(1/3.0)
    return 2.0*dia

def get_cube_max_val(dual_pt): #actually min val
	inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
	max_val=100
	#flag=0
	for inc in inc_arr:
		curr_pt=np.add(dual_pt,inc)			#flag=1
		curr_val=msc.vert_func(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
		if(curr_val<max_val):
			max_val=curr_val
	#print(flag)
	return(max_val)


def get_sad_val(sad_id):
    sad_pt=[x/2.0 for x in msc.cp_cellid(sad_id)]
    fac_arr=[]
    max_val=0
    poss_coords=defaultdict(list)
    for i in range(3):
        if(sad_pt[i]%1!=0):
            poss_coords[i].append(sad_pt[i]+0.5)
            poss_coords[i].append(sad_pt[i]-0.5)
        else:
            poss_coords[i].append(sad_pt[i])
    fac_arr=list(itertools.product(poss_coords[0],poss_coords[1],poss_coords[2]))
    for fac in fac_arr:
        curr_val=image.GetPixel(int(fac[0]),int(fac[1]),int(fac[2]))
        if(curr_val>max_val):
            max_val=curr_val
    return(max_val)

def get_range_vol(max_id):
    inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
    max_val=0
    min_val=10000000
    avg_val=0
    count=0
    des_geom=msc.des_geom(max_id)
    for cube_id in des_geom:
        curr_cube=np.array(dual_pts[cube_id])
        for inc in inc_arr:
            curr_pt=np.add(curr_cube,inc)
            curr_val=image.GetPixel(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
            avg_val+=curr_val
            count+=1
            if(curr_val>max_val):
                max_val=curr_val
            if(curr_val<min_val):
                min_val=curr_val
    print(max_val)
    print(min_val)
    print(avg_val/count)
    return (max_val-min_val)

def get_avg_val_vol(max_id,test_val):
    inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
    avg_val=0
    count=0
    des_geom=msc.des_geom(max_id)
    for cube_id in des_geom:
        curr_cube=np.array(dual_pts[cube_id])
        for inc in inc_arr:
            curr_pt=np.add(curr_cube,inc)
            curr_val=image.GetPixel(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
            if(curr_val>test_val):
                avg_val+=curr_val
                count+=1
    return (avg_val/count)


def get_max_val(max_id):
    max_pt=[x/2.0 for x in msc.cp_cellid(max_id)]
    inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
    max_val=0
    for inc in inc_arr:
        curr_pt=np.add(max_pt,inc)
        curr_val=image.GetPixel(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
        if(curr_val>max_val):
            max_val=curr_val
    return(max_val)

def get_angle(a,b,c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return(np.degrees(angle))

def get_max_dist(max_id,sad_id):
    max_pt=np.array([x/2.0 for x in msc.cp_cellid(max_id)])
    sad_pt=np.array([x/2.0 for x in msc.cp_cellid(sad_id)])
    des_geom=msc.des_geom(max_id)
    max_dist=0
    for cube_id in des_geom:
        curr_pt=np.array(dual_pts[cube_id])
        if(get_cube_max_val(curr_pt)>0):
            if(89<get_angle(curr_pt,max_pt,sad_pt)<91):
                curr_dist=np.linalg.norm(max_pt-dual_pts[cube_id])
                if(curr_dist>max_dist):
                    max_dist=curr_dist
    return max_dist

def get_max_rad(max_id):
    max_pt=np.array([x/2.0 for x in msc.cp_cellid(max_id)])
    des_geom=msc.des_geom(max_id)
    max_dist=0
    for i in range(len(des_geom)):
        cube_id=des_geom[i]
        curr_pt=np.array(dual_pts[cube_id])
        if(get_cube_max_val(curr_pt)<0):
            continue
        #print(i)
        for j in range(i+1,len(des_geom)):
            alt_cube_id=des_geom[j]
            alt_pt=np.array(dual_pts[alt_cube_id])
            if(get_cube_max_val(alt_pt)<0):
                continue
            curr_dist=np.linalg.norm(curr_pt-alt_pt)
            if(curr_dist>max_dist):
                max_dist=curr_dist
    return max_dist/2.0



parser=argparse.ArgumentParser()
parser.add_argument('data_file',type=str,help='data file name')
parser.add_argument('img_file',type=str,help='data file name')
parser.add_argument('output_path',type=str,help='output file path')

args=parser.parse_args()
msc_file_name=args.data_file
img_file_name=args.img_file
output_path_name=args.output_path

reader=sitk.ImageFileReader()
reader.SetImageIO("MetaImageIO")
reader.SetFileName(img_file_name)
image=reader.Execute()

otsuFilter=sitk.OtsuThresholdImageFilter()
otsuFilter.Execute(image)
otsu_thresh=otsuFilter.GetThreshold()

print(otsu_thresh)

msc=pyms3d.mscomplex()
msc.load(msc_file_name)
msc.collect_geom(dim=3,dir=0)
dual_pts=msc.dual_points()

max_data_dict=defaultdict(list)

dict_file=open(output_path_name+"max_data_dict.pkl","rb")
max_data_dict=pickle.load(dict_file)

particle_segmentation=np.loadtxt(output_path_name+'segmentation')

surv_sads=np.loadtxt(output_path_name+'surviving_sads')

remove_max=[]
max_removal_dict={}
dual_max_vis=[]

for sad in surv_sads:
    sad=int(sad)
    flag=0
    max_val=0
    sad_dist_val=msc.cp_func(int(sad))
    for m,k in msc.asc(int(sad)):
        m=int(m)
        if(flag!=0):
            if(msc.cp_func(m)<msc.cp_func(prev_max)):
                max_val=msc.cp_func(m)
                sad_val=msc.cp_func(sad)
		print(sad_val/max_val)
                #print(max_val-sad_val)
                #print(sad_val-otsu_thresh)
                if(sad_val/max_val>0.75):
                    #max_dist=get_max_rad(m)
                    #print(sad_val/max_dist)
                    #if(sad_val/max_dist>0.5):
                    remove_max.append(m)
                    max_removal_dict[m]=prev_max
                        #dual_max_vis.append(prev_max)
            else:
                max_val=msc.cp_func(prev_max)
                sad_val=msc.cp_func(sad)
		print(sad_val/max_val)
                if(sad_val/max_val>0.75):
                    #max_dist=get_max_rad(prev_max)
                    #print(sad_val/max_dist)
                    #if(sad_val/max_dist>0.5):
                    remove_max.append(prev_max)
                    max_removal_dict[prev_max]=m
                    #dual_max_vis.append(m)
        prev_max=m
        flag=1

print(remove_max) 
print(np.size(remove_max))
new_dict_file=open(output_path_name+"remove_max_dict.pkl","wb")
pickle.dump(max_removal_dict,new_dict_file)

