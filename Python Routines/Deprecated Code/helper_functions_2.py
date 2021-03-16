import pyms3d
import itertools
import numpy as np
from collections import defaultdict
from python_algorithms.basic.union_find import UF
import math

msc=pyms3d.mscomplex()
#msc.load('msc_otsu_Depositional_Steel_Air_downsampled_mean_213_256_251_pers_2.0')
msc.load('msc_chamfer_dist_pers_2.0')
msc.collect_geom(dim=3,dir=0)
cps_2sad=msc.cps(2)
max_conn_dict=defaultdict(list)
dual_pts=msc.dual_points()

for s in cps_2sad:
	if(msc.cp_func(s)<0):
		continue
	max_list=[]
    	for m,k in msc.asc(s):
        	max_list.append(m)
    	max_pairs=list(itertools.combinations(max_list,2))
    	for pair in max_pairs:
            	max_conn_dict[tuple(pair)].append(s)


def get_cube_max_val(dual_pt):
	inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
	max_val=0
	for inc in inc_arr:
		curr_pt=np.add(dual_pt,inc)
		curr_val=msc.vert_func(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
		if(curr_val>max_val):
			max_val=curr_val
	return(max_val)
			
def get_saddle_cofacets(s):
	cofac_arr=[]
	dual_pt=[x/2.0 for x in msc.cp_cellid(s)]
	inc_arr=[[0.5,0,0],[0,0.5,0],[0,0,0.5]]
	for i in range(3):	
		if(dual_pt[i]%1==0):
			cofac_arr.append(np.add(dual_pt,inc_arr[i]))
			cofac_arr.append(np.subtract(dual_pt,inc_arr[i]))
	return cofac_arr

def pos_des_mfold(m):
	mfold_pt_list=[]
	for k in msc.des_geom(m):
		curr_pt=dual_pts[k]
		if(get_cube_max_val(curr_pt)>0):
			mfold_pt_list.append(curr_pt)
			mfold_dict[tuple(curr_pt)]=m
	print(np.size(mfold_pt_list))
	return mfold_pt_list

def get_nbs_cell_id(cell_id):
	nbd_list=[]
	nbd_arr=list(itertools.product([1,-1,0],repeat=3))
	nbd_arr.remove((0,0,0))
	for arr in nbd_arr:
		inc=(arr[0]+250*arr[1]+255*250*arr[2])
		if((cell_id+inc)>0 and (cell_id+inc)<255*250*212):
			nbd_list.append(cell_id+inc)
	return nbd_list	

def get_sad_facs(sad_id):
	sad_pt=[x/2.0 for x in msc.cp_cellid(sad_id)]
	fac_arr=[]
	dim=[251,256,213]
	poss_coords=defaultdict(list)
	for i in range(3):
		if(sad_pt[i]%1!=0):
			poss_coords[i].append(sad_pt[i]+0.5)
			poss_coords[i].append(sad_pt[i]-0.5)
		else:
			poss_coords[i].append(sad_pt[i])
	fac_arr=list(itertools.product(poss_coords[0],poss_coords[1],poss_coords[2]))
	return(fac_arr)

def get_sad_nbs(sad_id):
	sad_pt=[x/2.0 for x in msc.cp_cellid(sad_id)]
	fac_arr=[]
	nbd_arr=[]
	dim=[251,256,213]
	poss_coords=defaultdict(list)
	nbd_coords=defaultdict(list)
	for i in range(3):
		if(sad_pt[i]%1!=0):
			poss_coords[i].append(sad_pt[i]+0.5)
			poss_coords[i].append(sad_pt[i]-0.5)
		else:
			poss_coords[i].append(sad_pt[i])
	#print(poss_coords)
	for i in range(3):
		if(len(poss_coords[i])>1):
			if(min(poss_coords[i])-1>0):
				nbd_coords[i].append(min(poss_coords[i])-1)
			if(max(poss_coords[i])+1<dim[i]):
				nbd_coords[i].append(max(poss_coords[i])+1)
			if(i==0):
				nbd_arr.extend(list(itertools.product(nbd_coords[0],poss_coords[1],poss_coords[2])))
			if(i==1):
				nbd_arr.extend(list(itertools.product(poss_coords[0],nbd_coords[1],poss_coords[2])))
			if(i==2):
				nbd_arr.extend(list(itertools.product(poss_coords[0],poss_coords[1],nbd_coords[2])))		
	#fac_arr=list(itertools.product(poss_coords[0],poss_coords[1],poss_coords[2]))
	return nbd_arr

def check_act_sad(sad_id):
	sad_pt=[x/2.0 for x in msc.cp_cellid(sad_id)]
	sad_fn_val=msc.cp_func(sad_id)
	nbd_arr=get_sad_nbs(sad_id)
	for nbd_pt in nbd_arr:
		if(msc.vert_func(int(nbd_pt[0]),int(nbd_pt[1]),int(nbd_pt[2]))>sad_fn_val):
			return False
	return True

def check_sad_inside(sad_id):
	for fac in get_sad_facs(sad_id):
		if(msc.vert_func(int(fac[0]),int(fac[1]),int(fac[2]))<0):
			return False
	return True

def show_sad_nbs(sad_id):
	print(msc.cp_func(sad_id))
	sad_pt=[x/2.0 for x in msc.cp_cellid(sad_id)]
	print(sad_pt)	
	for nb in get_sad_nbs(sad_id):
		print(nb)
		print(msc.vert_func(int(nb[0]),int(nb[1]),int(nb[2])))

def show_sad_facs(sad_id):
	print(msc.cp_func(sad_id))
	sad_pt=[x/2.0 for x in msc.cp_cellid(sad_id)]
	print(sad_pt)	
	for nb in get_sad_facs(sad_id):
		print(nb)
		print(msc.vert_func(int(nb[0]),int(nb[1]),int(nb[2])))

total_count=0
for key in max_conn_dict.keys():
	sad_count=len(max_conn_dict[key])
	if(sad_count>1):
		total_count+=1
			
act_count=0		
for key in max_conn_dict.keys():
	sad_count=len(max_conn_dict[key])
	if(sad_count>1):
		for sad in max_conn_dict[key]:
			if(check_act_sad(sad)==False):
				sad_count-=1
		#if(sad_count==0):
			#print(key)
			#print(sad_count)
		if(sad_count>1):
			act_count+=1
			#print(key)
			#print(sad_count)

inside_count=0			
for key in max_conn_dict.keys():
	sad_count=len(max_conn_dict[key])
	if(sad_count>1):
		for sad in max_conn_dict[key]:
			if(check_sad_inside(sad)==False):
				sad_count-=1
		#if(sad_count==0):
		#	print(key)
		#	print(sad_count)
		if(sad_count>1):
			inside_count+=1
			print(key)
			print(sad_count)
			print(max_conn_dict[key])
			

inside_act_count=0		
for key in max_conn_dict.keys():
	sad_count=len(max_conn_dict[key])
	if(sad_count>1):
		for sad in max_conn_dict[key]:
			act_flag=0
			if(check_sad_inside(sad)==False):
				sad_count-=1
			elif(check_act_sad(sad)==False):
				sad_count-=1
				act_flag=1
		if(sad_count==0 and act_flag==1):
			print(key)
			print(sad_count)
			print(max_conn_dict[key])
		if(sad_count>1):
			inside_act_count+=1
			#print(key)
			#print(max_conn_dict[key])
			#print(sad_count)


def check_bd_cell_id(cell_id,mfold_1,mfold_2):
	nb_list=[]
	inc_arr=[1,250,255*250]
	if(cell_id not in mfold_1):
		return False
	if(get_cube_max_val(dual_pts[cell_id])<0):
		return False
	for inc in inc_arr:
		if(cell_id+inc in mfold_2):
			if(get_cube_max_val(dual_pts[cell_id+inc])>0):
				return True
		if(cell_id-inc in mfold_2):
			if(get_cube_max_val(dual_pts[cell_id-inc])>0):
				return True
	return False

def get_saddle_cofac_cell_id(s):
	cofac_arr=[]
	int_dual_pt=[x/2 for x in msc.cp_cellid(s)]
	cell_id=int_dual_pt[0]+250*int_dual_pt[1]+255*250*int_dual_pt[2]
	dual_pt=[x/2.0 for x in msc.cp_cellid(s)]	
	for i in range(3):	
		if(dual_pt[i]%1==0):
			if(i==0):
				cofac_arr.append(cell_id+1)
				cofac_arr.append(cell_id-1)
			if(i==1):
				cofac_arr.append(cell_id+250)
				cofac_arr.append(cell_id-250)
			if(i==2):
				cofac_arr.append(cell_id+255*250)
				cofac_arr.append(cell_id-255*250)
	return cofac_arr

def uf_add(uf_struct):
	uf_struct._id.append(len(uf_struct._id))
	uf_struct._count+=1
	uf_struct._rank.append(0)
	return uf_struct


#check for redundant ones
def init_uf_sads(uf_struct,init_list,uf_id_dict):
	for elem in init_list:
		if(elem in uf_id_dict.keys()):
			continue
		count=len(uf_struct._id)
		uf_id_dict[elem]=count
		uf_add(uf_struct)
		if(count%2==1):
			uf_struct.union(count,count-1)
	return uf_struct

def uf_add_node(node,elem,uf_struct,uf_id_dict):
	if elem in uf_id_dict.keys():
		uf_struct.union(uf_id_dict[elem],uf_id_dict[node])
 	else:
		uf_id_dict[elem]=len(uf_struct._id)
		uf_add(uf_struct)
		uf_struct.union(uf_id_dict[elem],uf_id_dict[node])
	return uf_struct
	

def bfs_intersection_ids(m1,m2):
	visited=[]
	uf_struct=UF(0)
	uf_id_dict={}
	init_list=[]
	for sad in max_conn_dict[(m1,m2)]:
		init_list.extend(get_saddle_cofac_cell_id(sad))
	init_uf_sads(uf_struct,init_list,uf_id_dict)
	queue=init_list
	bd_int=[]
	explored=set()
	mfold_1=msc.des_geom(m1)
	mfold_2=msc.des_geom(m2)
	while queue:
		if(uf_struct.count()==1):
			return 1
		node=queue.pop(0)
		if node not in explored:
			explored.add(node)
			neighbours=get_nbs_cell_id(node)			
			for neighbour in neighbours:
				if(check_bd_cell_id(neighbour,mfold_1,mfold_2)):
					queue.append(neighbour)
					uf_add_node(node,neighbour,uf_struct,uf_id_dict)
	#print(uf_struct.count())
	return uf_struct.count()



def get_neighbours(dual_pt):
	nbd_list=[]
	nbd_arr=list(itertools.product([1,-1,0],repeat=3))
	nbd_arr.remove((0,0,0))
	for arr in nbd_arr:
		nbd_list.append(np.add(dual_pt,list(arr)))
	return nbd_list

def check_bd(dual_pt,mfold_1,mfold_2):
	nb_list=[]
	nbd_arr=[[1,0,0],[0,1,0],[0,0,1]]
	if(tuple(dual_pt) not in mfold_dict.keys()):
		return False
	if(mfold_dict[tuple(dual_pt)]!=mfold_1):
		return False
	for arr in nbd_arr:
		if(tuple(np.add(dual_pt,arr)) not in mfold_dict.keys()):
			continue
		if(tuple(np.subtract(dual_pt,arr)) not in mfold_dict.keys()):
			continue
		if(mfold_dict[tuple(np.add(dual_pt,arr))]==mfold_2):
			return True
		if(mfold_dict[tuple(np.subtract(dual_pt,arr))]==mfold_2):
			return True
	return False

	
#have to either fix m1,m2 here or call accordingly.

def bfs_intersection(init_list,m1,m2):
	visited=[]
	queue=init_list
	bd_int=[]
	explored=set()
	while queue:
		node=queue.pop(0)
		if tuple(node) not in explored:
			explored.add(tuple(node))
			neighbours=get_neighbours(node)			
			for neighbour in neighbours:
				if(check_bd(neighbour,m1,m2)):
					queue.append(tuple(neighbour))
	return explored
				

def clean_saddles(max_pair,sad_list):
	populate_max_mfold_flags(max_pair[0])
	populate_max_mfold_flags(max_pair[1])
	for sad in sad_list:
		populate_sad_cofac_flags(sad)
	bd_comps=[]
	for sad in sad_list:
		
def get_num_comps(sad_list):
	des_set_list=[]
	uf_struct=UF(len(sad_list))
	sad_count=0
	for sad in sad_list:
		des_set=set()
		for id,k in msc.des(sad):
			des_set.add(id)
		for i in range(len(des_set_list)):
			intersect_set=des_set.intersection(des_set_list[i])
			if(len(intersect_set)!=0):
				for one_sad in intersect_set:
					if(msc.cp_func(one_sad)>=0):
						uf_struct.union(i,sad_count)
		des_set_list.append(des_set)
		sad_count+=1
	return (uf_struct.count())	

def get_des_vol(max_id):
	des_geom=msc.des_geom(max_id)
	count=0
	for cube_id in des_geom:
		if(get_cube_max_val(dual_pts[cube_id])>0):
			count+=1
	return count
	
def dia_from_vol(vol):
	dia=(vol*3.0)/(4.0*math.pi)
	dia=dia**(1/3.0)
	return 2.0*dia

max_cps=msc.cps(3)


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
	curr_max_data.append(max_id)	
	max_data_array.append(curr_max_data)

avg_dia=0
num_spheres=0

for max_id in max_cps:
	if(msc.cp_func(max_id)>0):
		dia=dia_from_vol(get_des_vol(max_id))
		avg_dia+=dia
		num_spheres+=1

avg_dia=avg_dia/num_spheres

count=0
for key in max_conn_dict.keys():
	if(len(max_conn_dict[key])>1):
		if(count==0):
			print(key)
		count+=1	

for key in max_conn_dict.keys():
	if(len(max_conn_dict[key])>1):
		clean_saddles(key,max_conn_dict[key])
#np.savetxt('des_mfold'+str(m)+'.csv',mfold_pt_list,delimiter=',')

