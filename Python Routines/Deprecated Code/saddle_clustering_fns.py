import pyms3d
import itertools
import numpy as np
from collections import defaultdict
from python_algorithms.basic.union_find import UF
import argparse
import vtk.util.numpy_support as nps
import vtk


def check_bd_cell_id(cell_id,mfold_1,mfold_2,min_sad_val):
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

def get_cube_max_val(dual_pt):
    inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
    max_val=-1
    #flag=0
    for inc in inc_arr:
        curr_pt=np.add(dual_pt,inc)
        if(curr_pt[0]<x_max and curr_pt[1]<y_max and curr_pt[2]<z_max):
            #flag=1
            curr_val=msc.vert_func(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
            if(curr_val>max_val):
                max_val=curr_val
    #print(flag)
    return(max_val)


def get_nbs_cell_id(cell_id):
    nbd_list=[]
    nbd_arr=list(itertools.product([1,-1,0],repeat=3))
    nbd_arr.remove((0,0,0))
    for arr in nbd_arr:
        inc=(arr[0]+x_max*arr[1]+y_max*x_max*arr[2])
        if((cell_id+inc)>0 and (cell_id+inc)<x_max*y_max*z_max):
            nbd_list.append(cell_id+inc)
    return nbd_list

def get_sad_facs(sad_id):
    sad_pt=[x/2.0 for x in msc.cp_cellid(sad_id)]
    fac_arr=[]
    poss_coords=defaultdict(list)
    for i in range(3):
        if(sad_pt[i]%1!=0):
            poss_coords[i].append(sad_pt[i]+0.5)
            poss_coords[i].append(sad_pt[i]-0.5)
        else:
            poss_coords[i].append(sad_pt[i])
    fac_arr=list(itertools.product(poss_coords[0],poss_coords[1],poss_coords[2]))
    return(fac_arr)


def check_sad_inside(sad_id):
    avg_val=0
    for fac in get_sad_facs(sad_id):
        avg_val+=msc.vert_func(int(fac[0]),int(fac[1]),int(fac[2]))
        #if(msc.vert_func(int(fac[0]),int(fac[1]),int(fac[2]))<0):
        #    return False
    if(avg_val<0):
        return False
    return True


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


def bfs_intersection_ids(m1,m2,uf_id_dict):
    visited=[]
    uf_struct=UF(0)
    init_list=[]
    min_sad_val=100
    for sad in max_conn_dict[(m1,m2)]:
        init_list.extend(get_saddle_cofac_cell_id(sad))
        if(msc.cp_func(sad)<min_sad_val):
            min_sad_val=msc.cp_func(sad)
    init_uf_sads(uf_struct,init_list,uf_id_dict)
    queue=init_list
    bd_int=[]
    explored=set()
    mfold_1=msc.des_geom(m1)
    mfold_2=msc.des_geom(m2)
    while queue:
        if(uf_struct.count()==1):
            return uf_struct
        node=queue.pop(0)
        if node not in explored:
            explored.add(node)
            neighbours=get_nbs_cell_id(node)
            for neighbour in neighbours:
                if(check_bd_cell_id(neighbour,mfold_1,mfold_2,min_sad_val)):
                    queue.append(neighbour)
                    uf_add_node(node,neighbour,uf_struct,uf_id_dict)
    #print(uf_struct.count())
    return uf_struct

def get_surviving_sads(m1,m2):
    uf_id_dict={}
    uf_struct=bfs_intersection_ids(m1,m2,uf_id_dict)
    sad_list=max_conn_dict[(m1,m2)]
    sad_comp_dict=defaultdict(list)
    surviving_sads=[]
    for sad in sad_list:
        comp_id=(uf_struct.find(uf_id_dict[get_saddle_cofac_cell_id(sad)[0]]))
        sad_comp_dict[comp_id].append(sad)
    for key in sad_comp_dict.keys():
        max_sad_val=0
        max_sad_index=0
        for sad in sad_comp_dict[key]:
            if(msc.cp_func(sad)>max_sad_val):
                max_sad_val=msc.cp_func(sad)
                max_sad_index=sad
        surviving_sads.append(max_sad_index)
    return surviving_sads


def cluster_saddles(msc_init, dim):
    global msc, max_conn_dict,x_max,y_max,z_max, dual_pts, max_list, max_pairs, sad_list, surviving_sads
    x_max=dim[0]-1
    y_max=dim[1]-1
    z_max=dim[2]-1
    msc = msc_init
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

    for key in max_conn_dict.keys():
        sad_list=max_conn_dict[key]
        #print(len(max_conn_dict[key]))
        for sad in sad_list:
            if(check_sad_inside(sad)==False):
                max_conn_dict[key].remove(sad)
        #print(len(max_conn_dict[key]))

    for key in max_conn_dict.keys():
        sad_list=max_conn_dict[key]
        if(len(sad_list)>1):
            sad_list=get_surviving_sads(key[0],key[1])
        max_conn_dict[key]=sad_list

    surviving_sads=set()
    for key in max_conn_dict.keys():
        for sad in max_conn_dict[key]:
            surviving_sads.add(sad)
    np.savetxt('surviving_sads_1',list(surviving_sads))
    return list(surviving_sads)


#np.savetxt(output_path_name+'surviving_sads',list(surviving_sads))

#write_extremum_graph(msc,surviving_sads)

