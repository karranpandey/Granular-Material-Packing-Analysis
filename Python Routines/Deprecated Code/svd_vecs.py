import numpy as np
from collections import defaultdict

def get_normal_vec(arr):
    G=arr.sum(axis=0)/arr.shape[0]
    u,s,vh = np.linalg.svd(arr - G)
    u_norm = vh[2,:]
    return u_norm

sad_surf_dict=defaultdict(list)
int_surfs=np.loadtxt('intersection_surfaces_2')
for elem in int_surfs:
    pt=[elem[0],elem[1],elem[2]]
    sad_surf_dict[elem[3]].append(pt)

vec_arr=[]

for key in sad_surf_dict.keys():
    normal_vec=get_normal_vec(sad_surf_dict[key])
    normal_vec.append(key)
    vec_arr.append(normal_vec)
