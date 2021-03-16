import pyms3d
import numpy as np
from collections import defaultdict
import itertools

msc_file_name='msc_otsu_Depositional_Steel_Air_downsampled_mean_213_256_251_pers_2.0'
msc=pyms3d.mscomplex()
msc.load(msc_file_name)
max_conn_dict=defaultdict(list)
cps_2sad=msc.cps(2)
for s in cps_2sad:
    if(msc.cp_func(s)<0):
    	    continue
    max_list=[]
    for m,k in msc.asc(s):
            max_list.append(m)
    max_pairs=list(itertools.combinations(max_list,2))
    for pair in max_pairs:
            max_conn_dict[tuple(pair)].append(s)

dist_tuples={}
	
for max_ids,sads in max_conn_dict.iteritems():
    if(len(sads)>1):
            sad_pairs=itertools.combinations(sads,2)
            max_dist=0
            for pair in sad_pairs:
                    dist=np.linalg.norm(np.array(msc.cp_cellid(pair[1]))-np.array(msc.cp_cellid(pair[0])))
                    if(dist>max_dist):
                            max_dist=dist
            print(max_dist)
            dist_tuples[max_ids]=(max_dist)


