import numpy as np
from collections import defaultdict
import itertools

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



for max_ids,sads in max_conn_dict.iteritems():
    if(len(sads)>1):
            max_dist=0
            final_dist=[0,0,0]
            sad_pairs=itertools.combinations(sads,2)
            for pair in sad_pairs:
                    dist=np.array(msc.cp_cellid(pair[1]))-np.array(msc.cp_cellid(pair[0]))
                    if(np.linalg.norm(dist)>max_dist):
                            max_dist=np.linalg.norm(dist)
                            final_dist=dist
            print(final_dist)
            print(max_dist)


def get_dual_max_val(dual_pt):
	inc_array=[[0.5,0.5,0.5],[-0.5,0.5,0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5],[-0.5,-0.5,0.5],[0.5,-0.5,-0.5],[-0.5,0.5,-0.5],[-0.5,-0.5,-0.5]]
	max_val=-100000
	for i in inc_array:
		curr_pt=np.add(dual_pt,i)
		#print(curr_pt)
		curr_val=msc.vert_func(int(curr_pt[0]),int(curr_pt[1]),int(curr_pt[2]))
		#print(curr_val)
		if(curr_val>max_val):
			max_val=curr_val
	return max_val
		

def get_saddle_cofacets(saddle_cell_id):
	saddle_pt=[saddle_cell_id[0]/2.0,saddle_cell_id[1]/2.0,saddle_cell_id[2]/2.0]
	cofac_arr=[]
	arr_inc=[[0.5,0,0],[0,0.5,0],[0,0,0.5]]
	for i in range(3):
		if(saddle_pt[i]%1==0):
			#print(str(i) +' directional')
			cofac_arr.append(np.add(saddle_pt,arr_inc[i]))
			cofac_arr.append(np.add(saddle_pt,[-x for x in arr_inc[i]]))
	return(cofac_arr)

count=0
only_flat=0
for max_ids,sads in max_conn_dict.iteritems():
	flag=1
	for sad in sads:
		sad_val=msc.cp_func(sad)
            	cofacs=get_saddle_cofacets(msc.cp_cellid(sad))
            	cofac_val_one=get_dual_max_val(cofacs[0])
	    	cofac_val_two=get_dual_max_val(cofacs[1])
		if(cofac_val_one==cofac_val_two):
			if(cofac_val_one==sad_val):
				count+=1
			else:
				flag=0
		else:
			flag=0
	if(flag==1 and len(sads)>1):
		only_flat+=1
		print([msc.cp_func(x) for x in sads])
		print(sads)
				



init_set=get_saddle_cofacets(msc.cp_cellid(10977))
bd_set=[]
nbd_array=[[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]]
visited_dict={}

while len(init_set)!=0:
	curr_item=init_set.pop()
	visited_dict[tuple(curr_item)]=1
	for inc in nbd_array:
		test_item=np.add(curr_item,inc)
		if(tuple(test_item) in visited_dict.keys()):
			continue
		if(tuple(test_item) not in m_fold_dict.keys()):
			continue	
		flag=0
		for incr in nbd_array:
			rel_item=np.add(test_item,incr)
			if(tuple(rel_item) in m_fold_dict.keys()):
				if(m_fold_dict[tuple(rel_item)]!=m_fold_dict[tuple(test_item)]):
					flag=1
		if(flag==1):
			init_set.append(test_item)
	bd_set.append(curr_item)


