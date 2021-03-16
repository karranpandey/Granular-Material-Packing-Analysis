import pyms3d
import numpy as np
import os

msc=pyms3d.mscomplex()
msc_file_name='test2/msc_level_set_chamf_distance_testData2_pers_0.0'
base_name=os.path.basename(msc_file_name)
base_name=os.path.splitext(base_name)[0]
msc.load(msc_file_name)
cps_2sad=msc.cps(2)
cps_max=msc.cps(3)
cps_fun_vals=msc.cps_func()
msc.simplify_pers(thresh=1,is_nrm=True)
cp_pairs=msc.cps_pairid()

p_diagm_list=[]

for m in cps_max:
    sad=cp_pairs[m]
    b_value=cps_fun_vals[sad]
    d_value=cps_fun_vals[m]
    pers=d_value-b_value
    p_diagm_list.append([b_value,d_value,pers,sad,m])

np.savetxt('pers_diagm_data_'+base_name,p_diagm_list)
np.savetxt('pers_diagm_csv_'+base_name+'.csv',p_diagm_list,delimiter=',')
