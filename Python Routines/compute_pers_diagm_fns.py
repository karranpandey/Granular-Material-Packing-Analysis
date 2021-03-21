import pyms3d
import vtk.util.numpy_support as nps
import numpy as np
import vtk
import os
import argparse
import matplotlib.pyplot as plt

def display_pcurve(p_diagm):
    sorted_arr =sorted(p_diagm, key = lambda x:x[2])
    pers_vals = []
    num_remaining = []
    list_len = len(sorted_arr)
    for i in range(len(sorted_arr)):
        pers_vals.append(sorted_arr[i][2])
        num_remaining.append(list_len - i)
    plt.plot(pers_vals,num_remaining)
    plt.show()
    return
        
def compute_pers_diagm(msc):
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
    display_pcurve(p_diagm_list)
    return

#np.savetxt(output_path_name+'pers_diagm_data_'+base_name,p_diagm_list)
#np.savetxt(output_path_name+'pers_diagm_csv_'+base_name+'.csv',p_diagm_list,delimiter=',')
