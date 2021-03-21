import pyms3d
import numpy as np
import vtk
import os 
import argparse
import vtk.util.numpy_support as nps
import itertools


def get_saddles(msc,surv_sads):
    cpList = []
    cps_2sad = msc.cps(2)
    for s in surv_sads:
        s = int(s)
        maxList = []
        co_ords = msc.cp_cellid(s)
        index = 2
        val = msc.cp_func(s)    
        for m,k in msc.asc(s):
            maxList.append(m)
        cpList.append((co_ords, index, val, s, maxList[0], maxList[1]))    
    cpList = list(set(cpList))        
    
    # create the vtk objects
    pa = vtk.vtkPoints()
    ia = vtk.vtkIntArray()
    max_1 = vtk.vtkIntArray()
    max_2 = vtk.vtkIntArray()
    ca = vtk.vtkCellArray()
    cp_ids = vtk.vtkIntArray()
    cp_ids.SetName("CP ID")
    ia.SetName("Index")
    fa = vtk.vtkFloatArray()
    fa.SetName("Val")
    max_1.SetName("Max 1")
    max_2.SetName("Max 2")
    for i,(p,idx,val,cpidx, m1idx, m2idx) in enumerate(cpList):
        pa.InsertNextPoint(np.array(p,np.float32)/2)
        ia.InsertNextValue(idx)
        fa.InsertNextValue(val)
        max_1.InsertNextValue(m1idx)
        max_2.InsertNextValue(m2idx)
        cp_ids.InsertNextValue(cpidx)
        ca.InsertNextCell(1)
        ca.InsertCellPoint(i)

    #Set the outputs

    pd=vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetVerts(ca)
    pd.GetPointData().AddArray(ia)
    pd.GetPointData().AddArray(fa)
    pd.GetPointData().AddArray(cp_ids)
    pd.GetPointData().AddArray(max_1)
    pd.GetPointData().AddArray(max_2)
    return pd    

def get_max(msc):
    cps_max = msc.cps(3)
    cpList = []
    for m in cps_max:
        if(msc.cp_func(m)>=0):
            cpList.append((msc.cp_cellid(m),3,msc.cp_func(m),m))    
    cpList = list(set(cpList))        
    
    # create the vtk objects
    pa = vtk.vtkPoints()
    ia = vtk.vtkIntArray()
    ca = vtk.vtkCellArray()
    cp_ids = vtk.vtkIntArray()
    cp_ids.SetName("CP ID")
    ia.SetName("Index")
    fa = vtk.vtkFloatArray()
    fa.SetName("Val")

    for i,(p,idx,val,cpidx) in enumerate(cpList):
        pa.InsertNextPoint(np.array(p,np.float32)/2)
        ia.InsertNextValue(idx)
        fa.InsertNextValue(val)
        cp_ids.InsertNextValue(cpidx)
        ca.InsertNextCell(1)
        ca.InsertCellPoint(i)

    #Set the outputs

    pd=vtk.vtkPolyData()

    pd.SetPoints(pa)
    pd.SetVerts(ca)
    pd.GetPointData().AddArray(ia)
    pd.GetPointData().AddArray(fa)
    pd.GetPointData().AddArray(cp_ids)
    return pd
    
    
def get_extremum_graph(msc,surviving_sads):
    msc.collect_geom(dim=2,dir=1)
    dp = msc.dual_points()
    pa = vtk.vtkPoints()
    pa.SetData(nps.numpy_to_vtk(dp,"Pts"))
    cps_2sad = msc.cps(2)
    ca = vtk.vtkCellArray()
    ia = vtk.vtkIntArray()
    ia.SetName("SaddleIndex")
    fa = vtk.vtkFloatArray()
    fa.SetName("SaddleVal")
    for s in surviving_sads:
        s = int(s)
        gm = msc.asc_geom(s)
        for a,b in gm:
            ca.InsertNextCell(2)
            ca.InsertCellPoint(a)
            ca.InsertCellPoint(b)
            ia.InsertNextValue(s)
            fa.InsertNextValue(msc.cp_func(s))
    pd=vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetLines(ca)
    pd.GetCellData().AddArray(ia)
    pd.GetCellData().AddArray(fa)
    return pd
    
    
def get_segmentation_primal(msc):
    msc.collect_geom(dim=3, dir = 0)
    dp = msc.dual_points()
    cps_max=msc.cps(3)
    inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
    pa = vtk.vtkPoints()
    ca = vtk.vtkCellArray()
    cp_ids = vtk.vtkIntArray()
    cp_ids.SetName("CP ID")
    count = 0
    for m in cps_max:
        if(msc.cp_func(m)<0):
            continue     
        des_geom=msc.des_geom(m)
        for cube_id in des_geom:
            dual_pt = dp[cube_id]
            for inc in inc_arr:
                curr_pt=np.add(dual_pt,inc)
                if(msc.vert_func(int(curr_pt[0]),int(curr_pt[1]), int(curr_pt[2]))<0):
                    continue
                pa.InsertNextPoint(curr_pt)
                ca.InsertNextCell(1)
                ca.InsertCellPoint(count)
                cp_ids.InsertNextValue(m)
                count+=1
    pd = vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetVerts(ca)
    pd.GetPointData().AddArray(cp_ids)
    return pd
         
def get_segmentation_dual(msc):
    msc.collect_geom(dim=3, dir = 0)
    dp = msc.dual_points()
    cps_max=msc.cps(3)
    inc_arr=list(itertools.product([0.5,-0.5],repeat=3))
    pa = vtk.vtkPoints()
    ca = vtk.vtkCellArray()
    cp_ids = vtk.vtkIntArray()
    cp_ids.SetName("CP ID")
    count = 0
    for m in cps_max:
        if(msc.cp_func(m)<0):
            continue     
        des_geom=msc.des_geom(m)
        for cube_id in des_geom:
            dual_pt = dp[cube_id]
            pa.InsertNextPoint(dual_pt)
            ca.InsertNextCell(1)
            ca.InsertCellPoint(count)
            cp_ids.InsertNextValue(m)
            count+=1
    pd = vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetVerts(ca)
    pd.GetPointData().AddArray(cp_ids)
    return pd
                
def get_segmentation_index_dual(msc, img):
    msc.collect_geom(dim=3, dir = 0)
    dp = msc.dual_points()
    pp = msc.primal_points()
    cps_max=msc.cps(3)
    pa = vtk.vtkPoints()
    ca = vtk.vtkCellArray()
    cp_ids = vtk.vtkIntArray()
    cp_ids.SetName("CP ID")
    #val = vtk.vtkDoubleArray()
    #val.SetName('Distance Val')
    count = 0
    for m in cps_max:
        if(msc.cp_func(m)<0):
            continue     
        des_geom=msc.des_geom(m)
        for cube_id in des_geom:
            dual_pt = dp[cube_id]
            #val.InsertNextValue(img[int(dual_pt[0]),int(dual_pt[1]),int(dual_pt[2])])
            if(img[int(dual_pt[0]),int(dual_pt[1]),int(dual_pt[2])] < 0):
                continue
            #if(img[cube_id] < 0): 
            #    continue
            pa.InsertNextPoint(dual_pt)
            ca.InsertNextCell(1)
            ca.InsertCellPoint(count)
            cp_ids.InsertNextValue(m)
            count+=1
    pd = vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetVerts(ca)
    pd.GetPointData().AddArray(cp_ids)
    #pd.GetPointData().AddArray(val)
    return pd
    
    

