import pyms3d
import numpy as np
import vtk
import os 
import argparse

parser=argparse.ArgumentParser()

parser.add_argument('data_file',type=str,help='data file name')

parser.add_argument('pers',type=float,help='persistence threshold for simplification')

parser.add_argument('output_path',type=str,help='output file path')

args=parser.parse_args()

msc_file_name=args.data_file
percent_pers=args.pers
output_path = args.output_path
#base_name=os.path.basename(data_file_name)
#base_name=os.path.splitext(base_name)[0]

#msc_file_name='msc_'+ base_name + '_pers_'+str(percent_pers)
crit_points_file_name='crit_points_pers_'+str(percent_pers)+'.vtp'

msc=pyms3d.mscomplex()
msc.load(msc_file_name)
cps_2sad = msc.cps(2)

cpList = []
for s in cps_2sad:
#remove 2-sads not connected to 2 max
    if(len(msc.asc(s))!=2):
	continue
#end change
    cpList.append((msc.cp_cellid(s),2,msc.cp_func(s),s))
    for m,k in msc.asc(s):
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
pd.SetLines(ca)
pd.GetPointData().AddArray(ia)
pd.GetPointData().AddArray(fa)
pd.GetPointData().AddArray(cp_ids)

print(pd)

writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(output_path+crit_points_file_name)
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(pd)
else:
    writer.SetInputData(pd)
writer.Write()
