import pyms3d
import vtk.util.numpy_support as nps
import numpy as np
import vtk
import os
import argparse

parser=argparse.ArgumentParser()

parser.add_argument('data_file',type=str,help='data file name')

parser.add_argument('dim',nargs=3,type=int,help='dimensions')

parser.add_argument('output_path',type=str,help='output file path')

parser.add_argument('pers',type=float,help='persistence threshold for simplification')

args=parser.parse_args()

data_file_name=args.data_file
Dim = tuple(args.dim)
percent_pers=args.pers

output_path_name=args.output_path

base_name=os.path.basename(data_file_name)
base_name=os.path.splitext(base_name)[0]

output_poly_data_file_name='extremum_graph_'+base_name + '_pers_'+str(percent_pers)+'.vtp'
msc_file_name='msc_'+ base_name + '_pers_'+str(percent_pers)


print(output_poly_data_file_name)

DataFile = data_file_name

msc_computed=0

if(msc_computed==0):
    print pyms3d.get_hw_info()
    msc = pyms3d.mscomplex()
    msc.compute_bin(DataFile,Dim)
    msc.simplify_pers(thresh=percent_pers,is_nrm=False)
    msc.collect_geom(dim=2,dir=1)
    msc.save(output_path_name+msc_file_name)

if(msc_computed==1):
    msc=pyms3d.mscomplex()
    msc.load(msc_file_name)
    msc.simplify_pers(thresh=percent_pers,is_nrm=1)
    msc.save(output_path_name+msc_file_name)


dp = msc.dual_points()
pa = vtk.vtkPoints()
pa.SetData(nps.numpy_to_vtk(dp,"Pts"))

cps_2sad = msc.cps(2)

# create a vtk CellArray for the line segments
ca = vtk.vtkCellArray()
ia = vtk.vtkIntArray()
ia.SetName("SaddleIndex")
fa = vtk.vtkFloatArray()
fa.SetName("SaddleVal")
for s in cps_2sad:
#remove 2-sads not connected to 2 max
    if(len(msc.asc(s))!=2):
        continue
#change end
    gm = msc.asc_geom(s)
    for a,b in gm:
	print("Adding arc to ext graph")
        ca.InsertNextCell(2)
        ca.InsertCellPoint(a)
        ca.InsertCellPoint(b)
        ia.InsertNextValue(s)
        fa.InsertNextValue(msc.cp_func(s))

print(ca)


#Set the outputs
#pd = self.GetOutput()
pd=vtk.vtkPolyData()
pd.SetPoints(pa)
pd.SetLines(ca)
pd.GetCellData().AddArray(ia)
pd.GetCellData().AddArray(fa)

print(pd)

writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(output_path_name+output_poly_data_file_name)
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(pd)
else:
    writer.SetInputData(pd)
writer.Write()
