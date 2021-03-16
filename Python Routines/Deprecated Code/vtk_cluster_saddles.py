import pyms3d
import vtk
import vtk.util.numpy_support as nps
import SimpleITK as sitk

def write_polydata(pd, file_name):
    writer = vtk.vtkXMLPolyDataWriter()
    output_poly_data_file_name=file_name
    writer.SetFileName(output_path_name+output_poly_data_file_name)
    writer.SetInputData(pd)
    writer.Write()
    return

def addConnectivityData(dataset):
    connectivity_filter=vtk.vtkConnectivityFilter()
    connectivity_filter.SetInputData(dataset)
    connectivity_filter.SetExtractionModeToAllRegions()
    connectivity_filter.ColorRegionsOn()
    connectivity_filter.Update()
    return (connectivity_filter.GetOutput())

def extract_surviving_sads(des_man):
    region_val_dict = {}
    region_cp_dict = {}
    num_cells = des_man.GetNumberOfCells()
    regions = des_man.GetCellData().GetArray('RegionId')    
    cp_ids = des_man.GetCellData().GetArray('CP ID')
    for i in range(num_cells):
        region_id = regions.GetTuple1(i)
        cp_id = cp_ids.GetTuple1(i)
        cp_val = msc.cp_func(int(cp_id))
        if(region_id not in region_val_dict.keys()):
            region_val_dict[region_id] = cp_val
            region_cp_dict[region_id] = cp_id
        else:
            if(cp_val>region_val_dict[region_id]):
                region_val_dict[region_id] = cp_val
                region_cp_dict[region_id] = cp_id
    return list(set(region_cp_dict.values()))

msc = pyms3d.mscomplex()
msc.load('../Outputs/msc_chamf_distance_Sand_particles_pers_0.7')
output_path_name = '../Outputs/'

reader=sitk.ImageFileReader()
reader.SetImageIO("MetaImageIO")
reader.SetFileName('../Outputs/chamf_distance_Sand_particles.mhd')
image=reader.Execute()
image = sitk.GetArrayFromImage(image)
image = image.flatten()

msc.collect_geom(dim = 2, dir = 0)

primal_pts = msc.primal_points()

cps_2sad = msc.cps(2)

cp_ids = vtk.vtkIntArray()
val = vtk.vtkDoubleArray()
cp_ids.SetName('CP ID')
val.SetName('Val')

print('Data setting')

des_man_pts = vtk.vtkPoints()
des_man_quads = vtk.vtkCellArray()
des_man_pts.SetData(nps.numpy_to_vtk(primal_pts,"Pts"))

val = nps.numpy_to_vtk(image)
val.SetName('Distance Val')

print('init process start')

for s in cps_2sad:
    if(msc.cp_func(s) < 0 or len(msc.asc(s))!=2):
        continue
    des_man = msc.des_geom(s)
    for elem in des_man:
        if(min([image[elem[0]],image[elem[1]], image[elem[2]], image[elem[3]]])<0):
            continue
        cp_ids.InsertNextValue(int(s))
        des_man_quads.InsertNextCell(4)
        des_man_quads.InsertCellPoint(elem[0])
        des_man_quads.InsertCellPoint(elem[1])
        des_man_quads.InsertCellPoint(elem[3])
        des_man_quads.InsertCellPoint(elem[2])

des_man = vtk.vtkPolyData()
des_man.SetPoints(des_man_pts)
des_man.SetPolys(des_man_quads)
des_man.GetCellData().AddArray(cp_ids)
des_man.GetPointData().AddArray(val)

des_man = addConnectivityData(des_man)

surv_sads = extract_surviving_sads(des_man)

print(len(surv_sads))

write_polydata(des_man,'des_man.vtp')




    
