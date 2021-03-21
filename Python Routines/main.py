import pyms3d
import numpy as np
import vtk
import os
import argparse
import SimpleITK as sitk
from scipy.ndimage.filters import convolve
from scipy.ndimage import maximum_filter, minimum_filter
from vtk_cluster_saddle_fns import compute_contact_regions
from crit_point_fns import get_saddles, get_max, get_extremum_graph, get_segmentation_index_dual
from compute_pers_diagm_fns import compute_pers_diagm

def write_polydata(pd, file_name):
    writer = vtk.vtkXMLPolyDataWriter()
    output_poly_data_file_name=file_name
    writer.SetFileName(output_path_name+output_poly_data_file_name)
    writer.SetInputData(pd)
    writer.Write()
    return

def read_img(img_file_name):
    reader=sitk.ImageFileReader()
    reader.SetImageIO("MetaImageIO")
    reader.SetFileName(img_file_name)
    image=reader.Execute()
    image = sitk.GetArrayFromImage(image)
    return image

def read_msc_to_img(msc, dim):
    np_arr = np.zeros(dim)
    for x in range(dim[0]):
        for y in range(dim[1]):
            for z in range(dim[2]):
                np_arr[x,y,z] = msc.vert_func(x,y,z)
    return np_arr #flatten order - F
   
parser=argparse.ArgumentParser()

parser.add_argument('data_file',type=str,help='data file name')

parser.add_argument('dim',nargs=3,type=int,help='dimensions')

args=parser.parse_args()

data_file_name=args.data_file
dim = tuple(args.dim)

output_path_name='../Outputs/'

base_name=os.path.basename(data_file_name)
base_name=os.path.splitext(base_name)[0]

msc_file_name='msc_'+ base_name + '_unsimplified'

#print(pyms3d.get_hw_info())

print("Computing initial morse-smale complex")

msc = pyms3d.mscomplex()
msc.compute_bin(data_file_name,dim)
msc.simplify_pers(thresh=0.0,is_nrm=True)

img = read_msc_to_img(msc, dim)
raw_writer = sitk.ImageFileWriter()
output_path_name='../Outputs/'
raw_writer.SetFileName(output_path_name+'msc_image.mhd')
raw_writer.Execute(sitk.GetImageFromArray(img))


while(True):
    val = input("1. Display Persistence Curve \n 2. Simplify Morse-Smale Complex \n 3. Compute Contact Information \n 4. Compute Segmentation \n 5. Load Stored MSC \n 6. Exit \n")
    val = int(val)
    if(val == 1):
        compute_pers_diagm(msc)
    if(val == 2):
        percent_pers = float(input("Enter Persistence Threshold: \n"))
        msc_file_name+='_pers_'+str(percent_pers)
        msc.simplify_pers(thresh=percent_pers,is_nrm=False)
        msc.save(output_path_name+msc_file_name)
        print('MSC Simplified')
    if(val == 3):
        print('Computing Contact Regions')
        des_man, surv_sads = compute_contact_regions(msc, img)
        print('Contact Regions Extracted')
        contacts = get_saddles(msc,surv_sads)
        print('Contacts Computed')
        grain_centres = get_max(msc)
        print('Grain Centres Computed')
        connectivity_network = get_extremum_graph(msc,surv_sads)
        print('Connectivity Network Computed')
        write_polydata(grain_centres,base_name + 'grain_centres.vtp')
        write_polydata(contacts,base_name + 'contacts.vtp')
        write_polydata(des_man,base_name + 'contact_regions.vtp')
        write_polydata(connectivity_network,base_name + 'connectivity_network.vtp')
    if(val == 4):
        #weights = [[[1/8,1/8],[1/8,1/8]],[[1/8,1/8],[1/8,1/8]]]
        #dual_img = convolve(img, weights, mode = 'constant')
        foot_print = [[[0,0,0],[0,0,0],[0,0,0]],[[0,0,0],[0,1,1],[0,1,1]],[[0,0,0],[0,1,1],[0,1,1]]]
        min_img = minimum_filter(img, footprint = foot_print)
        segmentation = get_segmentation_index_dual(msc, min_img)
        print('Segmentation Computed')
        write_polydata(segmentation, base_name + 'segmentation.vtp')
    if(val == 5):
        file_name = input("Enter MSC Path")
        msc=pyms3d.mscomplex()
        msc.load(file_name)    
    if(val == 6):
        break



