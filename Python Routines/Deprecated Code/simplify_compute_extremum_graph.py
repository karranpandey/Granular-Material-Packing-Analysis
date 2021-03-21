import pyms3d
import numpy as np
import vtk
import os
import argparse
import SimpleITK as sitk
from scipy.ndimage.filters import convolve
from scipy.ndimage import maximum_filter
from vtk_cluster_saddle_fns import compute_contact_regions
from crit_point_fns import get_saddles, get_max, get_extremum_graph, get_segmentation_index_dual

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

#parser.add_argument('data_file',type=str,help='data file name')

parser.add_argument('msc_file', type=str, help = 'msc file name')

parser.add_argument('dim',nargs=3,type=int,help='dimensions')

#parser.add_argument('output_path',type=str,help='output file path')

parser.add_argument('--pers',type=float,help='persistence threshold for simplification')

parser.add_argument('--contacts',type=float,help='persistence threshold for simplification')

parser.add_argument('--network',type=float,help='persistence threshold for simplification')

parser.add_argument('--segmentation',type=float,help='persistence threshold for simplification')

args=parser.parse_args()

#img_file_name = args.data_file
msc_file_name=args.msc_file
dim = tuple(args.dim)
percent_pers=args.pers
contacts = args.contacts
network = args.network
seg = args.segmentation

output_path_name='../Outputs/'

base_name=os.path.basename(msc_file_name)
base_name=os.path.splitext(base_name)[0]


msc=pyms3d.mscomplex()
msc.load(msc_file_name)
img = read_msc_to_img(msc, dim)

if(percent_pers is not None):
    msc_file_name+=_file_name + '_pers_'+str(percent_pers)
    msc.simplify_pers(thresh=percent_pers,is_nrm=False)
    msc.save(output_path_name+msc_file_name)
    print('MSC Simplified')

if(contacts is not None):
    print('Clustering Saddles')
    des_man, surv_sads = compute_contact_regions(msc, img)
    print('Saddles Clustered, Contact Regions Extracted')
    contacts = get_saddles(msc,surv_sads)
    print('Contacts Written')
    grain_centres = get_max(msc)
    print('Grain Centres Written')
    write_polydata(grain_centres,base_name + 'grain_centres.vtp')
    write_polydata(contacts,base_name + 'contacts.vtp')
    write_polydata(des_man,base_name + 'contact_regions.vtp')

if(network is not None):
    connectivity_network = get_extremum_graph(msc,surv_sads)
    print('Connectivity Network Computed')
    write_polydata(connectivity_network,base_name + 'connectivity_network.vtp')

if(seg is not None):
    #weights = [[[1/8,1/8],[1/8,1/8]],[[1/8,1/8],[1/8,1/8]]]
    #dual_img = convolve(img, weights, mode = 'constant')
    max_img = maximum_filter(img, size = 2)
    segmentation = get_segmentation_index_dual(msc, max_img)
    print('Segmentation Computed')
    write_polydata(segmentation, base_name + 'segmentation.vtp')

