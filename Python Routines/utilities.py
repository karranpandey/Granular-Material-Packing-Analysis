# PYMS3d modules
import pyms3d

# Reading files: module
import h5py
from scipy import io

# Imaging modules
import itk
import SimpleITK as sitk
from skimage import filters
from skimage.measure import block_reduce
from skimage.segmentation import morphological_chan_vese

# Input output modules
import os
import re
import vtk
import vtk.util.numpy_support as nps

# Plotting modules
from matplotlib import pyplot as plt

# Other modules
import itertools
from numba import jit
import numba
import numpy as np


def addConnectivityData(dataset):
    """Extract cells that share common points

    Args:
        dataset (vtk polydata): vtk polydata with points and cells

    Returns:
        vtk polydata: polydata with connectivity
    """
    connectivity_filter = vtk.vtkConnectivityFilter()
    connectivity_filter.SetInputData(dataset)
    connectivity_filter.SetExtractionModeToAllRegions()
    connectivity_filter.ColorRegionsOn()
    connectivity_filter.Update()
    return (connectivity_filter.GetOutput())


@jit(nopython=True)
def adaptive_thresh(image, T=0.15):
    """Adaptive thresholding

    Args:
        image (uint16 np array): two-dimensional image
        T (float, optional): background pixel is less than mean*(1-T).
        Defaults to 0.15.

    Reference:
    Bradley, D. and Roth, G., 2007. Adaptive thresholding using the integral
    image. Journal of graphics tools, 12(2), pp.13-21.

    Returns:
        out_image: thresholded image
    """
    num_row, num_col = image.shape
    win_size = max(num_row, num_col) // 8
    half_win_size = win_size // 2

    # initialize integral image and output image
    int_image = image.copy().astype(numba.uint64)
    out_image = np.zeros((num_row, num_col), dtype=numba.uint8)

    # integral image
    for col in range(num_col):
        temp = 0
        for row in range(num_row):
            temp = temp + image[row, col]
            if col == 0:
                int_image[row, col] = temp
            else:
                int_image[row, col] = int_image[row, col-1] + temp

    for col in range(num_col):
        for row in range(num_row):
            # win_size x win_size region around row and col
            irow = max(row-half_win_size, 0)
            frow = min(row+half_win_size, num_row-1)
            icol = max(col-half_win_size, 0)
            fcol = min(col+half_win_size, num_col-1)

            region_size = (frow-irow)*(fcol-icol)

            sum_ = int_image[frow, fcol]-int_image[irow, fcol] - \
                int_image[frow, icol]+int_image[irow, icol]

            if image[row, col]*region_size <= sum_*(1.0-T):
                out_image[row, col] = 0
            else:
                out_image[row, col] = 255
    return out_image


def bd_extraction(arr, slicewise=True, filterName='InterMode', ace=False,
                  visualize=False):
    """Boundary extraction from image using different filters

    Args:
        arr (np array): attenuation density image
        slicewise (bool, optional): slicewise if true; globally if false. Defaults to True.
        filterName (str, optional): many thresholding filters are available. Defaults to 'InterMode'.
        ace (bool, optional): active countour as boundary surface if true. Defaults to False.
        visualize (bool, optional): visualize the result if true. Defaults to False.

    Returns:
        numpy array: Binary volume after filtering and active contour
    """
    print('Commencing Boundary Extraction')
    # dictionary of filters
    # Minimum filter is from itk
    # Adaptive filter is manually implemented
    filterTypeList = dict(
        Huang=sitk.HuangThresholdImageFilter,
        Isodata=sitk.IsoDataThresholdImageFilter,
        InterMode=sitk.IntermodesThresholdImageFilter,
        Kitt=sitk.KittlerIllingworthThresholdImageFilter,
        Li=sitk.LiThresholdImageFilter,
        MaxEnt=sitk.MaximumEntropyThresholdImageFilter,
        Min=itk.IntermodesThresholdImageFilter,
        Moments=sitk.MomentsThresholdImageFilter,
        Otsu=sitk.OtsuThresholdImageFilter,
        Renyi=sitk.RenyiEntropyThresholdImageFilter,
        Shanbhag=sitk.ShanbhagThresholdImageFilter,
        Triangle=sitk.TriangleThresholdImageFilter,
        Yen=sitk.YenThresholdImageFilter,
        Adaptive=adaptive_thresh,
    )
    filterType = filterTypeList[filterName]

    # slicewise thresholding
    if slicewise:
        ls = np.zeros_like(arr, dtype=np.float32)
        Thresh = np.zeros(arr.shape[0])
        for ii, slice in enumerate(arr):

            if np.max(slice) == 0:
                continue

            # thresholding
            if filterName == 'Min':
                # minimum filter
                filter = filterType.New(itk.GetImageFromArray(slice))
                filter.SetInput(itk.GetImageFromArray(slice))
                filter.UseInterModeOff()
                filter.Update()
                thresh = filter.GetThreshold()

            elif filterName == 'Adaptive':
                # adaptive thresholding
                slice = filterType(slice, T=0.10)
                thresh = 100

            else:
                # sitk filter
                filter = filterType()
                filter.Execute(sitk.GetImageFromArray(slice))
                thresh = filter.GetThreshold()

            Thresh[ii] = thresh

            # active contour
            if ace:
                # active contour
                lambda_1, lambda_2 = 1, 1
                active_contour = morphological_chan_vese(slice, 200, init_level_set=(slice-thresh),
                                                         lambda1=abs(lambda_1), lambda2=abs(lambda_2), smoothing=0)
                ls[ii, :, :] = (active_contour).astype(np.float32)
            else:
                ls[ii, :, :] = (slice > thresh).astype(np.float32)

            if visualize:
                fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(18, 9))
                ax0.imshow(arr[ii, :, :], cmap="Greys")
                ax1.imshow(slice > thresh, cmap="Greys")
                plt.savefig(f'./Dump/slice_{ii}.png', dpi=600)
                plt.close(fig)

        # fill enclosed voids in the binary volume
        fill = sitk.BinaryFillholeImageFilter()
        ls = fill.Execute(sitk.GetImageFromArray(ls.astype(np.ubyte)))
        ls = sitk.GetArrayFromImage(ls)
        np.savetxt('Threshold.txt', Thresh, fmt='%10.4f')

    else:
        # thresholding
        if filterName == 'Min':
            filter = itk.IntermodesThresholdImageFilter.New(
                itk.GetImageFromArray(arr))
            filter.SetInput(itk.GetImageFromArray(arr))
            filter.UseInterModeOff()
            filter.Update()
            thresh = filter.GetThreshold()

        elif filterName == 'Adaptive':
            print('Adaptive thresholding is not available for 3d')
            raise NotImplementedError

        else:
            filter = filterType()
            filter.Execute(sitk.GetImageFromArray(arr))
            thresh = filter.GetThreshold()
            print("The threshold is: ", thresh)

        # active contour
        if ace:
            # active contour
            lambda_1, lambda_2 = 1, 1
            active_contour = (morphological_chan_vese(arr, 200,
                                                      init_level_set=(
                                                          arr-thresh),
                                                      lambda1=abs(lambda_1), lambda2=abs(lambda_2), smoothing=0))
            ls = (active_contour).astype(np.float32)
        else:
            ls = (arr > thresh).astype(np.float32)

    print('Boundary Extracted.')
    return ls


def bimode_log_min(vols):
    """Automatic thresholding -- volume cuoff

    Args:
        vols (list): list of volume of particles

    Returns:
        float: volume cutoff
    """
    plt.figure()
    plt.hist(np.log(vols))
    plt.show()
    return np.exp(filters.threshold_otsu(np.log(vols)))


def check_segmentation(raw_file_name, centers):
    """Write the centers on slices

    Args:
        raw_file_name (str): raw file with original data
        centers (list): list of particle centers

    Returns:
        None: None
    """
    reader = sitk.ImageFileReader()
    reader.SetImageIO("MetaImageIO")
    reader.SetFileName(raw_file_name)
    arr = sitk.GetArrayFromImage(reader.Execute())
    arr = arr.transpose(2, 1, 0)
    # arrmin, arrmax = np.min(arr), np.max(arr)

    s = 2
    for center in centers:
        center = center.astype(int)
        try:
            arr[center[0]-s:center[0]+s, center[1]-s:center[1] +
                s, center[2]-s:center[2]+s] = 0
        except:
            arr[center[0]-s:center[0]+s, center[1] -
                s:center[1]+s, center[2]] = 0

    for ii in range(arr.shape[2]):
        slice = arr[:, :, ii]
        # slice = arrmin + (arrmax - arrmin) * slice / (np.max(slice) - np.min(slice))
        fig, ax = plt.subplots(1, 1)
        ax.imshow(slice, interpolation="none", cmap="Greys_r")
        plt.savefig(f'./Dump/check_slice_{ii:04d}')
        plt.close(fig)

    return None


def collect_neighbours(points):
    """This function collects the 8-neighbours for midpoints in a grid

    Args:
        points (np array): array of midpoints in a grid

    Returns:
        np array: 8-neighbour grid points
    """
    inc = np.array(list(itertools.product((-0.5, 0.5), repeat=3)))
    pts = inc.reshape(1, inc.shape[1], inc.shape[0])\
        + points.reshape(points.shape[0], points.shape[1], 1)
    pts = pts.transpose([0, 2, 1]).reshape(-1, 3)
    return np.unique(pts.astype(int), axis=0)


def compute_contact_regions(msc, image):
    """Get the contact regions

    Args:
        msc (msc object): Morse Complex object
        image (np aray): distance field

    Returns:
        vtk polydata: vtkpolydata with contact points and cells
    """
    # get the contact region -- descending manifold of 2-saddle
    msc.collect_geom(dim=2, dir=0)
    # ''' critical points type
    # dim: Critical point type \n"\
    #     "   dim=-1      --> All (default)\n"\
    #     "   dim=0,1,2,3 --> Minima, 1-saddle,2-saddle,Maxima \n"\
    # dir: Geometry type \n"\
    #     "   dir=0 --> Descending \n"\
    #     "   dir=1 --> Ascending \n"\
    #     "   dir=2 --> Both (default) \n"\
    # '''

    # get the coordinates of primal points
    primal_pts = msc.primal_points()

    # get the 2 saddle points
    cps_2sad = msc.cps(2)

    # initialize vtk data type
    cp_ids = vtk.vtkIntArray()
    cp_ids.SetName('CP ID')
    val = vtk.vtkDoubleArray()
    val.SetName('Val')
    des_man_pts = vtk.vtkPoints()
    des_man_pts.SetData(nps.numpy_to_vtk(primal_pts, "Pts"))
    des_man_quads = vtk.vtkCellArray()

    surv_sads = []
    for s in cps_2sad:
        # ignore the saddles in background
        # or saddle which is connected to just one maxima
        if (msc.cp_func(s) < 0) or (len(msc.asc(s)) != 2):
            # print("Saddle belongs to background")
            continue

        # if (msc.cp_func(msc.asc(s)[0, 0]) < 0) or (msc.cp_func(msc.asc(s)[1, 0]) < 0):
        #     print("This maxima point lies in the background.")

        # Descending manifold geometry of 2-saddle point
        des_man = msc.des_geom(s)
        surv_sads.append(int(s))  # added inplace of exract surv saddle
        for elem in des_man:
            # check if elem in descending manifold belongs to background
            ind = np.ravel_multi_index(
                primal_pts[elem].astype(int).transpose(), image.shape)
            dist_vals = image.ravel()[ind]
            if (np.min(dist_vals) <= 0):
                continue

            # for a  correct elem stores the cp_id and des_man_quad
            cp_ids.InsertNextValue(int(s))
            des_man_quads.InsertNextCell(4)
            des_man_quads.InsertCellPoint(elem[0])
            des_man_quads.InsertCellPoint(elem[1])
            des_man_quads.InsertCellPoint(elem[3])
            des_man_quads.InsertCellPoint(elem[2])

    # stores the descending manifolds
    des_man = vtk.vtkPolyData()
    des_man.SetPoints(des_man_pts)
    des_man.SetPolys(des_man_quads)
    des_man.GetCellData().AddArray(cp_ids)
    # des_man.GetPointData().AddArray(val)
    des_man = addConnectivityData(des_man)
    # des_man, surv_sads = extract_surviving_sads(des_man, msc)
    return des_man, surv_sads

def compute_pers_diagm(data_file_name, dim):
    """Comput the persistence diagram

    Args:
        data_file_name (str): raw file with distance field
        dim (tuple): dimensions of distance field

    Returns:
        None: None
    """
    # Comput msc
    msc = pyms3d.mscomplex()
    msc.compute_bin(data_file_name, dim)
    # simplify for base case
    msc.simplify_pers(thresh=0.0, is_nrm=True)
    # get the critical points
    cps_max = msc.cps(3)
    # cps_min, cps_1sad, cps_2sad = msc.cps(0), msc.cps(1), msc.cps(2)
    # value of scalar function at the critical point
    cps_fun_vals = msc.cps_func()

    # simplify for highest persistence (normalized value = 1)
    msc.simplify_pers(thresh=1, is_nrm=True)
    # get the critical point pairs that got cancelled (all got cancelled)
    # saddle--maximum pairs cancelled in simplification
    cp_pairs = msc.cps_pairid()

    # persistence diagram between birth and death
    sad = cp_pairs[cps_max, np.newaxis]
    b_value = cps_fun_vals[sad]
    d_value = cps_fun_vals[cps_max, np.newaxis]
    pers = d_value - b_value
    # p_diagm_list = np.concatenate((b_value, d_value, pers, sad,
    #                               cps_max[:, np.newaxis]), axis=1)

    print('Saddle values (quantiles - 0.25 spacing) ',
          np.quantile(b_value, [0, 0.25, 0.50, 0.75, 1.0]))
    print('Max values (quantiles - 0.25 spacing) ',
          np.quantile(d_value, [0, 0.25, 0.50, 0.75, 1.0]))

    plt.figure()
    plt.hist(b_value, stacked=True, label="saddle value", rwidth=0.1,
             cumulative=True, density=True, histtype='bar')
    plt.hist(d_value, stacked=True, label="Max value", rwidth=0.1,
             cumulative=True, density=True, histtype='bar')
    plt.xlabel('func value')
    plt.legend()
    plt.show()

    print("Take a note of the persistence value at the knee!")
    plt.figure()
    plt.plot(np.sort(pers, axis=0)[::-1], np.arange(pers.shape[0]))
    plt.xlabel("Persistence")
    plt.ylabel("Survived critical points")
    plt.title("Persistence curve")
    plt.show()

    plt.figure()
    plt.plot(b_value, d_value, 'r.')
    plt.plot([0, max(max(b_value), max(d_value))],
             [0, max(max(b_value), max(d_value))])
    plt.xlabel("Birth")
    plt.ylabel("Death")
    plt.title("Persistence diagram")
    plt.show()
    return None


def dist_field_comp(ls):
    """Get the chamfer distance field from the binary volume

    Args:
        ls (numpy array): binary volume array

    Returns:
        numpy array: distance field as numpy array
    """
    print('Commencing Distance Field Computation')
    itk_image = itk.GetImageFromArray(ls.astype(np.float32))

    antialiasfilter = itk.AntiAliasBinaryImageFilter.New(itk_image)
    antialiasfilter.SetInput(itk_image)
    antialiasfilter.Update()
    antialias_image = antialiasfilter.GetOutput()

    isoContourFilter = itk.IsoContourDistanceImageFilter.New(antialias_image)
    isoContourFilter.SetLevelSetValue(0.5)
    isoContourFilter.SetFarValue(100)
    isoContourFilter.SetInput(antialias_image)
    isoContourFilter.Update()
    isoContour_image = isoContourFilter.GetOutput()

    chamferFilter = itk.FastChamferDistanceImageFilter.New(isoContour_image)
    chamferFilter.SetMaximumDistance(50.0)
    chamferFilter.Update()
    chamferFilter.SetInput(isoContour_image)
    chamferFilter.Update()
    chamf_image = chamferFilter.GetOutput()
    dist_field = itk.GetArrayFromImage(chamf_image)
    print('Distance Field Computed')
    return dist_field.astype(np.float32)


def extract_surviving_sads(des_man, msc):
    """Extract surviving saddles from descending manifold after contact computation

    Args:
        des_man (vtk polydata): contact regions
        msc (msc object): Morse Complex object

    Returns:
        list: list of surviving saddles indices
    """
    # region_val_dict = {}
    # region_cp_dict = {}
    # num_cells = des_man.GetNumberOfCells()
    # regions = des_man.GetCellData().GetArray('RegionId')
    # cp_ids = des_man.GetCellData().GetArray('CP ID')
    # for i in range(num_cells):
    #     region_id = regions.GetTuple1(i)
    #     cp_id = cp_ids.GetTuple1(i)
    #     cp_val = msc.cp_func(int(cp_id))
    #     if(region_id not in region_val_dict.keys()):
    #         region_val_dict[region_id] = cp_val
    #         region_cp_dict[region_id] = cp_id
    #     else:
    #         if(cp_val > region_val_dict[region_id]):
    #             region_val_dict[region_id] = cp_val
    #             region_cp_dict[region_id] = cp_id
    # for i in range(num_cells):
    #     region_id = regions.GetTuple1(i)
    #     cp_ids.SetTuple1(i, region_cp_dict[region_id])
    # return des_man, list(set(region_cp_dict.values()))
    pass


def get_dims(filename):
    """This function extracts the DimSize as tuple from the mhd file.

    Args:
        filename (str): mhd file with the dimensions of 3-d volume

    Returns:
        tuple: dimensions
    """
    path, file = os.path.split(filename)
    filename = os.path.join(path, os.path.splitext(file)[0]+'.mhd')
    with open(filename, mode='r') as f:
        text = f.read().split('\n')

    for line in text:
        if re.search('DimSize', line):
            txt = line.split(' ')
            dims = (int(txt[2]), int(txt[3]), int(txt[4]))
    return dims


def get_extremum_graph(msc, surviving_sads):
    """Get the extremum graph of surviving saddles

    Args:
        msc (msc object): Morse Smale object
        surviving_sads (list): list of surviving saddle indices

    Returns:
        vtk polydata: vtk polydata with connectivity network
    """
    '''
    Collect the geometry of all survivng critical points
        Parameters:
            dir: Geometry type
            dir=0 --> Descending
            dir=1 --> Ascending
            dir=2 --> Both (default)
            dim: Critical point type
            dim=-1 --> All (default)
            dim=0,1,2,3 --> Minima, 1-saddle,2-saddle,Maxima \n"\
    '''
    # Ascending manifold of 2-saddle
    msc.collect_geom(dim=2, dir=1)
    # coordinates of critical points
    dp = msc.dual_points()

    pa, ia = vtk.vtkPoints(), vtk.vtkIntArray()
    fa, ca = vtk.vtkFloatArray(), vtk.vtkCellArray()
    pa.SetData(nps.numpy_to_vtk(dp, "Pts"))
    ia.SetName("SaddleIndex")
    fa.SetName("SaddleVal")
    for s in surviving_sads:
        s = int(s)
        # collect the ascending geom
        gm = msc.asc_geom(s)
        for a, b in gm:
            ca.InsertNextCell(2)
            ca.InsertCellPoint(a)
            ca.InsertCellPoint(b)
            ia.InsertNextValue(s)
            fa.InsertNextValue(msc.cp_func(s))
    pd = vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetLines(ca)
    pd.GetCellData().AddArray(ia)
    pd.GetCellData().AddArray(fa)
    return pd


def get_max(msc):
    """Get the information about the maximas

    Args:
        msc (msc objec): Morse Complex object

    Returns:
        vtk polydata: coordinates, index type, function value, max index
    """
    cps_max = msc.cps(3)
    cpList = []
    for m in cps_max:
        if(msc.cp_func(m) > 0):
            cpList.append((msc.cp_cellid(m), 3, msc.cp_func(m), m))
    cpList = list(set(cpList))

    # create the vtk objects
    pa, ia = vtk.vtkPoints(), vtk.vtkIntArray()
    fa, cp_ids = vtk.vtkFloatArray(), vtk.vtkIntArray()
    ia.SetName("Index")
    fa.SetName("Val")
    cp_ids.SetName("CP ID")

    ca = vtk.vtkCellArray()
    for i, (p, idx, val, cpidx) in enumerate(cpList):
        pa.InsertNextPoint(np.array(p, np.float32)/2)
        ia.InsertNextValue(idx)
        fa.InsertNextValue(val)
        cp_ids.InsertNextValue(cpidx)
        ca.InsertNextCell(1)
        ca.InsertCellPoint(i)

    # Set the outputs
    pd = vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetVerts(ca)
    pd.GetPointData().AddArray(ia)
    pd.GetPointData().AddArray(fa)
    pd.GetPointData().AddArray(cp_ids)
    return pd


def get_saddles(msc, surv_sads):
    """From indices of surviving saddles get further information

    Args:
        msc (msc object): Morse complex object
        surv_sads (list): indices of saddle points that survived

    Returns:
        vtk polydata: coords, index type, function at the saddle points, saddle index, max1, max2
    """
    cpList = []
    for s in surv_sads:
        s, index = int(s), 2  # saddle index - 2
        # maxList is the max critical points connected with s
        co_ords, val, maxList = msc.cp_cellid(
            s), msc.cp_func(s), msc.asc(s)[:, 0]
        cpList.append((co_ords, index, val, s, maxList[0], maxList[1]))

    # cpList = list(set(cpList))

    # create the vtk objects
    # (pa => coords(), ia => index, fa => val)
    pa, ia = vtk.vtkPoints(), vtk.vtkIntArray()
    fa, cp_ids = vtk.vtkFloatArray(), vtk.vtkIntArray()
    max_1, max_2 = vtk.vtkIntArray(), vtk.vtkIntArray()
    ia.SetName("Index")
    fa.SetName("Val")
    cp_ids.SetName("CP ID")
    max_1.SetName("Max 1")
    max_2.SetName("Max 2")
    maxs = []

    ca = vtk.vtkCellArray()
    for i, (p, idx, val, cpidx, m1idx, m2idx) in enumerate(cpList):
        pa.InsertNextPoint(np.array(p, np.float32)/2)
        ia.InsertNextValue(idx)
        fa.InsertNextValue(val)
        max_1.InsertNextValue(m1idx)
        max_2.InsertNextValue(m2idx)
        cp_ids.InsertNextValue(cpidx)
        ca.InsertNextCell(1)
        ca.InsertCellPoint(i)
        maxs.append(m1idx)
        maxs.append(m2idx)

    # Set the outputs
    pd = vtk.vtkPolyData()
    pd.SetPoints(pa)
    pd.SetVerts(ca)
    pd.GetPointData().AddArray(ia)
    pd.GetPointData().AddArray(fa)
    pd.GetPointData().AddArray(cp_ids)
    pd.GetPointData().AddArray(max_1)
    pd.GetPointData().AddArray(max_2)
    return pd, maxs


def get_segmentation_index_dual(msc, img, rtype="VTP"):
    """Get the segmentation from morse smale complex

    Args:
        msc (msc object): Morse Complex object
        img (np array): Distance field
        rtype (str, optional): "VTP" or "NP". Defaults to "VTP".

    Raises:
        TypeError: "VTP" and "NP" are allowed as output

    Returns:
        np array: numpy array as segmentation or vtp as segmenation
    """
    # descending manifold of maxima
    msc.collect_geom(dim=3, dir=0)
    # ''' critical points type
    # dim: Critical point type \n"\
    #     "   dim=-1      --> All (default)\n"\
    #     "   dim=0,1,2,3 --> Minima, 1-saddle,2-saddle,Maxima \n"\
    # dir: Geometry type \n"\
    #     "   dir=0 --> Descending \n"\
    #     "   dir=1 --> Ascending \n"\
    #     "   dir=2 --> Both (default) \n"\
    # '''

    # point coordinates
    dp = msc.dual_points()
    cps_max = msc.cps(3)
    if rtype == "VTP":
        pa, cp_ids = vtk.vtkPoints(), vtk.vtkIntArray()
        ca, val = vtk.vtkCellArray(), vtk.vtkDoubleArray()
        cp_ids.SetName("CP ID")
        val.SetName('Distance Val')
        count = 0
        for m in cps_max:
            if(msc.cp_func(m) <= 0):
                continue
            des_geom = msc.des_geom(m)
            for cube_id in des_geom:
                dual_pt = dp[cube_id]
                if(img[int(dual_pt[0]), int(dual_pt[1]), int(dual_pt[2])] < 0):
                    continue
                val.InsertNextValue(
                    img[int(dual_pt[0]), int(dual_pt[1]), int(dual_pt[2])])
                pa.InsertNextPoint(dual_pt)
                ca.InsertNextCell(1)
                ca.InsertCellPoint(count)
                cp_ids.InsertNextValue(m)
                count += 1
        pd = vtk.vtkPolyData()
        pd.SetPoints(pa)
        pd.SetVerts(ca)
        pd.GetPointData().AddArray(cp_ids)
        pd.GetPointData().AddArray(val)
        return pd
    elif rtype == "NP":
        count = 0
        seg_img = np.full(img.shape, 0, dtype=np.uint32, order='C')
        centers, maxima, labs, vols = [], [], [], []
        for m in cps_max:
            if (msc.cp_func(m) <= 0):
                continue
            des_geom = msc.des_geom(m)

            # we loose information by converting to int <--- issure here
            points = surv_voxs(dp[des_geom].astype(int), img)
            if len(points) == 0:
                pts = collect_neighbours(dp[des_geom])
                points = surv_voxs(pts, img)
                if len(points) == 0:
                    continue
            points_ind = np.ravel_multi_index(points.transpose(), img.shape)
            np.put(seg_img, points_ind, m)
            centers.append(np.mean(points, axis=0))
            maxima.append(np.array(msc.cp_cellid(m), dtype=np.float)/2)
            labs.append(m)
            vols.append(points.shape[0])
            count += 1
        print(f'Number of particles segmented: {count}')
        centers, maxima, labs, vols = np.array(centers), np.array(
            maxima), np.array(labs), np.array(vols)
        return seg_img, centers, maxima, labs, vols
    else:
        raise TypeError("Return type not configured.")
    return None


def read_mat_file(filename, factor=1):
    """Read and downsample mat file into numpy array

    Args:
        filename (str): Name of mat file with volm array
        factor (int, optional): Downsampling factor. Defaults to 1.

    Raises:
        Exception: Not able to read mat file

    Returns:
        numpy array: downsampled array
    """
    try:
        # if saved wih -v7.3 flag
        f = h5py.File(f'{filename}.mat', 'r')
        sci = False
    except OSError:
        # else
        f = io.loadmat(f'{filename}.mat')
        sci = True
    except:
        print("Something went wrong")
        raise Exception('Not able to read the mat file.')

    data = f.get('volm')
    arr = np.array(data)
    if sci:
        arr = np.transpose(arr, axes=[2, 1, 0])

    if (factor > 1):
        arr = block_reduce(arr, (factor, factor, factor),
                           func=np.mean).astype(np.uint16)
    return arr


def read_msc_to_img(msc, dim):
    """This function extracts the Morse Smale complex vertices data to a numpy array.

    Args:
        msc (msc): morse smale complex object
        dim (tuple): dimensions of the image

    Returns:
        np array: numpy array of the function value at msc vertices
    """
    np_arr = np.zeros(dim)
    for x in range(dim[0]):
        for y in range(dim[1]):
            for z in range(dim[2]):
                # scalar value at the vertex coordinate
                np_arr[x, y, z] = msc.vert_func(x, y, z)
    return np_arr  # flatten order - F


def surv_voxs(points, img):
    """This function removes the voxels in background from the list of voxels.

    Args:
        points (np array): array of voxel coordinates
        img (np array): image

    Returns:
        np array: surviving voxel coordinates
    """
    points_ind = np.ravel_multi_index(points.transpose(), img.shape)
    points_val = img.ravel()[points_ind]
    points = points[points_val > 0, ...]
    return points
