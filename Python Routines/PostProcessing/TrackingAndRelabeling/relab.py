# Run Particle Tracking before running this file
#%% 
from re import A
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import SimpleITK as sitk

#%%
filename = "../chamf_distance_CementedSand_12_"

#%% read contact data
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(filename + "0_Cropped_" + "contacts.vtp")
reader.Update()

# %% connectivity in neighbor variable neigh
PolydataOutput = reader.GetOutput()
max1 = vtk_to_numpy(PolydataOutput.GetPointData().GetArray("Max 1"))
max2 = vtk_to_numpy(PolydataOutput.GetPointData().GetArray("Max 2"))
labs = np.unique(np.concatenate((max1, max2), axis=0), axis=0)
neigh = dict()

CN = np.zeros(len(labs))
for ii, l in enumerate(labs):
    # label l particle can be in both max1 and max2
    # so we collect all the neighbors of l from max1 and max2
    first = max2[np.where(max1 == l)]
    second = max1[np.where(max2 == l)]
    neigh[l] = np.unique(np.concatenate((first, second)))
    CN[ii] = len(neigh[l])

#%% depth first search
def dfs(neigh, start, visited, count, maxlab):
    """Recursive depth first search to relabel the particles

        Args:
            neigh (dict): dictionary of neighbors for each particle
            start (int): first particle to start the depth first search
            visited (dict): dictionary of particles that have been visited if not visited then False
            count (int): counter for labeling particles
            maxlab (int): maximum value of labels
    """
    def countadd(x, maxlab):
        x += 1
        if x > maxlab:
            return 0
        else:
            return x

    count = countadd(count, maxlab)
    
    iter = maxlab + 10
    # pick the count which does not belong to neighboring particles
    while count in [visited[l] for l in neigh[start]]:
        count = countadd(count, maxlab)
        iter -= 1
        if iter == 0:
            raise Exception("Exceeds number of iterations")
    
    # count is unique now so we can assign it
    visited[start] = count

    # Here we go over all the neighbors of start and do it recursively
    # this will explore the connected component of the start
    for n in neigh[start]:
        if not visited[n]:
            dfs(neigh, n, visited, count, maxlab)

    # for other connected components
    for key in visited.keys():
        if visited[key] == False:
            dfs(neigh, key, visited, count, maxlab)
    
visited = dict.fromkeys(neigh.keys(), False)
dfs(neigh, np.random.choice(labs, 1)[0], visited, 0, np.max(CN))

#%% read raw file
arr = sitk.GetArrayFromImage(sitk.ReadImage(filename + "0_Cropped_" + "Segmentation.mhd"))


data = np.loadtxt("SegmentedVolumes_" + str(0) + '.csv', dtype='str', delimiter=',', skiprows=1)
data = data.astype('float')
plabs, z = data[:, 0], data[:, 4]

alabs = np.unique(arr)
count = 0
for al in alabs:
    if z[plabs==al] < 1.5:
        arr[arr == al]  = 0
        continue
    if al in labs:
        arr[arr == al] = visited[al]
    elif al == 0:
        continue
    else:
        count += 1
        arr[arr == al] = np.max(CN) + 1

sitk.WriteImage(sitk.GetImageFromArray(arr), filename + "0_Cropped_" + "Relabeled_Segmentation.mhd")

#%% read raw file
visited_fix = visited.copy()
for ii in range(4):
    arr = sitk.GetArrayFromImage(sitk.ReadImage(filename + str(ii+1) + "_Cropped_" + "Segmentation.mhd"))

    data = np.loadtxt("SegmentedVolumes_" + str(ii+1) + '.csv', dtype='str', delimiter=',', skiprows=1)
    data = data.astype('float')
    labs, z = data[:, 0], data[:, 4]

    pairs = np.load(f"pairs_{ii}.npy")
    visited_new = dict()
    for pair in pairs:
        try:
            visited_new[pair[1]] = visited_fix[pair[0]]
        except:
            print(f"In {ii+1} there is no {pair[0]}")
    
    visited_fix = visited_new.copy()

    alabs = np.unique(arr)
    count = 0
    for al in alabs:
        if z[labs==al] < 1.5:
            arr[arr == al]  = 0
            continue
        if al in pairs[:, 1]:
            try:
                arr[arr == al] = visited_fix[al]
            except:
                print(f"In {ii+1} there is no {al}")
                count += 1
                arr[arr == al] = np.max(CN) + 1
        elif al == 0:
            continue
        else:
            count += 1
            arr[arr == al] = np.max(CN) + 1
    
    print("Count: ", count)
    sitk.WriteImage(sitk.GetImageFromArray(arr), filename + str(ii+1) + "_Cropped_" + "Relabeled_Segmentation.mhd")

# %%
# count = 0
# for al in alabs:
#     if al not in labs:
#         print(al)
#         count += 1