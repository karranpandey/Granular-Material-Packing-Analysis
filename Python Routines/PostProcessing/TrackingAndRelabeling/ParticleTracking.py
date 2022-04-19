# %% import modules
import os
import glob
import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment as match
from scipy.spatial.distance import cdist

# %% directory for files
directory = '.'
csv_files = sorted(glob.glob(os.path.join(directory, 'SegmentedVolumes_*.csv')))

# %% read the files for centers, vol
features = pd.DataFrame()
for num, csv_file in enumerate(csv_files):
    data = pd.read_csv(csv_file)
    lab, vol, x, y, z = data["Label"], data["Volume"], data["Points:0"], data["Points:1"], data["Points:2"]

    for jj in range(len(vol)):
        if z[jj] < 1.5:
            continue
        features = features.append([{'x': x[jj],
                                     'y': y[jj],
                                     'z': z[jj],
                                     'lab': lab[jj],
                                     'vol': vol[jj],
                                     'frame': num,
                                     }, ])

# %% Cost matrix
def index_match(features, frame_1, frame_2):
    """[summary]

    Args:
        features (pandas dataframe): with center, lab, vol, and frame number
        frame_1 (int): first frame number
        frame_2 (int): second frame number

    Returns:
        [(int array, int array)]: pairs of matched index from frames
    """
    vol_1, vol_2 = features[features['frame'] == frame_1]['vol'], features[features['frame'] == frame_2]['vol']
    lab_1, lab_2 = features[features['frame'] == frame_1]['lab'], features[features['frame'] == frame_2]['lab']
    x_1, x_2 = features[features['frame'] == frame_1]['x'], features[features['frame'] == frame_2]['x']
    y_1, y_2 = features[features['frame'] == frame_1]['y'], features[features['frame'] == frame_2]['y']
    z_1, z_2 = features[features['frame'] == frame_1]['z'], features[features['frame'] == frame_2]['z']

    center_1 = np.concatenate((x_1.to_numpy()[:, np.newaxis], y_1.to_numpy()[:, np.newaxis], z_1.to_numpy()[:, np.newaxis]), axis=1)
    center_2 = np.concatenate((x_2.to_numpy()[:, np.newaxis], y_2.to_numpy()[:, np.newaxis], z_2.to_numpy()[:, np.newaxis]), axis=1)
    rad_1, rad_2 = ((3/(4*np.pi))*vol_1.to_numpy())**(1/3), ((3/(4*np.pi))*vol_2.to_numpy())**(1/3)
    rad_1, rad_2 = rad_1[:, np.newaxis], rad_2[:, np.newaxis]
    lab_1, lab_2 = lab_1.to_numpy(), lab_2.to_numpy()

    rad_diff = np.abs(rad_1-np.transpose(rad_2))
    rad_diff = rad_diff / np.amax(rad_diff, axis=0)
    ind_rad = rad_diff > 0.3

    pair_dist = cdist(center_1, center_2)
    ind_pair = pair_dist > 40
    pair_dist = pair_dist / np.amax(pair_dist, axis=0)

    weight_dist, weight_rad = 0.9, 0.3
    cost = weight_dist * pair_dist + weight_rad * rad_diff

    cost[ind_rad] = 999999
    cost[ind_pair] = 999999

    # bipartite graph matching
    match_1, match_2 = match(cost)
    lab_1, lab_2 = lab_1[match_1], lab_2[match_2]

    return lab_1, lab_2


for ii in range(len(csv_files)-1):
    l1, l2 = index_match(features, ii, ii+1)
    l1, l2 = l1[:, np.newaxis], l2[:, np.newaxis]
    pairs = np.concatenate((l1, l2), axis=1)
    np.save('pairs_'+str(ii), pairs)
    # np.savetxt(f'pairs_{ii+1:04d}.txt', pairs, fmt='%d')
    # print(pairs)