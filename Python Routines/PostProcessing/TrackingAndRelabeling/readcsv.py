import numpy as np
from matplotlib import pyplot as plt

def vol(filename):
    data = np.loadtxt(filename, dtype='str', delimiter=',', skiprows=1)
    data = data.astype('float')
    label, volume, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    volfilter = volume[z>0.5]
    return volfilter

vol0 = vol('SegmentedVolumes_0.csv')
vol1 = vol('SegmentedVolumes_1.csv')
vol2 = vol('SegmentedVolumes_2.csv')
vol3 = vol('SegmentedVolumes_3.csv')

vol = np.array([vol0.sum(), vol1.sum(), vol2.sum(), vol3.sum()])
time = np.arange(0, 4, 1)

plt.figure()
plt.plot(time, vol, 'ko')
plt.xlabel('Stage')
plt.ylabel('Volume')
plt.show()