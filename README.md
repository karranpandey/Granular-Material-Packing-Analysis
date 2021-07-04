# Granular-Material-Packing-Analysis

## Description 

Repository containing code for a discrete Morse theory based method to analyse 3-D CT scans of granular material packings. 

(e-mail: karran13 [AT] gmail.com)

---

## Requirements

To compile and execute this code, you will need Cmake > 3.1.8, Boost > 1.58, OpenCL 1.1 (implicitly available along with CUDA > 10.1), OpenMP and python with the following packages: 

1. vtk
2. scikit-image
3. SimpleITK
4. ITK
5. numpy/scipy

---

## Instructions

### Installing pyms3d

The MS-Complex computation is done using pyms3d. The pyms3d source code is in the ms complex software folder in the repository. To install pyms3d, follow the basic process outlined in the <a href = "https://bitbucket.org/vgl_iisc/mscomplex-3d/wiki/Installation"> the pyms3d installation wiki </a>. 

Additionally, make sure to set the PYTHON_SITE_PACKAGES_DIR in cmake to the 'Python Routines' folder in the repository. This should ideally put the 'pyms3d.so' file in the 'Python Routines' folder as required.  

### Running the Pipeline

The python scripts to run the pipeline can be found in the Python Routines folder. The scripts and their input formats are described below: 

#### bd_extraction.py 

This script takes as input the raw CT image and outputs the distance field based on the extracted boundary. It takes as input a MetaImage file (.mhd + .raw). To run the script, execute the following command in the terminal: 

`python bd_extraction.py [Path to .mhd file]` 

This will store the computed distance field in MetaImage format (.mhd + .raw) in the 'Outputs' folder in the repository.

#### downsample_skimage.py 

This script is a useful downsampling routine if the original CT image doesn't fit in memory for the MS-Complex computation. It takes as input the raw CT image and downsamples it based on the given factor. It takes as input a MetaImage file (.mhd + .raw). To run the script, execute the following command in the terminal: 

`python downsample_skimage.py [Path to .mhd file] [factor]` 

This will store the downsampled CT image in MetaImage format (.mhd + .raw) in the 'Outputs' folder in the repository.

#### main.py 

This script is the main interface to run the Morse-Smale Complex computation and extract relevant geometric and topological structures for analysis. It takes the computed distance field as input and returns the structures selected from the in-program menu. The program allows for the visualization of the persistence curve, computation/simplification of the MS-Complex and extraction of the segmentation, connectivity network and contact regions in the granular material packing.  

`python main.py [Path to .raw file of distance field] [dimension 1] [dimension 2] [dimension 3]` 

Running this will store the selected structures in '.vtp' format (accessible through VTK/ParaView) in the 'Outputs' folder in the repository.

Few notes here: Use the 'knee' in the persistence curve to select a simplification threshold. Secondly, Make sure the dimension order is as per the mhd file description of the distance field.

---
