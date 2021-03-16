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

Additionally, make sure to set the PYTHON_SITE_PACKAGES_DIR in cmake to the 'Python Routines' folder in the repository. 

### Running the Pipeline

The python scripts to run the pipeline can be found in the Python Routines folder. The scripts and their input formats are described below: 

#### bd_extraction.py 

This script takes as input the raw CT image and outputs the distance field based on the extracted boundary. It takes as input a MetaImage file (.mhd + .raw). To run the script, execute the following command in the terminal: 

`python bd_extraction.py [Path to .mhd file]` 

This will store the computed distance field in MetaImage format (.mhd + .raw) in the 'Outputs' folder in the repository.

---
