# import modules
import SimpleITK as sitk
import vtk

# Here we can also add functions to write polydata from numpy array instead of writing it in main.py


def write_img_from_arr(arr, file_name):
    """Write a numpy array to sitk image and save it to file.

    Args:
        arr (numpy array): array to be written
        file_name (str): file name to be written
    """
    # conventinally numpy -- z, y, x and sitk -- x, y, z
    img = sitk.GetImageFromArray(arr.transpose(2, 1, 0))
    writer = sitk.ImageFileWriter()
    writer.SetFileName(file_name + '.mhd')
    writer.Execute(img)


def write_polydata(pd, file_name):
    """Write vtk polydata to file.

    Args:
        pd (vtk polydata): polydata to be written
        file_name (str): file name to be written
    """
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(file_name)
    writer.SetInputData(pd)
    writer.Write()
    return
