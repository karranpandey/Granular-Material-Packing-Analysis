import vtk
import SimpleITK as sitk

spheres=[]

def getSphere(x,y,z,i,d):
	sphere = vtk.vtkSphereSource()
	centre=[0,0,0]
	centre[i]+=d*16	
	sphere.SetCenter(centre[0],centre[1],centre[2])
	sphere.SetRadius(10)
	sphere.Update()
	return sphere.GetOutput()

def genSphere(x,y,z,r):
	sphere = vtk.vtkSphereSource()
	centre=[x,y,z]
	sphere.SetCenter(centre[0],centre[1],centre[2])
	sphere.SetRadius(r)
	sphere.Update()
	return sphere.GetOutput()

def genPacking(num_x,num_y,num_z,r,close_frac):
	sphere_packing=[]
	for i in range(num_x):
		for j in range(num_y):
			for k in range(num_z):
				curr_centre=[2*r*i,2*r*j,2*r*k]
				sphere_packing.append(genSphere(curr_centre[0],curr_centre[1],curr_centre[2],r*close_frac))
	return sphere_packing	

def init_Spheres():
	spheres.append(getSphere(0,0,0,0,0))
	spheres.append(getSphere(0,0,0,0,1))
	spheres.append(getSphere(0,0,0,1,1))
	spheres.append(getSphere(0,0,0,2,1))
	spheres.append(getSphere(0,0,0,0,-1))
	spheres.append(getSphere(0,0,0,1,-1))
	spheres.append(getSphere(0,0,0,2,-1))
	return



spheres=genPacking(1,1,2,15,1.01)
#init_Spheres()
pd=vtk.vtkPolyData()



appendFilter=vtk.vtkAppendPolyData()

#print(spheres[0])

pd=spheres[0]

print(len(spheres))

for i in range(len(spheres)):
	appendFilter.AddInputData(spheres[i])
	#boolOp=vtk.vtkBooleanOperationPolyDataFilter()
	#boolOp.SetOperationToUnion()
	#boolOp.SetInputData(0,pd)
	#boolOp.SetInputData(1,spheres[i+1])
	#boolOp.Update()
	#pd=boolOp.GetOutput()
	#print(i)

appendFilter.Update()
pd=appendFilter.GetOutput()



whiteImage = vtk.vtkImageData()    
bounds=pd.GetBounds()
spacing=[1,1,1]
whiteImage.SetSpacing(spacing)
dim=[0,0,0]

for i in range(3):
    dim[i] = int((bounds[i * 2 + 1] - bounds[i * 2])/spacing[i])
 
whiteImage.SetDimensions(dim);
whiteImage.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1)

origin=[0,0,0]

origin[0] = bounds[0] + spacing[0] / 2
origin[1] = bounds[2] + spacing[1] / 2
origin[2] = bounds[4] + spacing[2] / 2
  
whiteImage.SetOrigin(origin);

whiteImage.AllocateScalars(vtk.VTK_FLOAT,1)

#print(whiteImage)

inval = 255.0
outval = 0.0
count = whiteImage.GetNumberOfPoints();

for i in range(count):
    whiteImage.GetPointData().GetScalars().SetTuple1(i, inval)
 

pol2stenc = vtk.vtkPolyDataToImageStencil()

pol2stenc.SetInputData(pd)
pol2stenc.SetOutputOrigin(origin)
pol2stenc.SetOutputSpacing(spacing)
pol2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
pol2stenc.Update()

imgstenc = vtk.vtkImageStencil()

imgstenc.SetInputData(whiteImage)
imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
imgstenc.ReverseStencilOff()
imgstenc.SetBackgroundValue(outval)
imgstenc.Update()

image_data=(imgstenc.GetOutput())

filename = "spheres.vti"

writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName(filename)
writer.SetInputData(image_data)
writer.Write()

meta_writer=vtk.vtkMetaImageWriter()
meta_writer.SetFileName("spheres.mhd")
meta_writer.SetRAWFileName("spheres.raw")
meta_writer.SetCompression(False)
meta_writer.SetInputData(image_data)
meta_writer.Write()


#print("hello")

reader=sitk.ImageFileReader()
reader.SetImageIO("MetaImageIO")
reader.SetFileName("spheres.mhd")

image=reader.Execute()

thresholdFilter=sitk.BinaryThresholdImageFilter()

thresh_image=thresholdFilter.Execute(image,0,100,1,0)

signedDistFilter=sitk.SignedMaurerDistanceMapImageFilter()

signed_dist_image=signedDistFilter.Execute(thresh_image)

raw_writer = sitk.ImageFileWriter()
raw_writer.SetFileName("spheres_dist.mhd")
raw_writer.Execute(signed_dist_image)


