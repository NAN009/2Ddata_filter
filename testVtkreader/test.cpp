#include "vtkReader.h"
#include "vtk.h"
#include <iostream>
 VTK vtk;
 vtkReader vr;
 void test2d()
 {
	 std::string path = "D:/newData/taiwan_2/taiwan_2_double.nak";
	 vtk.Load(path, 0);

	 vtk.GaussianFilter3D(3,0.5,0.01);
 }
 void testVr()
 {
	 std::string path = "D:/2Ddata/2d.vtk";
	 vtk.Load(path, 0);

	 vtk.dataFilter();
 }
 void testFourier()
 {
	 std::string path = "D:/2Ddata/2d.vtk";
	 vtk.Load(path, 0);

	 vtk.FourierValue2D(3.1415926, 3.1415926);
 }