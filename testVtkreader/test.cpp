#include "vtkReader.h"
#include "vtk.h"
#include <iostream>
 VTK vtk;
 vtkReader vr;
 void test2d()
 {
	 std::string path = "D:/1data/iran.vtk";
	 vtk.Load(path, 0);

	 vtk.dataFilter();
 }
 void testVr()
 {
	 vr.loadFile("D:\\2Ddata\\2d.vtk");
	

	 vtk.GaussianFilter(3, 0.5, 0.05);
 }