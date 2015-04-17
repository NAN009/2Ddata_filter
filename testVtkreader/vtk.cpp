#include <math.h>
#include "vtk.h"
#include <iostream>
#include <fstream>

using namespace std;
#define DATA(x,y) (data[(x)+dimX*(y)])


void VTK::Load(std::string filename, double threshold)
{
	filebuf fb;
	double val;
	if (fb.open(filename, ios::in))
	{
		istream is(&fb);
		string info;
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		is >> info >> dimX >> dimY;
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		getline(is, info);
		unsigned int nelems = dimX*dimY;
		data = new double[nelems];
		minVal = 1e10; maxVal = -1e10;
		for (unsigned int i = 0; i < nelems; ++i)
		{
			is >> val;
			if (val < threshold)
				val = threshold;
			data[i] = (val - threshold);
			if (data[i] < 0)data[i] = 0;
			if (data[i] < minVal)minVal = data[i];
			if (data[i]>maxVal)maxVal = data[i];
		}
		cout << "Read VTK file completed" << endl;
	}


}
void VTK::GaussianFilter(int radius, double sigma, double ratio)
{
	int i, j, x, y;
	const int diameter = 7;
	cout << data[0] << " " << data[1] << " " << data[10] << endl;
	cout << data[2150] << " " << data[2151] << " " << data[2152] << " " << data[2153] << " " << data[2154] << endl;
	cout << data[2364] << " " << data[2365] << " " << data[2366] << " " << data[2367] << " " << data[2368] << endl;
	int *ix, *iy;
	double *fileter, temVal, *temData;

	fileter = new double[diameter*diameter];
	ix = new int[dimX + 2 * radius];
	iy = new int[dimY + 2 * radius];
	for (i = 0; i < radius; ++i)
	{
		ix[i] = iy[i] = 0;
		ix[dimX + radius + i] = dimX - 1;
		iy[dimY + radius + i] = dimY - 1;
	}
	for (i = 0; i < dimX; ++i)
		ix[i + radius] = i;
	for (i = 0; i < dimY; ++i)
		iy[i + radius] = i;
	for (i = 0; i <= radius; ++i)
	{
		for (j = 0; j <= radius; ++j)
		{
			temVal = sqrt(i*i + j*j);
			temVal = ratio*sqrt(exp(-0.5*(temVal / sigma)*(temVal / sigma)));
			fileter[(radius - j) + diameter*(radius - i)] = temVal;
			fileter[(radius - j) + diameter*(radius + i)] = temVal;
			fileter[(radius + j) + diameter*(radius - i)] = temVal;
			fileter[(radius + j) + diameter*(radius + i)] = temVal;
		}
	}
	temVal = -fileter[radius + diameter*radius];
	for (i = 0; i < diameter; ++i)
	{
		for (j = 0; j < diameter; ++j)
			temVal += fileter[j + diameter*i];
	}
	ofstream value("D:\\2Ddata\\2d_filter.txt");
	fileter[radius + diameter*radius] = 1.0 - temVal;
	temData = new double[dimX*dimY];
	for (y = 0; y < dimY; ++y)
	{
		for (x = 0; x < dimX; ++x)
		{
			temVal = 0;
			for (i = 0; i < diameter; ++i)
			{
				for (j = 0; j < diameter; ++j)
					temVal += fileter[j + diameter*i] * DATA(ix[x+j],iy[y+i]);
			}
			temData[x + dimX*y] = temVal;
			value << temVal << " ";
		}
		value << endl;
	}
	delete[] data;
	data = temData;
	cout << data[0] <<" "<<data[1]<<" "<<data[10]<< endl;
	cout << data[2150] << " " << data[2151] << " " << data[2152] << " " << data[2153] << " " << data[2154] << endl;
	cout << data[2364] << " " << data[2365] << " " << data[2366] << " " << data[2367] << " " << data[2368] << endl;
	delete[] iy;
	delete[] ix;
	delete[] fileter;
	cout << "data points are smoothened by Gaussian filter" << endl;
}

void VTK::dataFilter()
{
	ofstream value("D:\\1data\\iran_dataFilter.txt");
	double *temData;
	temData = new double[dimX*dimY];
	double lacalPoint1, loaclPoint2;

	

	for (int i = 0; i < dimY; ++i)
	{
		for (int j = 0; j < dimX; ++j)
		{
			if (i<1 || i>=dimX - 1 || j<1 || j>=dimY - 1)
				temData[j + i*dimX] = DATA(j,i);
			else
			{
				lacalPoint1 = 3 * DATA(j, i)/2 ;

				loaclPoint2 = DATA(j + 1, i);
				loaclPoint2 += DATA(j - 1, i);

				
				loaclPoint2 += DATA(j, i + 1);
				loaclPoint2 += DATA(j, i - 1);
				
			


				temData[j + i*dimX] = lacalPoint1 - 1 * loaclPoint2 / 8;
			}
			
			value << temData[j + i*dimX] << " ";
		}
		value << endl;
	}
	delete[] data;
	data = temData;
}

#undef DATA
#undef INDEX    
