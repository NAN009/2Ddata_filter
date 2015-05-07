#include <math.h>
#include "vtk.h"
#include <iostream>
#include <fstream>
#include <iomanip>
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
	const int diameter =11;
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
	ofstream value("D:\\newData\\taiwan_2\\taiwan_2_gray_filter.nak");
	value << "#vtk DataFile Version 3.0" << endl;
	value << "Created by ZhangNan!" << endl;
	value << "ASCII" << endl;
	value << "DATASET STRUCTURED_POINTS" << endl;
	value << "DIMENSIONS "<<dimX<<" "<<dimY << endl;
	value << "ORIGIN -3 -3 0" << endl;
	value << "SPACING 1 1 1" << endl;
	value << "POINT_DATA " <<dimX*dimY<< endl;
	value << "SCALARS image_data double" << endl;
	value << "LOOKUP_TABLE default" << endl;
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
			value << setiosflags(ios::fixed) << setprecision(16) << temVal << " ";

		}
		value << endl;
	}
	delete[] data;
	data = temData;
	cout << setiosflags(ios::fixed) << setprecision(16) << data[0] << " " << data[1] << " " << data[10] << endl;
	cout << data[2150] << " " << data[2151] << " " << data[2152] << " " << data[2153] << " " << data[2154] << endl;
	cout << data[2364] << " " << data[2365] << " " << data[2366] << " " << data[2367] << " " << data[2368] << endl;
	delete[] iy;
	delete[] ix;
	delete[] fileter;
	cout << "data points are smoothened by Gaussian filter" << endl;
}

void VTK::dataFilter()
{
	ofstream value("D:\\2Ddata\\2d_dataFilter.vtk");
	value << "#vtk DataFile Version 3.0" << endl;
	value << "Created by ZhangNan!" << endl;
	value << "ASCII" << endl;
	value << "DATASET STRUCTURED_POINTS" << endl;
	value << "DIMENSIONS " << 4*dimX-3 << " " << 4*dimY-3<< endl;
	value << "ORIGIN -3 -3 0" << endl;
	value << "SPACING 1 1 1" << endl;
	value << "POINT_DATA " << (4 * dimX - 3) * 4 * dimY - 3 << endl;
	value << "SCALARS image_data double" << endl;
	value << "LOOKUP_TABLE default" << endl;
	vector <vector< double > > temData(dimY*4-3);
	for (int i = 0; i < dimY * 4 -3; ++i)
		temData[i].resize(dimX * 4 -3);
	//double **temData = new double[(dimY - 1) * 4 + 1][(dimX - 1) * 4 + 1];
	double lacalPoint1, loaclPoint2;
	int i, j, j1, i1;
	for (i = 0; i < dimY-1; ++i)
	{
		for (j = 0; j < dimX-1; ++j)
		{
			
			temData[i*4][j*4] = DATA(j, i);
			temData[i*4][j*4+1] = DATA(j, i) * 3 / 4 + DATA(j+1, i) * 1 / 4;
			temData[i*4][j*4+2] = DATA(j, i) * 1 / 2 + DATA(j+1, i) * 1 / 2;
			temData[i*4][j*4+3] = DATA(j, i) * 1 / 4 + DATA(j+1, i) * 3 / 4;
			
			temData[i*4+1][j*4] = DATA(j, i) * 3 / 4 + DATA(j , i+1) * 1 / 4;
			temData[i*4 + 1][j*4 + 1] = temData[i*4 + 1][j*4] * 3 / 4 + (DATA(j + 1, i) * 3 / 4 + DATA(j + 1, i+1) * 1 / 4) * 1 / 4;
			temData[i*4 + 1][j*4 + 2] = temData[i*4 + 1][j*4] * 1 / 2 + (DATA(j + 1, i) * 3 / 4 + DATA(j + 1, i + 1) * 1 / 4) * 1 / 2;
			temData[i*4 + 1][j*4 + 3] = temData[i*4 + 1][j*4] * 1 / 4 + (DATA(j + 1, i) * 3 / 4 + DATA(j + 1, i + 1) * 1 / 4) * 3 / 4;
			
			temData[i*4 + 2][j*4] = DATA(j, i) * 1 / 2 + DATA(j, i + 1) * 1 / 2;
			temData[i*4 + 2][j*4 + 1] = temData[i*4 + 2][j*4] * 3 / 4 + (DATA(j + 1, i) * 1 / 2 + DATA(j + 1, i + 1) * 1 / 2) * 1 / 4;
			temData[i*4 + 2][j*4 + 2] = temData[i*4 + 2][j*4] * 1 / 2 + (DATA(j + 1, i) * 1 / 2 + DATA(j + 1, i + 1) * 1 / 2) * 1 / 2;
			temData[i*4 + 2][j*4 + 3] = temData[i*4 + 2][j*4] * 1 / 4 + (DATA(j + 1, i) * 1 / 2 + DATA(j + 1, i + 1) * 1 / 2) * 3 / 4;

			temData[i*4 + 3][j*4] = DATA(j, i) * 1 / 4 + DATA(j, i + 1) * 3 / 4;
			temData[i*4 + 3][j*4 + 1] = temData[i*4 + 3][j*4] * 3 / 4 + (DATA(j + 1, i) * 1 / 4 + DATA(j + 1, i + 1) * 3 / 4) * 1 / 4;
			temData[i*4 + 3][j*4 + 2] = temData[i*4 + 3][j*4] * 1 / 2 + (DATA(j + 1, i) * 1 / 4 + DATA(j + 1, i + 1) * 3 / 4) * 1 / 2;
			temData[i*4 + 3][j*4 + 3] = temData[i*4 + 3][j*4] * 1 / 4 + (DATA(j + 1, i) * 1 / 4 + DATA(j + 1, i + 1) * 3 / 4) * 3 / 4;


		}
		if (j == dimX - 1)
		{
			temData[i * 4][j * 4] = DATA(j, i);
			temData[i * 4 + 1][j * 4] = DATA(j, i) * 3 / 4 + DATA(j, i + 1) * 1 / 4;
			temData[i * 4 + 2][j * 4] = DATA(j, i) * 1 / 2 + DATA(j, i + 1) * 1 / 2;
			temData[i * 4 + 3][j * 4] = DATA(j, i) * 1 / 4 + DATA(j, i + 1) * 3 / 4;
		}	
	}
	for (int i = 0; i < dimX-1; ++i)
	{
		temData[4 * (dimY - 1)][i * 4] = DATA(i, dimY - 1);
		temData[4 * (dimY - 1)][i * 4 + 1] = DATA(i, dimY - 1) * 3 / 4 + DATA(i+1, dimY - 1) * 1 / 4;
		temData[4 * (dimY - 1)][i * 4 + 2] = DATA(i, dimY - 1) * 1 / 2+ DATA(i + 1, dimY - 1) * 1 / 2;
		temData[4 * (dimY - 1)][i * 4 + 3] = DATA(i, dimY - 1) * 1 / 4 + DATA(i + 1, dimY - 1) * 3 / 4;
	}
	temData[4 * (dimY - 1)][4 * (dimX - 1)] = DATA(dimX-1,dimY-1);
	for (int i = 0; i < dimY*4-3 ; ++i)
	{
		for (int j = 0; j < dimX*4-3 ; ++j)
		{
			value << setiosflags(ios::fixed) << setprecision(9) << temData[i][j] << " ";
		}
		value << endl;
	}
	
}

void VTK::FourierValue2D(double w1,double w2)
{
	ofstream value("D:\\2Ddata\\2d_FourierFilter.vtk");
	value << "#vtk DataFile Version 3.0" << endl;
	value << "Created by ZhangNan!" << endl;
	value << "ASCII" << endl;
	value << "DATASET STRUCTURED_POINTS" << endl;
	value << "DIMENSIONS " << dimX << " " << dimY << endl;
	value << "ORIGIN -3 -3 0" << endl;
	value << "SPACING 1 1 1" << endl;
	value << "POINT_DATA " << dimX*dimY << endl;
	value << "SCALARS image_data double" << endl;
	value << "LOOKUP_TABLE default" << endl;

	double value1 = 0,value2=0;
	double *cosW1 = new double[dimX];
	double *sinW1 = new double[dimX];
	double *cosW2 = new double[dimY];
	double *sinW2 = new double[dimY];
	for (int i = 0; i < dimX; ++i)
	{
		sinW1[i] = sin(i*w1);
		cosW1[i] = cos(i*w1);
	}
	for (int i = 0; i < dimY; ++i)
	{
		sinW2[i] = sin(i*w2);
		cosW2[i] = cos(i*w2);
	}
	for (int i = 0; i < dimX; ++i)
	{
		for (int j = 0; j < dimY; ++j)
		{
			double value1 = 0, value2 = 0;
			value1 += DATA(i, j)*(cosW1[i] * cosW2[j] - sinW1[i] * sinW2[j]);
			value2 += DATA(i, j)*(sinW1[i] * cosW2[j] - cosW1[i] * sinW2[j]);
			value << setiosflags(ios::fixed) << setprecision(16) << abs(value1 + value2) << " ";
		}
		value << endl;
	}
}
#undef DATA

