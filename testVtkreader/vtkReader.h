#ifndef _VTK_READER_H
#define _VTK_READER_H
#include "minmax.h"
#include <fstream>
#include <iostream>
#include<math.h>
#pragma warning (disable : 4996)
using namespace std;

class vtkReader {
public:
	vtkReader() : dataArr(NULL) {}
	~vtkReader() {
		if ( dataArr != NULL ) {
			delete []dataArr;
		}
		if ( upperArr != NULL ) {
			delete []upperArr;
		}
		if ( lowerArr != NULL ) {
			delete []lowerArr;
		}
	}
	const int getSize() const {
		return size;
	}
	const int* getDimension() const {
		return dim;
	}
	const double* getDataArray() const {
		return dataArr;
	}
	const double getData(int x, int y) const {
		return dataArr[ getIndex(x,y) ];
	}
	const double getUpper(int x, int y) const {
		return upperArr[ getIndex(x,y) ];
	}
	const double getLower(int x, int y) const {
		return lowerArr[ getIndex(x,y) ];
	}
	int getIndex(int x, int y) const {
		//先读列
		return max(min(x,dim[0]-1),0) + max(min(y,dim[1]-1),0) * dim[0];
	}//返回（i,j）为数组中第几个数
	bool loadFile(const char* filename);
	int dim[2];
private:
	
	double *dataArr, *upperArr, *lowerArr;
	int size;
	
};

bool vtkReader::loadFile(const char* filename) {
	char cachedFilename[1024];
	strcpy(cachedFilename, filename);
	strcat(cachedFilename, ".cache");

	FILE *fp = fopen(cachedFilename, "rb");
	if ( fp!=NULL ) 
	{
		printf("Reading from cache file.\n");
		fread(dim, sizeof(dim), 2, fp);
		printf("Dimension: %d * %d\n", dim[0], dim[1]);
		size = dim[0]*dim[1];
		dataArr = new double[size];
		upperArr = new double[ size ];
		lowerArr = new double[ size ];
		fread(dataArr, sizeof(double), size, fp);
		fread(upperArr, sizeof(double), size, fp);
		fread(lowerArr, sizeof(double), size, fp);
		fclose(fp);
		return true;
	}
	else 
	{
		fp = fopen(filename, "r");
		if ( fp==NULL ) {
			fprintf(stderr, "File %s not found.\n", filename);
			return false;
		}
	}
	char ignored[1024];
	for(int i=0; i<4; i++) {
		fgets(ignored, sizeof(ignored), fp);
		printf("Ignored line: %s", ignored);
	}
	fscanf(fp, "%s%d%d", ignored, &dim[0], &dim[1]);
	//fin >> ignored >> dim[0] >> dim[1] >> dim[2];
	printf("Dimension information: %d * %d\n", dim[0], dim[1]);
	for(int i=0; i<6; i++) {
		fgets(ignored, sizeof(ignored), fp);
		printf("Ignored line: %s", ignored);
	}
	size = dim[0]*dim[1];
	dataArr = new double[ size ];
	upperArr = new double[ size ];
	lowerArr = new double[ size ];
	for(int i=0; i<size; i++) {
		fscanf(fp, "%lf", dataArr+i);
	}
	fclose(fp);
	// Calculating upper and lower bound
	for(int i=0; i<dim[0]; i++) 
	{
		for(int j=0; j<dim[1]; j++) 
		{
			//for(int k=0; k<dim[2]; k++) {
				int index = getIndex(i,j);
				upperArr[index] = -1e300;
				lowerArr[index] = 1e300;
				for(int a=-2; a<=2; a++)
					for(int b=-2; b<=2; b++)
					{
						double data = getData(i + a, j + b);
						if ( upperArr[index] < data )
						{
							upperArr[index] = data;
						}
						if ( lowerArr[index] > data )
						{
							lowerArr[index] = data;
						}
						
				   }
		}
	}
	cout << dataArr[499] << " " << dataArr[500] << " " << dataArr[501] << " " << dataArr[502] << " " << dataArr[503] << endl;

	//set filter
	int radius = 3, ratio = 0.5, sigma = 0.01;
	int i, j, x, y;
	const int diameter = 7;
	int *ix, *iy;
	double *fileter, temVal, *temData;
	fileter = new double[diameter*diameter];
	ix = new int[dim[0] + 2 * radius];
	iy = new int[dim[1] + 2 * radius];
	for (i = 0; i < radius; ++i)
	{
		ix[i] = iy[i] = 0;
		ix[dim[0] + radius + i] = dim[0] - 1;
		iy[dim[1] + radius + i] = dim[1] - 1;
	}
	for (i = 0; i < dim[0]; ++i)
		ix[i + radius] = i;
	for (i = 0; i < dim[1]; ++i)
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
	fileter[radius + diameter*radius] = 1.0 - temVal;
	temData = new double[dim[0]*dim[1]];

	for (y = 0; y < dim[1]; ++y)
	{
		for (x = 0; x < dim[0]; ++x)
		{
			temVal = 0;
			for (i = 0; i < diameter; ++i)
			{
				for (j = 0; j < diameter; ++j)
					temVal += fileter[j + diameter*i] * getData(ix[x + j], iy[y + i]); 
			}
			temData[x + dim[0]*y] = temVal;
		}
	}
	
	dataArr = temData;
	cout << dataArr[499] << " " << dataArr[500] << " " << dataArr[501] << " " << dataArr[502] << " " << dataArr[503] << endl;
	delete[] iy;
	delete[] ix;
	delete[] fileter;

	// Writing to cache file
	fp = fopen(cachedFilename, "wb");
	fwrite(dim, sizeof(dim), 2, fp);
	fwrite(dataArr, sizeof(double), size, fp);
	fwrite(upperArr, sizeof(double), size, fp);
	fwrite(lowerArr, sizeof(double), size, fp);
	fclose(fp);
	return true;
}

#endif