#pragma once
#include <string>
#include <list>
#include <vector>

typedef struct DataIndex
{
	int x, y;
	DataIndex(int xx = -1, int yy = -1)
	{
		x = xx;
		y == yy;
	}
};
#define  FLAGBG 1
#define  FLAGBORDER -1

class VTK
{
public:
	int dimX, dimY;
	double *data;
	char *bgFlag;
	std::string file_path;
	double minVal, maxVal;
	VTK()
	{
		dimX = 0; dimY = 0;
		file_path = "";
		data = nullptr;
		bgFlag = nullptr;
	}
	~VTK()
	{
		if (data)
			delete[] data;
		if (bgFlag)
			delete[] bgFlag;
	}
	void Load(std::string filename, double threshold = 0);
	void GaussianFilter(int radius, double sigma, double ratio);

};