#pragma once
#include <opencv2\core\core.hpp>
#include <iostream>
#include "Tool.h"

using namespace std;
using namespace cv;

class Legendre
{
public:
	Legendre();
	~Legendre();
	void process(const Mat& src, Mat& moments, int p, int q);
	void reconstruct(const Mat& moments, Mat& dst, int p, int q, int N);

	//最最原始的计算矩的方法 不使用任何快速算法
	void rawPorcess(const Mat& src, Mat& moments, int p, int q);
	//最最原始的重建方法
	void rawReconstruct(const Mat& moments, Mat& dst, int p, int q, int N);

	//"fast computation of legendre and zernike moments"--R.MUKUNDAN and R.RAMAKRISHNAN
	void mukundanProcess(const Mat& src,const Mat& points, Mat& moments, int order,int N);
	void mukundanReconstruct(const Mat& moments,Mat& dst,int order,int N);

	void getPoly(Mat& poly,int,int);
private:
	double pnx(int n, float x);
};


