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

	//����ԭʼ�ļ���صķ��� ��ʹ���κο����㷨
	void rawPorcess(const Mat& src, Mat& moments, int p, int q);
	//����ԭʼ���ؽ�����
	void rawReconstruct(const Mat& moments, Mat& dst, int p, int q, int N);

	//"fast computation of legendre and zernike moments"--R.MUKUNDAN and R.RAMAKRISHNAN
	void mukundanProcess(const Mat& src,const Mat& points, Mat& moments, int order,int N);
	void mukundanReconstruct(const Mat& moments,Mat& dst,int order,int N);

	void getPoly(Mat& poly,int,int);
private:
	double pnx(int n, float x);
};


