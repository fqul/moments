#pragma once
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>

#include <iostream>
#include <vector>

using namespace std;
using namespace cv;

enum RECURSION_TYPE{
	BY_X, BY_N
};

class Tchebichef
{
public:
	Tchebichef();
	~Tchebichef();
	//求矩：
	//src: 求矩形的图片；t：存放矩；N图片的大小，因此矩的阶为N-1+N-1；rtype：按照x递推，或者按照n递推
	void moments(const Mat& src,Mat& t,Mat& tnp,int N,RECURSION_TYPE rtype);

	//制定阶
	void moments(const Mat& src, Mat& t, Mat& tnp, int N, int M,RECURSION_TYPE rtype);
	//重建
	void reconstruct(const Mat& moments, const Mat& tnp, Mat& dst,int N);
	//近似重建
	void reconstruct(const Mat& moments, const Mat& tnp, Mat& dst, int N, int M);
	//重建误差
	double error(const Mat& src, const Mat& dst);

	void poly(Mat& tnp, int N, RECURSION_TYPE rtype);
};

