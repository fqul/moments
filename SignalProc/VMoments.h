#pragma once
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>

#include <iostream>
#include <vector>

using namespace std;
using namespace cv;

class VMoments
{
public:
	VMoments();
	~VMoments();
	void polar();
	void transform(const Mat& src, Mat& dst, int U, int M);
	double hx(double x);
};

