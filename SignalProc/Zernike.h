#pragma once
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;
using namespace cv;

class Zernike
{
public:
	Zernike();
	~Zernike();

	void process(const Mat& src, Mat& moments, Mat& cp, Mat& sp,map< pair<int, int>, vector<double> >& rnm, int nmax, int N);
	//����任
	void transform(const Mat& src, Mat& dst);

	void transform_count(const Mat& src, Mat& dst);
	//���귴�任
	void inverse(const Mat& dst, Mat& src,int x,int y);
	//�������ʽ
	void fillRNM(map< pair<int, int>, vector<double> >& rnm, int nmax, int N);
	//�ؽ�
	void reconstruct(map< pair<int, int>, vector<double> >& rnm, const Mat& moments, const Mat& cnm, const Mat& snm, Mat& dst, int nmax, int N);

private:
	vector< pair<pair<int,int>,pair<int,int>> > point_pair;
};

