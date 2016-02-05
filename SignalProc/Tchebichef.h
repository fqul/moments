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
	//��أ�
	//src: ����ε�ͼƬ��t����žأ�NͼƬ�Ĵ�С����˾صĽ�ΪN-1+N-1��rtype������x���ƣ����߰���n����
	void moments(const Mat& src,Mat& t,Mat& tnp,int N,RECURSION_TYPE rtype);

	//�ƶ���
	void moments(const Mat& src, Mat& t, Mat& tnp, int N, int M,RECURSION_TYPE rtype);
	//�ؽ�
	void reconstruct(const Mat& moments, const Mat& tnp, Mat& dst,int N);
	//�����ؽ�
	void reconstruct(const Mat& moments, const Mat& tnp, Mat& dst, int N, int M);
	//�ؽ����
	double error(const Mat& src, const Mat& dst);

	void poly(Mat& tnp, int N, RECURSION_TYPE rtype);
};

