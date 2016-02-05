#include "Legendre.h"


Legendre::Legendre()
{
}


Legendre::~Legendre()
{

}

void Legendre::process(const Mat& src, Mat& moments,int p,int q)
{

}

void Legendre::reconstruct(const Mat& src, Mat& dst,int p, int q,int N)
{

}

void Legendre::rawPorcess(const Mat& src, Mat& moments,int p, int q)
{
	if (src.rows != src.cols)
	{
		cout << "请使用方形图片" << endl;
		return;
	}

	//int n = p, m = q; //矩形的阶 通常N==M
	int N = src.cols; //原始图片的大小 NxN
	double a = 0, all = 0; //累加后的系数 累加后的结果
	moments = Mat(p, q, CV_32FC1, Scalar(0)); //矩
	double ele = 0, sr = 0;
	double xi = 0, xj = 0;
	for (int n = 0; n < moments.rows; n++)
	{
		for (int m = 0; m < moments.cols; m++)
		{
			sr = 0;
			for (int i = 1; i <= N;i++)
			{
				for (int j = 1; j <= N; j++)
				{
					xi = (2 * i - N - 1) / (float)(N - 1);
					xj = (2 * j - N - 1) / (float)(N - 1);
					sr += src.at<uchar>(i-1, j-1)*pnx(n, xi)*pnx(m, xj);
				}
			}
			a = (2 * n + 1)*(2 * m + 1) / (float)pow(N - 1, 2); //归一化系数
			moments.at<float>(n, m) = ele*a;
		}
	}
}


double Legendre::pnx(int n, float x)
{
	//TODO
	float m = 1;
	double sr = 0;
	for (int k = 0; k <= n; k++) //< 还是 <=
	{
		if ((n - k) % 2 == 1) continue;
		int nmk = (n - k) / 2;
		int npk = (n + k) / 2;
		int a = nmk % 2 == 1 ? -1 : 1;
		double cnk= a*Tool::multiTo(n + k) / (float)(Tool::multiTo(nmk)*Tool::multiTo(npk)*Tool::multiTo(k))/pow(2,n);
		sr += cnk*pow(x, k);
	}
	return sr;
}


void  Legendre::rawReconstruct(const Mat& src, Mat& dst, int p, int q,int N)
{
	double ele = 0, sr = 0;
	double xi = 0, xj = 0;
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			sr = 0;
			for (int n = 0; n <src.rows; n++)
			{
				for (int m = 0; m < src.cols; m++)
				{
					xi = (2 * i - N - 1) / (float)(N - 1);
					xj = (2 * j - N - 1) / (float)(N - 1);
					sr += src.at<float>(n, m)*pnx(n, xi)*pnx(m, xj);
				}
			}
			dst.at<uchar>(i-1, j-1) = cvRound(sr)>255 ? 255 : cvRound(sr);
		}
	}
}

void Legendre::getPoly(Mat& poly,int order, int N)
{
	poly=Mat(order + 1, N+1, CV_64FC1);
	double x;
	for (int i = 0; i <= order; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			if (i == 0) 
			{
				poly.at<double>(i, j) = 1;
				continue;
			}
			x = (2 * j - N - 1) / (double)(N - 1);
			if (i == 1)
			{
				poly.at<double>(i, j) = x;
			}
			else
			{
				poly.at<double>(i, j) = ((2 * i - 1)*x*poly.at<double>(i - 1, j) - (i - 1)*poly.at<double>(i - 2, j)) / (double)i;
			}
		}
	}
}

//points: x1j,x2j,yj
void Legendre::mukundanProcess(const Mat& src, const Mat& points, Mat& moments, int order, int N)
{
	Mat poly;
	getPoly(poly,order, N);
	double s, tjp, tj;
	int x1j, x2j, yj;
	moments = Mat(order + 1, order + 1, CV_64FC1);
	for (int n = 0; n <= order; n++)
	{
		for (int m = 0; m <= order; m++)
		{
			s = 0;
			for (int j = 1; j <= points.rows; j++)
			{
				Point3i p=points.at<Point3i>(j);
				x1j = p.x; x2j = p.y; yj = p.z;
				tj = poly.at<double>(m, yj)*(x2j*poly.at<double>(n, x2j) - poly.at<double>(n - 1, x2j) - x1j*poly.at<double>(n, x1j) + poly.at<double>(n - 1, x1j));
				if (j > 1) s += (tj + tjp);
				tjp = tj;
			}
			moments.at<double>(n, m) = s*(2 * n + 1)*(2 * m + 1) / (4 * (n + 1)*(N - 1));
		}
	}
}

void Legendre::mukundanReconstruct(const Mat& moments, Mat& dst, int order, int N)
{
	double ele = 0, sr = 0;
	double xi = 0, xj = 0;
	Mat poly;
	getPoly(poly, order, N);
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			sr = 0;
			for (int n = 0; n <moments.rows; n++)
			{
				for (int m = 0; m < moments.cols; m++)
				{
					sr += moments.at<float>(n, m)*poly.at<double>(n, i)*poly.at<double>(m, j);
				}
			}
			dst.at<uchar>(i - 1, j - 1) = cvRound(sr)>255 ? 255 : cvRound(sr);
		}
	}
}