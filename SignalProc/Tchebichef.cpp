#include "Tchebichef.h"

Tchebichef::Tchebichef()
{

}

Tchebichef::~Tchebichef()
{

}


void Tchebichef::moments(const Mat& src, Mat& t, Mat& tnp, int N,int M, RECURSION_TYPE rtype)
{
	poly(tnp, N, rtype);
	Mat tn = tnp;
	t = Mat(N, N, CV_64FC1);
	double s;
	for (int m = 0; m <=M; m++)
	{
		for (int n = 0; n <=m; n++)
		{
			s = 0;
			for (int x = 0; x < N; x++)
			{
				for (int y = 0; y < N; y++)
				{
					s += tn.at<double>(m - n, x)*tn.at<double>(n, y)*src.at<uchar>(x, y);
				}
			}
			t.at<double>(m - n, n) = s;
		}
	}
}

void Tchebichef::moments(const Mat& src, Mat& t, Mat& tnp, int N, RECURSION_TYPE rtype)
{
	poly(tnp, N, rtype);
	Mat tn = tnp;
	t = Mat(N, N, CV_64FC1);
	double s;
	for (int n = 0; n < t.rows; n++)
	{
		for (int m = 0; m < t.cols; m++)
		{
			s = 0;
			for (int x = 0; x < N; x++)
			{
				for (int y = 0; y < N; y++)
				{
					s += tn.at<double>(n, x)*tn.at<double>(m, y)*src.at<uchar>(x, y);
				}
			}
			t.at<double>(n, m) = s;
		}
	}
}

double Tchebichef::error(const Mat& src, const Mat& dst)
{
	if (src.size() != dst.size())
	{
		cout << "size error" << endl;
		return 0.0;
	}
	if (src.type()!=dst.type())
	{
		cout << "type doesn't match" << endl;
		cout << "src type: " << src.type() << "; dst type: " << dst.type() << endl;
		return 0.0;
	}

	int C = src.rows*src.cols;
	Mat m=src - dst;
	pow(m, 2, m);
	
	double result = sum(m).val[0];
	result /= C;
	return result;
}

void Tchebichef::reconstruct(const Mat& moments, const Mat& tnp, Mat& dst, int N)
{
	double sm;
	dst = Mat(N, N, CV_8UC1);
	for (int x = 0; x < N; x++)
	{
		for (int y = 0; y < N; y++)
		{
			sm = 0;
			for (int m = 0; m < N;m++)
			{
				for (int n = 0; n < N;n++)
				{
					sm += moments.at<double>(m, n)*tnp.at<double>(m, x)*tnp.at<double>(n, y);
				}
			}
			dst.at<uchar>(x, y) = sm>255 ? 255 : cvFloor(sm);
		}
	}
}

void Tchebichef::reconstruct(const Mat& moments, const Mat& tnp, Mat& dst, int N, int M)
{
	double sm;
	dst = Mat(N, N, CV_8UC1);
	for (int x = 0; x < N; x++)
	{
		for (int y = 0; y < N; y++)
		{
			sm = 0;
			for (int m = 0; m <=M; m++)
			{
				for (int n = 0; n <=m; n++)
				{
					sm += moments.at<double>(m-n, n)*tnp.at<double>(m-n, x)*tnp.at<double>(n, y);
				}
			}
			dst.at<uchar>(x, y) = sm>255 ? 255 : cvRound(sm);
		}
	}
}


void Tchebichef::poly(Mat& tnp, int N, RECURSION_TYPE rtype)
{
	Mat tn = Mat(N, N, CV_64FC1, Scalar::all(0));
	if (rtype == RECURSION_TYPE::BY_X)
	{
		double a, b;
		tn.at<double>(0, 0) = 1.0 / (double)sqrt(N);
		for (int n = 1; n < N; n++)
		{
			tn.at<double>(n, 0) = -1 * sqrt((N - n) / (double)(N + n))*sqrt((2 * n + 1) / (double)(2 * n - 1))*tn.at<double>(n - 1, 0);
			tn.at<double>(n, N - 1) = (n % 2 == 0 ? 1 : -1)*tn.at<double>(n, 0);

		}

		for (int n = 0; n < N ; n++)
		{
			tn.at<double>(n, 1) = (1 + n*(1 + n) / (double)(1 - N))*tn.at<double>(n, 0);
			tn.at<double>(n, N - 2) = (n % 2 == 0 ? 1 : -1)*tn.at<double>(n, 1);
		}

		for (int n = 0; n < N; n++)
		{
			for (int x = 2; x < N/2; x++)
			{
				a = (-1 * n*(n + 1) - (2 * x - 1)*(x - N - 1) - x) / (double)(x*(N - x));
				b = (x - 1)*(x - N - 1) / (double)(x*(N - x));
				tn.at<double>(n, x) = a*tn.at<double>(n, x - 1) + b*tn.at<double>(n, x - 2);
				tn.at<double>(n, N - 1 - x) = (n % 2 == 0 ? 1 : -1)*tn.at<double>(n, x);
			}
		}
	}
	else if (rtype == RECURSION_TYPE::BY_N)
	{
		double a, b;
		//œ‡Õ¨X£¨µ›Õ∆N
		for (int x = 0; x < N/2; x++)
		{
			for (int n = 0; n < 2; n++)
			{
				if (n == 0)
				{
					tn.at<double>(n, x) = 1.0 / sqrt(N);
					tn.at<double>(n, N - 1 - x) = (n % 2 == 0 ? 1 : -1)*tn.at<double>(n, x);
				}
				else if (n == 1)
				{
					tn.at<double>(n, x) = (2 * x + 1 - N)*sqrt(3.0 / (N*(N*N - 1)));
					tn.at<double>(n, N - 1 - x) = (n % 2 == 0 ? 1 : -1)*tn.at<double>(n, x);
				}
			}
		}

		for (int x = 0; x < N/2; x++)
		{
			for (int n = 2; n < N; n++)
			{
				a = ((2 * x - N + 1) / (double)n)*sqrt((4 * n*n - 1) / (double)(N*N - n*n));
				b = -1 * ((n - 1) / (double)n)*sqrt((2 * n + 1) / (double)(2 * n - 3))*sqrt((N*N - (n - 1)*(n - 1)) / (double)(N*N - n*n));
				tn.at<double>(n, x) = a*tn.at<double>(n - 1, x) + b*tn.at<double>(n - 2, x);
				tn.at<double>(n, N - 1 - x) = (n % 2 == 0 ? 1 : -1)*tn.at<double>(n, x);
			}
		}
	}
	tnp = tn;
}