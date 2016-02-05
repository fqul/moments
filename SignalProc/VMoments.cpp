#include "VMoments.h"


VMoments::VMoments()
{
}


VMoments::~VMoments()
{
}

double VMoments::hx(double x)
{
	double b = abs(x);
	if (b<= 1)
	{
		return 1.5*pow(b, 3) - 2.5*pow(x, 2) + 1;
	}
	else if (b <= 2)
	{
		return -0.5*pow(b, 3) + 2.5*pow(b, 2) - 4 * b + 2;
	}
	else
	{
		return 0;
	}
}

void VMoments::transform(const Mat& src, Mat& dst, int U,int M)
{
	dst = Mat(M, M, CV_32FC3, Scalar(0, 0, 0));
	int bin = U / 2;
	float s = 0;
	for (int m = 0; m < M;m++)
	{
		for (int n = 0; n < m;n++)
		{
			s = 0;
			for (int u = 0; u < U; u++)
			{
				for (int v = 0; v < 4 * (2u - 1); v++)
				{
				}
			}
		}
	}
}
