#include "Zernike.h"


Zernike::Zernike()
{
}


Zernike::~Zernike()
{
}

void Zernike::inverse(const Mat& dst, Mat& src,int x,int y)
{
	src = Mat(x, y, CV_8UC1, Scalar::all(0));
	typedef vector < pair<pair<int, int>, pair<int, int> > > vppp;
	typedef pair<pair<int, int>, pair<int, int>> ppp;
	int p, q, i, j;
	for (vppp::iterator iter = point_pair.begin(); iter != point_pair.end(); iter++)
	{
		ppp ele = *iter;
		p = ele.first.first; q = ele.first.second; i = ele.second.first; j = ele.second.second;
		//cout << p << " " << q << "; " << i << " " << j << endl;
		src.at<uchar>(i, j) = dst.at<uchar>(p, q);
	}
}

void Zernike::fillRNM(map< pair<int, int>, vector<double> >& rnm,int nmax,int N)
{
	int m, s;
	double r, b0, brs;
	vector<double> pp(N+1, 0);

	for (int n = 0; n <=nmax; n++)
	{
		for (int m = 0; m <= n; m++)
		{
			//if ((n - m) % 2 == 0)
			//{
			rnm[make_pair(n, m)] = vector<double>(N + 1, 0);
			//}
		}
	}

	for (int p = 1; p <= N; p++)
	{
		r = 2 * p / (double)N;
		for (int n = 0; n <= nmax; n++)
		{
			m = n;
			s = n;
			b0 = pow(r, s);
			brs = b0;
			while (true)
			{
				rnm[make_pair(n, m)][p] += brs;
				if (s > m + 1)
				{
					brs *= -1 * (s*s - m*m) / (double)((s + n)*(n - s + 2));
					s -= 2;
				}
				else
				{
					b0 *= (n + m) / (double)(n - m + 2);
					if (m > 1)
					{
						m -= 2;
						s = n;
						brs = b0;
					}
					else
					{
						break;
					}
				}
			}
		}
	}
}

void Zernike::process(const Mat& src, Mat& moments, Mat& cp, Mat& sp, map< pair<int, int>, vector<double> >& rnmp, int nmax, int N)
{
	double r, factor, th;
	map< pair<int, int>, vector<double> > rnm; fillRNM(rnm, nmax, N);
	Mat c(nmax + 1, nmax + 1, CV_64FC1);
	Mat s(nmax + 1, nmax + 1, CV_64FC1);
	Mat m(nmax + 1, nmax + 1, CV_64FC1);
	Mat f; transform_count(src, f);
	for (int n = 0; n <= nmax; n++)
	{
		for (int p = 1; p <= N/2; p++)
		{
			r = 2 * p / (double)N;
			for (int q = 1; q <= 8 * p; q++)
			{
				th = CV_PI*q / (4 * p);
				for (int m = 1; m <= n; m++)
				{
					if ((n - m) % 2 == 0)
					{
						//cout << n << " " << m << ";" << p << " " << q << endl;
						c.at<double>(n, m) += rnm[make_pair(n, m)][p] * cos(m*th)*f.at<uchar>(p, q);
						s.at<double>(n, m) += rnm[make_pair(n, m)][p] * sin(m*th)*f.at<uchar>(p, q);
					}
				}
			}
		}

		factor = (2 * n + 2) / (double)(N*N);
		for (int m = 1; m <= n; m++)
		{
			if ((n - m) % 2 == 0)
			{
				c.at<double>(n, m) *= factor;
				s.at<double>(n, m) *= factor;
			}
		}
	}

	c.copyTo(cp);
	s.copyTo(sp);
	m = c + s;
	m.copyTo(moments);
	rnmp = rnm;
}



void Zernike::transform_count(const Mat& src, Mat& dst)
{
	assert(src.rows == src.cols);
	int N = src.rows; //为偶数
	point_pair.clear();
	dst = Mat(N / 2 + 1, 4 * N + 1, CV_8UC1, Scalar::all(0));
	int bin = N / 2;
	int x, y;
	int a, b;
	int dx, dy;
	int count;
	for (int p = 1; p<=bin ; p++)
	{
		count = 0;
		x = p; y = 0; dx = 0; dy = 1;
		bool first = true;
		while (!(x == p && y == 0) || first)
		{
			if (x*y != 0)
			{
				//平移 补空
				a = x; b = y;
				if (!(x < 0 && y>0))
				{
					if (y < 0) b = y + 1; 
					if (x > 0) a = x - 1;
				}
				b = -b; swap(a, b); a += bin; b += bin;
				count++;
				//cout << p << " " << count << "; " << a << " " << b << "; " << x << " " << y << endl;
				assert(a < N && b < N);
				dst.at<uchar>(p, count) = src.at<uchar>(a, b);
				point_pair.push_back(make_pair(make_pair(p, count), make_pair(a, b)));

			}
			//跳转到下一个位置；
			x += dx; y += dy;
			//调整方向
			if (dx == 0 && dy == 1 && y == p) { dx = -1; dy = 0; }
			if (dx == -1 && dy == 0 && x == -p) { dx = 0; dy = -1; }
			if (dx == 0 && dy == -1 && y == -p) { dx = 1; dy = 0; }
			if (dx == 1 && dy == 0 && x == p) { dx = 0; dy = 1; }
			if (first) first = false;
		}
	}
}

void Zernike::transform(const Mat& src, Mat& dst)
{
	assert(src.rows == src.cols);
	point_pair.clear();
	int N = src.rows; //得保证是偶数
	dst = Mat(N / 2 + 1, 4 * N + 1, CV_8UC1, Scalar::all(0));
	int bin = N / 2;
	int isEven = N % 2 == 0 ? true : false;
	int x, y;
	int p, q;
	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			if (isEven)
			{
				if (i >= bin)
					x = i - bin + 1;
				else
					x = i - bin;
				if (j >= bin)
					y = j - bin + 1;
				else
					y = j - bin;
			}
			else
			{
				//暂时不考虑
			}
			swap(x, y);
			y = -y;
			p = max(abs(x), abs(y));
			if (abs(x) == p)
				q = 2 * (p - x)*y / abs(y) + x*y / p;
			if (abs(y) == p)
				q = 2 * y - x*y / p;

			//q = abs(q);
			if (q < 0)
			{
				//cout << "<0: " << p << " " << q << "; " << i << " " << j << "; " << x << " " << y << endl;
				//continue;
				q += 4 * (2 * p - 1)+1;
			}
			//q = abs(q);
			//cout << p << " " << q << "; " << i << " " << j << "; " << x << " " << y << endl;
			dst.at<uchar>(p, q) = src.at<uchar>(i, j);
			point_pair.push_back(make_pair(make_pair(p, q), make_pair(i, j)));
		}
	}
}

void Zernike::reconstruct(map< pair<int, int>, vector<double> >& rnm, const Mat& moments, const Mat& cnm, const Mat& snm, Mat& dst, int nmax, int N)
{
	Mat pqm = Mat(N / 2 + 1, 4 * N + 1, CV_8UC1, Scalar::all(0));
	double s = 0, th;
	for (int p = 1; p <= N/2; p++)
	{
		for (int q = 1; q <= 8 * p; q++)
		{
			s = 0;
			th = CV_PI*q / (4 * p);
			for (int n = 0; n <= nmax;n++)
			{
				for (int m = 1; m <= n;m++)
				{
					if ((n - m) % 2 == 0)
					{
						s += (cnm.at<double>(n, m)*cos(m*th) + snm.at<double>(n, m)*sin(m*th))*rnm[make_pair(n, m)][p];
					}
				}
			}
			s = s*N / 2;
			pqm.at<uchar>(p, q) = s > 255 ? 255 : cvFloor(s);
		}
	}

	inverse(pqm, dst, N, N);
}