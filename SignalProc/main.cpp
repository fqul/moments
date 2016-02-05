#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "Legendre.h"
#include "Tchebichef.h"
#include "Zernike.h"

using namespace cv;
using namespace std;

void resize()
{
	for (char a = 'a'; a <= 'z'; a++)
	{
		string filename;
		stringstream ss(filename);
		ss << "characters\\" <<(char)a<< ".bmp";
		Mat src = imread(ss.str(), CV_LOAD_IMAGE_GRAYSCALE);
		if (src.data == NULL)
		{
			cout << "open file error: " << ss.str() << endl;
			continue;
		}

		string destFileName;
		stringstream sm(filename);
		sm << "characters\\" << a << "_60x60.bmp";
		cout << sm.str() << endl;
		imwrite(sm.str(), src(Range(0, 60), Range(0, 60)));
	}
}

void resize_bird()
{
	string filename = "bird.bmp";
	Mat src = imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
	imwrite("bird_172x172.bmp", src(Range(0, 172), Range(0, 172)));
}

void testCharacters()
{
	for (int i = 'a'; i <= 'z';i++)
	{
		string sr;
		stringstream sp(sr);
		sp << "characters\\" << (char)i << "_60x60.bmp";
		cout << sp.str() << endl;
		Mat src = cv::imread(sp.str(), CV_LOAD_IMAGE_GRAYSCALE);
		if (src.data == NULL)
		{
			cout << "open image file error" << endl;
			return;
		}
		cout << src.size() << endl;
		Tchebichef chebi;
		Mat moments, tnp;
		RECURSION_TYPE tt = RECURSION_TYPE::BY_N;
		int moment = 59;
		chebi.moments(src, moments, tnp, src.rows, moment, tt);
		cout << "get moments" << endl;

		Mat reconstruted_app;
		chebi.reconstruct(moments, tnp, reconstruted_app, src.rows, moment);
		cout << "get reconstructed image app" << endl;

		double e_rate = chebi.error(src, reconstruted_app);
		cout << "错误率：" << e_rate << endl;


		string ss;
		stringstream so(ss);
		so << i << "_60.bmp";
		cout << so.str() << endl;
		if (imwrite(so.str(), reconstruted_app))
		{
			cout << "写入成功" << endl;
		}
		else{
			cout << "写入失败" << endl;
		}
	}
}

void testTchebichef()
{
	//Mat src = cv::imread("characters\\c_60x60.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	Mat src = cv::imread("lena.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	if (src.data == NULL)
	{
		cout << "open image file error" << endl;
		return;
	}

	cout << src.size() << endl;
	Tchebichef chebi;
	Mat moments, tnp;
	//int moment = 80;
	for (int moment = 160; moment > 0;moment-=20)
	{
		//chebi.moments(src, moments, tnp, src.rows, RECURSION_TYPE::BY_N);
		RECURSION_TYPE tt = RECURSION_TYPE::BY_N;

		chebi.moments(src, moments, tnp, src.rows, moment, tt);
		cout << "get moments" << endl;

		//Mat reconstruted;
		//chebi.reconstruct(moments, tnp, reconstruted, src.rows);
		//cout << "get reconstructed image" << endl;


		Mat reconstruted_app;
		chebi.reconstruct(moments, tnp, reconstruted_app, src.rows, moment);
		cout << "get reconstructed image app" << endl;

		double e_rate = chebi.error(src, reconstruted_app);
		cout << "错误率：" << e_rate << endl;


		string ss;
		stringstream so(ss);
		so << "bird_" << moment << "_" << tt << "_" << e_rate << ".bmp";
		cout << so.str() << endl;
		if (imwrite(so.str(), reconstruted_app))
		{
			cout << "写入成功" << endl;
		}
		else{
			cout << "写入失败" << endl;
		}
	}
	//imshow("src", src);
	//imshow("reconstructed", reconstruted);
	//imshow("reconstructed app", reconstruted_app);

	//waitKey();
	//destroyAllWindows();
}

void testZernike()
{
	Mat src = imread("characters\\a_60x60.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	if (src.data == NULL)
	{
		cout << "image read error" << endl;
		return;
	}
	Zernike nike;
	Mat dst;
	//Mat dst, moments, cp, sp;
	nike.transform(src, dst);
	//map<pair<int, int>, vector<double> > rnm;
	//nike.process(src, moments, cp, sp, rnm, 30, 60);

	Mat src1;
	nike.inverse(dst, src1, 172, 172);

	Mat dst2, src2;
	nike.transform_count(src, dst2);
	nike.inverse(dst2, src2, 172, 172);

	Mat moments, cp, sp;
	map< pair<int, int>, vector<double> > rnm;
	nike.process(src, moments, cp, sp, rnm, 50, 60);
	Mat reconst;
	nike.reconstruct(rnm, moments, cp, sp, reconst, 50, 60);
	imshow("src2", src2);
	imshow("src1", src1);
	imshow("src", src);
	imshow("dst", dst);
	imshow("dst2", dst2);
	imshow("reconst", reconst);
	waitKey(0);
	destroyAllWindows();
}

void testSTL()
{
	vector<int> name(20, -1);
	name[10] = 10;
	for (vector<int>::iterator iter = name.begin(); iter != name.end(); iter++)
	{
		cout << *iter << endl;
	}

	cout << "--------------------" << endl;
	map<string, int> words;
	words["hello"] = 20;
	cout << words["hello"] << endl;
}

void testMat()
{
	Mat mat(20, 20, CV_64FC3, Scalar(0, 0, 0));
	Scalar s(1, 2, 3);
	cout << s[0] << " " << s[1] << endl;
	mat.at<Scalar>(1, 1)[1] = 20;
	cout << mat.at<Scalar>(1,1)[1] << endl;
}



int main(int argc, char** argv)
{
	int m;
	//testZernike();
	//testSTL();
	//testMat();
	//Mat image = imread("lena-gray.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	//resize(image, image, Size(256, 256));
	//imwrite("lena.bmp", image);
	//testTchebichef();
	testCharacters();
	cin >> m;
	return 0;
}

