#include <opencv2/opencv.hpp>
#include <cv.h>
#include "SubImageMatch.h"

using namespace std;
using namespace cv;

inline float mySqrt(float x)
{
	float a = x;
	unsigned int i = *(unsigned int *)&x;
	i = (i + 0x3f76cf62) >> 1;
	x = *(float *)&i;
	x = (x + a / x) * 0.5f;
	return x;
}

int ustc_ConvertBgr2Gray(Mat bgrImg, Mat& grayImg)
{
//	int a = 114, b = 587, c = 299;
	int i = 0, j = 0, row = bgrImg.rows, col = bgrImg.cols;
	grayImg = Mat(row, col, CV_8UC1, 1);
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			grayImg.at<uchar>(i, j) = bgrImg.at<Vec3b>(i, j)[0] * 114 / 1000 + bgrImg.at<Vec3b>(i, j)[1] * 587 / 1000 + bgrImg.at<Vec3b>(i, j)[2] * 299 / 1000;
		}
	}
	return 1;
}

int ustc_CalcGrad(Mat grayImg, Mat& gradImg_x, Mat& gradImg_y)
{
	int i = 0, j = 0, row = grayImg.rows, col = grayImg.cols;
	grayImg.convertTo(grayImg, CV_32FC1, 1.0 / 255);
	gradImg_x = Mat(row, col, CV_32FC1);
	gradImg_y = Mat(row, col, CV_32FC1);
	row--;
	col--;
	for (i = 1;i < row;i++)
	{
		for (j = 1;j < col;j++)
		{
			gradImg_x.at<float>(i, j) = grayImg.at<float>(i - 1, j + 1) + 2 * grayImg.at<float>(i, j + 1) + grayImg.at<float>(i + 1, j + 1) - grayImg.at<float>(i - 1, j - 1) - 2 * grayImg.at<float>(i, j - 1) - grayImg.at<float>(i + 1, j - 1);
			gradImg_y.at<float>(i, j) = grayImg.at<float>(i + 1, j - 1) + 2 * grayImg.at<float>(i + 1, j) + grayImg.at<float>(i + 1, j + 1) - grayImg.at<float>(i - 1, j - 1) - 2 * grayImg.at<float>(i - 1, j) - grayImg.at<float>(i - 1, j + 1);
		}
	}
	return 1;
}

int ustc_CalcAngleMag(Mat gradImg_x, Mat gradImg_y, Mat& angleImg, Mat& magImg)
{
	int i = 0, j = 0, row = gradImg_x.rows, col = gradImg_x.cols;
	float temp;
	gradImg_x.convertTo(gradImg_x, CV_8UC1, 255);
	gradImg_y.convertTo(gradImg_y, CV_8UC1, 255);
	angleImg = Mat(row, col, CV_32FC1);
	magImg = Mat(row, col, CV_8UC1);
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			magImg.at<uchar>(i, j) = mySqrt(gradImg_x.at<uchar>(i, j)*gradImg_x.at<uchar>(i, j) + gradImg_y.at<uchar>(i, j)*gradImg_x.at<uchar>(i, j));
			if (gradImg_x.at<uchar>(i, j) == 0)
			{
				angleImg.at<float>(i, j) = 90.0;
				continue;
			}
			temp = float(gradImg_y.at<uchar>(i, j)) / gradImg_x.at<uchar>(i, j);
			if (temp <= 0.5)
			{
				angleImg.at<float>(i, j) = (temp*(1.0f + temp*temp*(-1.0f / 3 + temp*temp*(1.0f / 5 - temp*temp / 7))))*57.295f;
			}
			else if (temp <= 2.0)
			{
				temp--;
				angleImg.at<float>(i, j) = (3.14159f / 4 + temp*(1.0f / 2 + temp*(-1.0f / 4 + temp*(1.0f / 12 + temp*temp*(1.0f / 40 + temp / 48)))))*57.295f;
			}
			else
			{
				angleImg.at<float>(i, j) = (3.14159f / 2 + (((((1.0f / (11 * temp*temp) - 1.0f / 9) / (temp*temp) + 1.0f / 7) / (temp*temp) - 1.0f / 5) / (temp*temp) + 1.0f / 3) / (temp*temp) - 1.0f) / temp)*57.295f;
			}
		}
	}
	return 1;
}

int ustc_Threshold(Mat grayImg, Mat& binaryImg, int th)
{
	int i = 0, j = 0, row = grayImg.rows, col = grayImg.cols;
	binaryImg = Mat(row, col, CV_8UC1, 1);
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			binaryImg.at<uchar>(i, j) = (grayImg.at<uchar>(i, j) > th) * 255;
		}
	}
	return 1;
}

int ustc_CalcHist(Mat grayImg, int* hist, int hist_len)
{
	int i = 0, j = 0, row = grayImg.rows, col = grayImg.cols;
	memset(hist, 0, hist_len*sizeof(int));
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			hist[grayImg.at<uchar>(i, j) * 256 / hist_len]++;
		}
	}
	return 1;
}

int ustc_SubImgMatch_gray(Mat grayImg, Mat subImg, int* x, int* y)
{
	int i = 0, j = 0, k = 0, l = 0, m = 0;
	int mx = 0, my = 0;
	uchar *pg(NULL), *ps(NULL);
	long long temp = 0;
	long long ans = 0xfffffffffffffff;
	int row = grayImg.rows, col = grayImg.cols;
	int rrr = subImg.rows, ccc = subImg.cols;
	row -= rrr;
	col -= ccc;
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			temp = 0;
			for (k = 0;k < rrr;k++)
			{
				pg = grayImg.ptr<uchar>(i + k);
				ps = subImg.ptr<uchar>(k);
				for (l = 0;l < ccc;l++)
				{
					m = pg[j + l] - ps[l];
					temp += m > 0 ? m : -m;
				}
			}
			if (temp < ans)
			{
				ans = temp;
				mx = i, my = j;
			}
		}
	}
	*x = mx;
	*y = my;
	return 1;
}

int ustc_SubImgMatch_bgr(Mat colorImg, Mat subImg, int* x, int* y)
{
	int i = 0, j = 0, k = 0, l = 0, m = 0;
	int mx = 0, my = 0;
	Vec3b *pc(NULL),*ps(NULL);
	long long temp = 0;
	long long ans = 0x7fffffffffffffff;
	int row = colorImg.rows, col = colorImg.cols;
	int rrr = subImg.rows, ccc = subImg.cols;
	row -= rrr;
	col -= ccc;
	for (i = 0;i < row;i++)
	{
//		cout << i << endl;
		for (j = 0;j < col;j++)
		{
			temp = 0;
			for (k = 0;k < rrr;k++)
			{
				pc = colorImg.ptr<Vec3b>(i + k);
				ps = subImg.ptr<Vec3b>(k);
				for (l = 0;l < ccc;l++)
				{
					m = pc[j + l][0] - ps[l][0];
					temp += m > 0 ? m : -m;
					m = pc[j + l][1] - ps[l][1];
					temp += m > 0 ? m : -m;
					m = pc[j + l][2] - ps[l][2];
					temp += m > 0 ? m : -m;
				}
			}
			if (temp < ans)
			{
				ans = temp;
				mx = i, my = j;
			}
		}
	}
	*x = mx;
	*y = my;
	return 1;
}

int ustc_SubImgMatch_corr(Mat grayImg, Mat subImg, int* x, int* y)
{
	int i, j, k, l, m, n;
	subImg.convertTo(subImg, CV_32FC1, 1.0 / 255);
	grayImg.convertTo(grayImg, CV_32FC1, 1.0 / 255);
	int grayrow = grayImg.rows, graycol = grayImg.cols;
	int subrow = subImg.rows, subcol = subImg.cols;
	Mat graysq = grayImg.clone(), subsq = subImg.clone();
	float *gs, *qs;
	float temp, ans=0.0, gg, ss;
	int mx = 255, my = 255;
	for (i = 0;i < grayrow;i++)
	{
		gs = graysq.ptr<float>(i);
		for (j = 0;j < graycol;j++)
		{
			gs[j] = gs[j] * gs[j];
		}
	}
	for (i = 0;i < subrow;i++)
	{
		gs = subsq.ptr<float>(i);
		for (j = 0;j < subcol;j++)
		{
			gs[j] = gs[j] * gs[j];
		}
	}
	grayrow -= subrow;
	graycol -= subcol;
	for (i = 0;i < grayrow;i++)
	{
		for (j = 0;j < graycol;j++)
		{
			temp = gg = ss = 0.0;
			for (k = 0;k < subrow;k++)
			{
				gs = graysq.ptr<float>(i + k);
				qs = subsq.ptr<float>(k);
				for (l = 0;l < subcol;l++)
				{
					temp += grayImg.at<float>(i + k, j + l)*subImg.at<float>(k, l);
					gg += gs[j+l];
					ss += qs[l];
				}
			}
			temp = temp*temp;
			temp /= (ss*gg);
			if (temp >= ans)
			{
				mx = i;
				my = j;
				ans = temp;
			}
		}
	}
	*x = mx;
	*y = my;
	return 1;
}

int ustc_SubImgMatch_angle(Mat grayImg, Mat subImg, int* x, int* y)
{
	int i = 0, j = 0, k = 0, l = 0;
	int mx = 0, my = 0;
	int row = grayImg.rows, col = grayImg.cols;
	int rrr = subImg.rows, ccc = subImg.cols;
	row -= rrr;
	col -= ccc;
	float *pg(NULL), *ps(NULL);
	float m, temp, ans=0;
	Mat gangle, gmag, sangle, smag;
	Mat ggradx, ggrady, sgradx, sgrady;
	bool flag = true;
	ustc_CalcGrad(grayImg, ggradx, ggrady);
	ustc_CalcGrad(subImg, sgradx, sgrady);
	ustc_CalcAngleMag(ggradx, ggrady, gangle, gmag);
	ustc_CalcAngleMag(sgradx, sgrady, sangle, smag);
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			temp = 0.0;
			for (k = 0;k < rrr;k++)
			{
				pg = gangle.ptr<float>(i + k);
				ps = sangle.ptr<float>(k);
				for (l = 0;l < ccc;l++)
				{
					m = pg[j + l] - ps[l];
					temp += m > 0 ? m : -m;
				}
			}
//			cout << temp << endl;
			if (temp < ans || flag)
			{
				ans = temp;
				mx = i, my = j;
				flag = false;
			}
		}
	}
	*x = mx;
	*y = my;
	return 1;
}

int ustc_SubImgMatch_mag(Mat grayImg, Mat subImg, int* x, int* y)
{
	int i = 0, j = 0, k = 0, l = 0;
	int mx = 0, my = 0;
	int row = grayImg.rows, col = grayImg.cols;
	int rrr = subImg.rows, ccc = subImg.cols;
	row -= rrr;
	col -= ccc;
	long long temp, ans = 0, m = 0;
	bool flag = true;
	uchar *pg(NULL), *ps(NULL);
	Mat gangle, gmag, sangle, smag;
	Mat ggradx, ggrady, sgradx, sgrady;
	ustc_CalcGrad(grayImg, ggradx, ggrady);
	ustc_CalcGrad(subImg, sgradx, sgrady);
	ustc_CalcAngleMag(ggradx, ggrady, gangle, gmag);
	ustc_CalcAngleMag(sgradx, sgrady, sangle, smag);
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			temp = 0.0;
			for (k = 0;k < rrr;k++)
			{
				pg = gmag.ptr<uchar>(i + k);
				ps = smag.ptr<uchar>(k);
				for (l = 0;l < ccc;l++)
				{
					m = pg[j + l] - ps[l];
					temp += m > 0 ? m : -m;
				}
			}
			if (temp < ans || flag)
			{
				ans = temp;
				mx = i, my = j;
				flag = false;
			}
		}
	}
	*x = mx;
	*y = my;
	return 1;
}

int ustc_SubImgMatch_hist(Mat grayImg, Mat subImg, int* x, int* y)
{
	int ghist[1000], shist[1000];
	int mx = 0, my = 0;
	int i = 0, j = 0, k = 0, l = 0, m = 0;
	int row = grayImg.rows, col = grayImg.cols;
	int rrr = subImg.rows, ccc = subImg.cols;
	int temp, ans = 0;
	bool flag = true;
	uchar *pg(NULL), *ps(NULL);
	row -= rrr;
	col -= ccc;
	grayImg.convertTo(grayImg, CV_8UC1);
	subImg.convertTo(subImg, CV_8UC1);
	for (i = 0;i < row;i++)
	{
		for (j = 0;j < col;j++)
		{
			temp = 0;
			memset(ghist, 0, 256*sizeof(int));
			memset(shist, 0, 256*sizeof(int));
			for (k = 0;k < rrr;k++)
			{
				pg = grayImg.ptr<uchar>(i + k);
				ps = subImg.ptr<uchar>(k);
				for (l = 0;l < ccc;l++)
				{
					ghist[pg[j + l]]++;
					shist[ps[l]]++;
				}
			}
			for (k = 0;k < 256;k++)
			{
//				cout << "  " << ghist[k];
				m = ghist[k] - shist[k];
				temp += m > 0 ? m : -m;
			}
			if (temp < ans|| flag)
			{
				ans = temp;
				mx = i, my = j;
				flag = false;
			}
		}
	}
	*x = mx;
	*y = my;
	return 1;
}
