#include "SubImageMatch.h"
#include<opencv2/opencv.hpp>
#include<cv.h>
using namespace cv;
using namespace std;
#include <iostream>
#include <math.h>
#include<stdio.h>
#define SUB_IMAGE_MATCH_OK 1
#define SUB_IMAGE_MATCH_FAIL -1

int ustc_ConvertBgr2Gray(Mat bgrImg, Mat& grayImg)
{
	if (NULL == bgrImg.data)
	{
		cout << "image is NULL." << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}

	int width = bgrImg.cols;
	int height = bgrImg.rows;
	for (int row_i = 0; row_i<height; row_i += 1)
	{
		for (int col_j = 0; col_j<width; col_j += 1)
		{
			int location = 3 * (row_i * width + col_j);
			int b = bgrImg.data[location + 0];
			int g = bgrImg.data[location + 1];
			int r = bgrImg.data[location + 2];
			

			int grayVal = (b * 19595 + g * 38469 + r * 7472) >> 16;
		
			grayImg.data[location / 3] = grayVal;
			
		}
	}
	return 0;
}

int ustc_CalcGrad(Mat grayImg, Mat& gradImg_x, Mat& gradImg_y)

{

	if (NULL == grayImg.data)

	{

		cout << "image is NULL." << endl;

		return SUB_IMAGE_MATCH_FAIL;

	}

	int width = grayImg.cols;

	int height = grayImg.rows;

	int temp0, temp1;

	int row_i, col_j;

	int grad_x, grad_y;

	for (row_i = 0; row_i < height; row_i++)

	{

		temp0 = row_i * width;

		gradImg_x.data[temp0] = 0;

		gradImg_x.data[temp0 + width - 1] = 0;

		gradImg_y.data[temp0] = 0;

		gradImg_y.data[temp0 + width - 1] = 0;

	}

	for (col_j = 1; col_j < width; col_j++)

	{

		gradImg_x.data[col_j] = 0;

		gradImg_x.data[temp0 + col_j] = 0;

		gradImg_y.data[col_j] = 0;

		gradImg_y.data[temp0 + col_j] = 0;

	}

	for (row_i = 1; row_i < height - 1; row_i++)

	{

		temp0 = row_i * width;

		for (col_j = 1; col_j < width - 1; col_j++)

		{

			temp1 = temp0 + col_j;

			grad_x =

				grayImg.data[temp1 - width + 1]

				+ 2 * grayImg.data[temp1 + 1]

				+ grayImg.data[temp1 + width + 1]

				- grayImg.data[temp1 - width - 1]

				- 2 * grayImg.data[temp1 - 1]

				- grayImg.data[temp1 + width - 1];

			((float*)gradImg_x.data)[temp1] = grad_x;



			grad_y =

				grayImg.data[temp1 + width + 1]

				+ 2 * grayImg.data[temp1 + width]

				+ grayImg.data[temp1 + width - 1]

				- grayImg.data[temp1 - width - 1]

				- 2 * grayImg.data[temp1 - width]

				- grayImg.data[temp1 - width + 1];

			((float*)gradImg_y.data)[temp1] = grad_y;

		}

	}



	return SUB_IMAGE_MATCH_OK;

}

int ustc_CalcAngleMag(Mat gradImg_x, Mat gradImg_y, Mat& angleImg, Mat& magImg)
{
	if (NULL == gradImg_x.data || NULL == gradImg_y.data)
	{
		cout << "image is NULL." << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}

	int width = gradImg_x.cols;
	int height = gradImg_x.rows;




	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j += 1)
		{
			float grad_x = ((float*)gradImg_x.data)[row_i * width + col_j];
			float grad_y = ((float*)gradImg_y.data)[row_i * width + col_j];
			float angle = atan2(grad_y, grad_x);


			((float*)angleImg.data)[row_i * width + col_j] = angle;
		}
	}

	for (int row_i = 0; row_i < height; row_i++)
	{
		for (int col_j = 0; col_j < width; col_j += 1)
		{
			float a = ((float*)gradImg_x.data)[row_i * width + col_j];
			float b = ((float*)gradImg_y.data)[row_i * width + col_j];
			float val = a*a + b*b;
			int t = *(int*)&val;
			t -= 0x3f800000; t >>= 1;
			t += 0x3f800000;

			magImg.data[row_i * width + col_j] = t;
		}
	}
	return 0;
}

int ustc_Threshold(Mat grayImg, Mat& binaryImg, int th)

{

	if (NULL == grayImg.data)

	{

		cout << "image is NULL." << endl;

		return SUB_IMAGE_MATCH_FAIL;

	}

	uchar transition[256];

	for (int i = 0; i <= th; i++)

	{

		transition[i] = 0;

	}

	for (int i = th; i <256; i++)

	{

		transition[i] = 255;

	}

	int width = grayImg.cols;

	int height = grayImg.rows;

	for (int row_i = 0; row_i < height; row_i++)

	{

		int temp0 = row_i * width;

		for (int col_j = 0; col_j < width; col_j++)

		{

			int temp1 = temp0 + col_j;

			int pixVal = grayImg.data[temp1];

			pixVal = transition[pixVal];

			binaryImg.data[temp1] = pixVal;

		}

	}

	return SUB_IMAGE_MATCH_OK;

}

int ustc_CalcHist(Mat grayImg, int* hist, int hist_len)

{

	if (NULL == grayImg.data)

	{

		cout << "image is NULL." << endl;

		return SUB_IMAGE_MATCH_FAIL;

	}

	int width = grayImg.cols;

	int height = grayImg.rows;

	uchar*p = grayImg.data;

	for (int i = 0; i < hist_len; i++)

	{

		hist[i] = 0;

	}

	for (int row_i = height; row_i >0; row_i--)

	{

		for (int col_j = width; col_j >0; col_j--)

		{

			int pixVal = *p;

			hist[pixVal]++;

			p++;

		}

	}

	return SUB_IMAGE_MATCH_OK;

}



int ustc_SubImgMatch_gray(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (grayImg.channels() != subImg.channels())
	{
		cout << "images have different channels" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;

	if (width < sub_width || height < sub_height)
	{
		cout << "the subimage is larger ,error input" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}



	Mat searchImg(height, width, CV_32FC1);

	searchImg.setTo(FLT_MAX);


	for (int i = 0; i <= height - sub_height; i++)
	{
		for (int j = 0; j <= width - sub_width; j++)
		{
			int total_diff = 0;

			for (int y = 0; y < sub_height; y++)
			{
				for (int x = 0; x < sub_width; x++)
				{

					int row_index = i + y;
					int col_index = j + x;
					int bigImg_pix = grayImg.data[row_index * width + col_index];

					int template_pix = subImg.data[y * sub_width + x];

					total_diff += abs(bigImg_pix - template_pix);
				}
			}

			((float*)searchImg.data)[i * width + j] = total_diff;
		}
	}


	int min = searchImg.data[0];
	*x = 0, *y = 0;
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			if (((float*)searchImg.data)[i * width + j] <= min)
			{
				min = ((float*)searchImg.data)[i * width + j];

				*y = i;
				*x = j;
			}
		}
	}



	return 0;
}
int ustc_SubImgMatch_bgr(Mat colorImg, Mat subImg, int* x, int* y)
{
	if (NULL == colorImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (colorImg.channels() != subImg.channels())
	{
		cout << "images have different channels" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}

	int width = colorImg.cols;
	int height = colorImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	if (width < sub_width || height < sub_height)
	{
		cout << "the subimage is larger ,error input" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}



	Mat searchImg(height, width, CV_32FC1);

	searchImg.setTo(FLT_MAX);


	for (int i = 0; i <= height - sub_height; i++)
	{
		for (int j = 0; j <= width - sub_width; j++)
		{
			int total_diffb = 0;
			int total_diffg = 0;
			int total_diffr = 0;
			int total_diff = 0;

			for (int y = 0; y < sub_height; y++)
			{
				for (int x = 0; x < sub_width; x++)
				{

					int row_index = i + y;
					int col_index = j + x;
					int  bigImg_b = colorImg.data[3 * (row_index * width + col_index) + 0];
					int  bigImg_g = colorImg.data[3 * (row_index * width + col_index) + 1];
					int  bigImg_r = colorImg.data[3 * (row_index * width + col_index) + 2];

					int template_b = subImg.data[3 * (y * sub_width + x) + 0];
					int template_g = subImg.data[3 * (y * sub_width + x) + 1];
					int template_r = subImg.data[3 * (y * sub_width + x) + 2];

					total_diffb += abs(bigImg_b - template_b);
					total_diffg += abs(bigImg_b - template_b);
					total_diffr += abs(bigImg_b - template_b);
					total_diff += (total_diffr + total_diffg + total_diffb);

				}
			}

			((float*)searchImg.data)[(i * width + j)] = total_diff;
		}
	}

	float min = ((float*)searchImg.data)[0];
	*x = 0, *y = 0;
	
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			if (((float*)searchImg.data)[(i * width + j)] < min)
			{
				min = ((float*)searchImg.data)[i * width + j];

				*y = i;
				*x = j;
			}
		}
	}




	return 0;
}

int ustc_SubImgMatch_corr(Mat grayImg, Mat subImg, int* x, int* y)

{

	if (NULL == grayImg.data || NULL == subImg.data)

	{

		cout << "image is NULL." << endl;

		return SUB_IMAGE_MATCH_FAIL;

	}

	int width = grayImg.cols;

	int height = grayImg.rows;

	int sub_width = subImg.cols;

	int sub_height = subImg.rows;

	int sub_i, sub_j, i, j;

	int resu1 = 0, resu2 = 0;

	float value = 0;

	float st1 = 0, st2 = 0, st3 = 0;

	for (sub_i = 0; sub_i< sub_height; sub_i++)

	{

		for (sub_j = 0; sub_j< sub_width; sub_j++)

		{

			int template_pix = subImg.data[sub_i * sub_width + sub_j];

			st3 += template_pix*template_pix;

		}

	}

	st3 = 1 / st3;

	for (i = 0; i <= (height - sub_height); i++)

	{

		for (j = 0; j <= (width - sub_width); j++)

		{

			st1 = 0;

			st2 = 0;

			for (sub_i = 0; sub_i< sub_height; sub_i++)

			{

				for (sub_j = 0; sub_j< sub_width; sub_j++)

				{

					int row_index = i + sub_i;

					int col_index = j + sub_j;

					int bigImg_pix = grayImg.data[row_index * width + col_index];

					int template_pix = subImg.data[sub_i * sub_width + sub_j];

					st1 += bigImg_pix*template_pix;

					st2 += bigImg_pix*bigImg_pix;

				}

			}

			st1 = st1*st1;

			st2 = 1 / st2;

			float rej = st1*st2*st3;

			if (rej > value)

			{

				value = rej;

				resu1 = i;

				resu2 = j;

			}

		}

	}

	*x = resu2;

	*y = resu1;

	return SUB_IMAGE_MATCH_OK;

}
int ustc_SubImgMatch_angle(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (grayImg.channels() != subImg.channels())
	{
		cout << "images have different channels" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	if (width < sub_width || height < sub_height)
	{
		cout << "the subimage is larger ,error input" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}



	Mat searchImg(height, width, CV_32FC1);

	searchImg.setTo(FLT_MAX);


	for (int i = 0; i <= height - sub_height - 1; i++)
	{
		for (int j = 0; j <= width - sub_width - 1; j++)
		{
			int total_diff = 0;

			for (int y = 1; y < sub_height - 1; y++)
			{
				for (int x = 1; x < sub_width - 1; x++)
				{

					int row_index = i + y;
					int col_index = j + x;


					int bigImg_grad_x =
						grayImg.data[(row_index - 1) * width + col_index + 1]
						+ 2 * grayImg.data[(row_index)* width + col_index + 1]
						+ grayImg.data[(row_index + 1)* width + col_index + 1]
						- grayImg.data[(row_index - 1) * width + col_index - 1]
						- 2 * grayImg.data[(row_index)* width + col_index - 1]
						- grayImg.data[(row_index + 1)* width + col_index - 1];

					int bigImg_grad_y =
						-grayImg.data[(row_index - 1) * width + col_index - 1]
						- 2 * grayImg.data[(row_index - 1)* width + col_index]
						- grayImg.data[(row_index - 1)* width + col_index + 1]
						+ grayImg.data[(row_index + 1) * width + col_index - 1]
						+ 2 * grayImg.data[(row_index + 1)* width + col_index]
						+ grayImg.data[(row_index + 1)* width + col_index + 1];

					float bigImg_angle = atan2(bigImg_grad_y, bigImg_grad_x);


					int template_grad_x =
						subImg.data[(y - 1) * width + x + 1]
						+ 2 * subImg.data[(y)* width + x + 1]
						+ subImg.data[(y + 1)* width + x + 1]
						- subImg.data[(y - 1) * width + x - 1]
						- 2 * subImg.data[(y)* width + x - 1]
						- subImg.data[(y + 1)* width + x - 1];

					int template_grad_y =
						-subImg.data[(y - 1) * sub_width + x - 1]
						- 2 * subImg.data[(y - 1)* sub_width + x]
						- subImg.data[(y - 1)* sub_width + x + 1]
						+ subImg.data[(y + 1) * sub_width + x - 1]
						+ 2 * subImg.data[(y + 1)* sub_width + x]
						+ subImg.data[(y + 1)* sub_width + x + 1];

					float template_angle = atan2(template_grad_y, template_grad_x);
					total_diff += abs(bigImg_angle - template_angle);
				}
			}
			((float*)searchImg.data)[i * width + j] = total_diff;
		}
	}

	float min = ((float*)searchImg.data)[0];
	*x = 0, *y = 0;
	for (int i = 0; i <= height - sub_height; i++)
	{
		for (int j = 0; j <= width - sub_width; j++)
		{
			if (((float*)searchImg.data)[i * width + j]<min)
			{
				min = ((float*)searchImg.data)[i * width + j];
				*y = i;
				*x = j;
			}
		}
	}


	return 0;

}
int ustc_SubImgMatch_mag(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (grayImg.channels() != subImg.channels())
	{
		cout << "images have different channels" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	if (width < sub_width || height < sub_height)
	{
		cout << "the subimage is larger ,error input" << endl;
		return SUB_IMAGE_MATCH_FAIL;
	}



	Mat searchImg(height, width, CV_32FC1);

	searchImg.setTo(FLT_MAX);


	for (int i = 0; i <= height - sub_height - 1; i++)
	{
		for (int j = 0; j <= width - sub_width - 1; j++)
		{
			int total_diff = 0;

			for (int y = 1; y < sub_height - 1; y++)
			{
				for (int x = 1; x < sub_width - 1; x++)
				{

					int row_index = i + y;
					int col_index = j + x;


					int bigImg_grad_x =
						grayImg.data[(row_index - 1) * width + col_index + 1]
						+ 2 * grayImg.data[(row_index)* width + col_index + 1]
						+ grayImg.data[(row_index + 1)* width + col_index + 1]
						- grayImg.data[(row_index - 1) * width + col_index - 1]
						- 2 * grayImg.data[(row_index)* width + col_index - 1]
						- grayImg.data[(row_index + 1)* width + col_index - 1];

					int bigImg_grad_y =
						-grayImg.data[(row_index - 1) * width + col_index - 1]
						- 2 * grayImg.data[(row_index - 1)* width + col_index]
						- grayImg.data[(row_index - 1)* width + col_index + 1]
						+ grayImg.data[(row_index + 1) * width + col_index - 1]
						+ 2 * grayImg.data[(row_index + 1)* width + col_index]
						+ grayImg.data[(row_index + 1)* width + col_index + 1];


					float bigImg_mag = sqrt(bigImg_grad_y*bigImg_grad_y + bigImg_grad_x*bigImg_grad_x);


					int template_grad_x =
						subImg.data[(y - 1) * width + x + 1]
						+ 2 * subImg.data[(y)* width + x + 1]
						+ subImg.data[(y + 1)* width + x + 1]
						- subImg.data[(y - 1) * width + x - 1]
						- 2 * subImg.data[(y)* width + x - 1]
						- subImg.data[(y + 1)* width + x - 1];

					int template_grad_y =
						-subImg.data[(y - 1) * sub_width + x - 1]
						- 2 * subImg.data[(y - 1)* sub_width + x]
						- subImg.data[(y - 1)* sub_width + x + 1]
						+ subImg.data[(y + 1) * sub_width + x - 1]
						+ 2 * subImg.data[(y + 1)* sub_width + x]
						+ subImg.data[(y + 1)* sub_width + x + 1];

					float template_mag = sqrt(template_grad_y*template_grad_y + template_grad_x*template_grad_x);
					total_diff += abs(bigImg_mag - template_mag);
				}
			}
			((float*)searchImg.data)[i * width + j] = total_diff;
		}
	}

	float min = ((float*)searchImg.data)[0];
	*x = 0, *y = 0;
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			if (((float*)searchImg.data)[i * width + j]<min)
			{
				min = ((float*)searchImg.data)[i * width + j];
				
				*y = i;
				*x = j;
			}
		}
	}

	return 0;

}
int ustc_SubImgMatch_hist(Mat grayImg, Mat subImg, int* x, int* y)

{

	if (NULL == grayImg.data || NULL == subImg.data)

	{

		cout << "image is NULL." << endl;

		return SUB_IMAGE_MATCH_FAIL;

	}

	int sign[2] = { 1,-1 };

	int width = grayImg.cols;

	int height = grayImg.rows;

	int sub_width = subImg.cols;

	int sub_height = subImg.rows;

	int subhist[256] = { 0 };

	ustc_CalcHist(subImg, subhist, 256);

	int bigImg_pix = 0, total_diff = 0;

	int hist[256] = { 0 };

	int hist1[256] = { 0 };

	int sub_i, sub_j, i, j;

	int var_i = 0, var_j = 0;

	int flag_0;

	for (sub_i = 0; sub_i < sub_height; sub_i++)

	{

		for (sub_j = 0; sub_j < sub_width; sub_j++)

		{

			bigImg_pix = grayImg.data[sub_i * width + sub_j];

			hist[bigImg_pix]++;

		}

	}

	for (int ii = 0; ii < 256; ii++)

	{

		hist1[ii] = hist[ii];

		int temp = hist[ii] - subhist[ii];

		if (temp < 0) temp = -temp;

		total_diff += temp;

	}

	flag_0 = total_diff;

	for (i = 0; i <= (height - sub_height); i++)

	{

		if (i != 0)

		{

			total_diff = 0;

			for (sub_i = 0; sub_i < sub_width; sub_i++)

			{

				bigImg_pix = grayImg.data[(i - 1) * width + sub_i];

				hist[bigImg_pix]--;

				bigImg_pix = grayImg.data[(i + sub_height - 1) * width + sub_i];

				hist[bigImg_pix]++;

			}

			for (int ii = 0; ii < 256; ii++)

			{

				hist1[ii] = hist[ii];

				int temp = hist[ii] - subhist[ii];

				if (temp < 0) temp = -temp;

				total_diff += temp;

			}

			if (flag_0 > total_diff)

			{

				flag_0 = total_diff;

				var_i = i;

				var_j = j;

			}

		}

		for (j = 1; j <= (width - sub_width); j++)

		{

			total_diff = 0;

			for (sub_j = 0; sub_j < sub_height; sub_j++)

			{

				bigImg_pix = grayImg.data[(i + sub_j) * width + j - 1];

				hist1[bigImg_pix]--;

				bigImg_pix = grayImg.data[(i + sub_j) * width + j + sub_width - 1];

				hist1[bigImg_pix]++;

			}

			for (int ii = 0; ii < 256; ii++)

			{

				int temp = hist1[ii] - subhist[ii];

				if (temp < 0) temp = -temp;

				total_diff += temp;

			}

			if (flag_0 > total_diff)

			{

				flag_0 = total_diff;

				var_i = i;

				var_j = j;

			}

		}

	}

	*x = var_j;

	*y = var_i;

	return SUB_IMAGE_MATCH_OK;

}
