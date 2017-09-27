#include "SubImageMatch.h"

// 作业1------------------------------------------------------------------
int ustc_ConvertBgr2Gray(Mat bgrImg, Mat& grayImg)
{
	if (NULL == bgrImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}

	int width = bgrImg.cols;
	int height = bgrImg.rows;

	for (int row_i = 0; row_i < height; row_i++)
	{
		for (int col_j = 0; col_j < width; col_j++)
		{
			int i = row_i * width + col_j;
			int b = bgrImg.data[3 * i + 0];
			int g = bgrImg.data[3 * i + 1];
			int r = bgrImg.data[3 * i + 2];
			int grayVal = (b * 114 + g * 587 + r * 229) / 1000;
			grayImg.data[i] = grayVal;
		}
	}

}



//-----------------------------------------------------------------------------------

//作业2

int ustc_CalcGrad(Mat grayImg, Mat& gradImg_x, Mat& gradImg_y)
{
	if (NULL == grayImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j += 1)
		{
			int grad_x =
				grayImg.data[(row_i - 1) * width + col_j + 1]
				+ 2 * grayImg.data[(row_i)* width + col_j + 1]
				+ grayImg.data[(row_i + 1)* width + col_j + 1]
				- grayImg.data[(row_i - 1) * width + col_j - 1]
				- 2 * grayImg.data[(row_i)* width + col_j - 1]
				- grayImg.data[(row_i + 1)* width + col_j - 1];
			((float*)gradImg_x.data)[row_i * width + col_j] = grad_x;
			int grad_y =
				grayImg.data[(row_i + 1) * width + col_j - 1]
				+ 2 * grayImg.data[(row_i + 1)* width + col_j]
				+ grayImg.data[(row_i + 1)* width + col_j + 1]
				- grayImg.data[(row_i - 1) * width + col_j - 1]
				- 2 * grayImg.data[(row_i - 1)* width + col_j]
				- grayImg.data[(row_i - 1)* width + col_j + 1];
			((float*)gradImg_y.data)[row_i * width + col_j] = grad_y;
		}
	}



}


//------------------------------------------------------------------

//3.

int ustc_CalcAngleMag(Mat gradImg_x, Mat gradImg_y, Mat& angleImg, Mat& magImg)
{
	if (NULL == gradImg_x.data || NULL == gradImg_y.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}

	int width = gradImg_x.cols;
	int height = gradImg_x.rows;
	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j++)
		{
			int grad_x = ((float*)gradImg_x.data)[row_i * width + col_j];
			int grad_y = ((float*)gradImg_y.data)[row_i * width + col_j];
			int angle = atan2(grad_y, grad_x);
	//		int x = grad_y * grad_y + grad_x * grad_x;
	//		int mag =(( x- 0x3f800000)>>1)+ 0x3f800000;
			int mag = sqrt(grad_y * grad_y + grad_x * grad_x);
			((float*)angleImg.data)[row_i * width + col_j] = angle;
			((float*)magImg.data)[row_i * width + col_j] = mag;
		}
	}

}


//--------------------------------------------------------------
//4.
int ustc_Threshold(Mat grayImg,Mat& binaryImg,int th)
{
	if (NULL == grayImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;

	for (int row_i = 0; row_i < height; row_i++)
	{
		int temp0 = row_i * width;
		for (int col_j = 0; col_j < width; col_j += 2)
		{
			//int pixVal = grayImg.at<uchar>(row_i, col_j);
			int temp1 = temp0 + col_j;
			int pixVal = grayImg.data[temp1];
			int dstVal = 0;
			if (pixVal > th)
			{
				dstVal = 255;
			}
			else 
			{
				dstVal = 0;
			}
			//binaryImg.at<uchar>(row_i, col_j) = dstVal;
			(binaryImg).data[temp1] = dstVal;
		}
	}

#ifdef IMG_SHOW
	namedWindow("binaryImg", 0);
	imshow("binaryImg", (binaryImg));
	waitKey();
#endif

	return MY_OK;
}




//-------------------------------------------------------
//5.
int ustc_CalcHist(Mat grayImg, int* hist,int hist_len)
{
	if (NULL == grayImg.data || NULL == hist)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;
	
	//计算直方图
	for (int row_i = 0; row_i < height; row_i++)
	{
		int temp0 = row_i * width;
		for (int col_j = 0; col_j < width; col_j += 1)
		{
			int pixVal = grayImg.data[temp0 + col_j];
			hist[pixVal]++;

		}
	}
	int max = 0;                   //确定显示直方图的高度
	for (int i = 0; i < hist_len; ++i)
	{
		if (max < hist[i])
		{
			max = hist[i];
		}
	}
	
	
	
#ifdef IMG_SHOW
	Mat showImg(max+100,256, CV_8UC1);
	showImg.setTo(0);
	for (int i = 0; i < hist_len; ++i)
	{
		
		line(showImg,Point(i, max + 100),Point(i, max + 100 -hist[i]), 255,1, 8, 0);
	}
	namedWindow("showImg", 0);
	imshow("showImg", showImg);
	waitKey();
#endif
}


//---------------------------------------------------------------------
//6.
int ustc_SubImgMatch_gray(Mat grayImg, Mat subImg, int *x, int *y)
{

	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;

	//该图用于记录每一个像素位置的匹配误差
	Mat searchImg(height, width, CV_32FC1);
	searchImg.setTo(FLT_MAX);
	//遍历大图每一个像素，注意行列的起始、终止坐标
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			int total_diff = 0;
			//遍历模板图上的每一个像素
			for (int x = 0; x < sub_height; x++)
			{
				for (int y = 0; y < sub_width; y++)
				{
					//大图上的像素位置
					int row_index = i + x;
					int col_index = j + y;
					int bigImg_pix = grayImg.data[row_index * width + col_index];
					//模板图上的像素
					int template_pix = subImg.data[x * sub_width + y];
					total_diff += abs(bigImg_pix - template_pix);
				}
			}
			//存储当前像素位置的匹配误差

			((float*)searchImg.data)[i * width + j] = total_diff;
		}
	}
	float total_diff = FLT_MAX;
	int min_width = 0;
	int min_height = 0;
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			if (((float*)searchImg.data)[i * width + j] < total_diff)
			{
				total_diff = searchImg.data[i * width + j];
				min_width = j;
				min_height = i;
				
			}
		}
	}

	*y = min_height;
	*x = min_width;
	


}

//----------------------------------------------------------
//7.
int ustc_SubImgMatch_bgr(Mat grayImg, Mat subImg,int* x,int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	
	Mat searchImg(height, width, CV_32FC1);
	searchImg.setTo(FLT_MAX);
	//遍历大图每一个像素，注意行列的起始、终止坐标
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			int total_diff = 0;
			//遍历模板图上的每一个像素
			for (int x = 0; x < sub_height; x++)
			{
				for (int y = 0; y < sub_width; y++)
				{
					//大图上的像素位置
					int row_index = i + x;
					int col_index = j + y;
					int bigImg_pix = grayImg.data[3 * (row_index * width + col_index) + 0] +
						grayImg.data[3 * (row_index * width + col_index) + 1] +
						grayImg.data[3 * (row_index * width + col_index) + 2];
					//模板图上的像素
					int template_pix = subImg.data[3 * (x * sub_width + y) + 0] +
						subImg.data[3 * (x * sub_width + y) + 1] +
						subImg.data[3 * (x * sub_width + y) + 2];
					total_diff += abs(bigImg_pix - template_pix);
				}
			}
			//存储当前像素位置的匹配误差
			((float*)searchImg.data)[i * width + j] = total_diff;
		}
	}
	float total_diff = ((float*)searchImg.data)[0];
	int min_width = 0;
	int min_height = 0;
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			if (((float*)searchImg.data)[i * width + j] < total_diff)
			{
				total_diff = ((float*)searchImg.data)[i * width + j];
				min_width = i;
				min_height = j;
			}
		}
	}
	*x = min_height;
	*y= min_width ;
}


//-------------------------------------------------------
//8.
int ustc_SubImgMatch_corr(Mat grayImg, Mat subImg,int* x,int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;

	Mat searchImg(height, width, CV_32FC1);
	searchImg.setTo(0);
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			float total_st = 0;
			float total_s = 0;
			float total_t = 0;
			for (int x = 0; x < sub_height; x++)
			{
				for (int y = 0; y < sub_width; y++)
				{

					int row_index = i + x;
					int col_index = j + y;
					int bigImg_pix = grayImg.data[row_index * width + col_index];
					int template_pix = subImg.data[x * sub_width + y];
					total_st += bigImg_pix * template_pix;
					total_s += bigImg_pix * bigImg_pix;
					total_t += template_pix * template_pix;
				}
			}

			((float*)searchImg.data)[i * width + j] = total_st / (sqrt(total_s) * sqrt(total_t));

		}
	}
	float total_diff = 0;
	int min_width = 0;
	int min_height = 0;
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			float test = ((float*)searchImg.data)[i * width + j];
			if (abs(1 - total_diff) > abs(1 - test))
			{
				total_diff = test;
				min_width = i;
				min_height = j;
			}
		}
	}
	*x = min_height; 
	*y = min_width;
}

//-------------------------------------------------------
//9.
int ustc_SubImgMatch_angle(Mat grayImg, Mat subImg,int* x,int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;
	Mat gradImg_x(height, width, CV_32FC1);
	Mat gradImg_y(height, width, CV_32FC1);
	gradImg_x.setTo(0);
	gradImg_y.setTo(0);
	int flag = ustc_CalcGrad(grayImg, gradImg_x, gradImg_y);
	Mat angleImg(height, width, CV_32FC1);
	angleImg.setTo(0);
	//计算角度图
	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j ++)
		{
			int grad_x = ((float*)gradImg_x.data)[row_i * width + col_j];
			int grad_y = ((float*)gradImg_y.data)[row_i * width + col_j];
			float angle = atan2(grad_y, grad_x);
			((float*)angleImg.data)[row_i * width + col_j] = angle * 180 / 3.14;
			if (((float*)angleImg.data)[row_i * width + col_j] < 0)
				((float*)angleImg.data)[row_i * width + col_j] = ((float*)angleImg.data)[row_i * width + col_j] + 360;
		}
	}

	int width_sub = subImg.cols;
	int height_sub = subImg.rows;
	Mat subImg_x(height_sub, width_sub, CV_32FC1);
	Mat subImg_y(height_sub, width_sub, CV_32FC1);
	subImg_x.setTo(0);
	subImg_y.setTo(0);
	int flag2 = ustc_CalcGrad(subImg, subImg_x, subImg_y);
	Mat angleImg_sub(height_sub, width_sub, CV_32FC1);
	angleImg_sub.setTo(0);
	for (int row_i = 1; row_i < height_sub - 1; row_i++)
	{
		for (int col_j = 1; col_j < width_sub - 1; col_j += 1)
		{
			int grad_x = ((float*)subImg_x.data)[row_i * width_sub + col_j];
			int grad_y = ((float*)subImg_y.data)[row_i * width_sub + col_j];
			float angle = atan2(grad_y, grad_x);
			((float *)angleImg_sub.data)[row_i * width_sub + col_j] = angle * 180 / 3.14;
			if (((float*)angleImg.data)[row_i * width + col_j] < 0)
				((float*)angleImg.data)[row_i * width + col_j] = ((float*)angleImg.data)[row_i * width + col_j] + 360;
		}
	}

	Mat searchImg(height, width, CV_32FC1);
	searchImg.setTo(0);
	for (int i = 0; i < height - height_sub; i++)
	{
		for (int j = 0; j < width - width_sub; j++)
		{
			int total_diff = 0;
			for (int x = 0; x < height_sub; x++)
			{
				for (int y = 0; y < width_sub; y++)
				{
					int row_index = i + x;
					int col_index = j + y;
					int bigImg_pix = ((float*)angleImg.data)[row_index * width + col_index];
					int template_pix = ((float*)angleImg_sub.data)[x * width_sub + y];
					int answer = bigImg_pix - template_pix;
					if (answer > 180)
					{
						answer = 360 - answer;
					}
					else if (answer < -180)
					{
						answer = 360 + answer;
					}
					total_diff += answer;
				}
			}

			((float*)searchImg.data)[i * width + j] = total_diff;
			waitKey(0);
		}
	}

	
	float total_diff = ((float*)searchImg.data)[0];
	int min_width = 0;
	int min_height = 0;
	for (int i = 0; i < height - height_sub; i++)
	{
		for (int j = 0; j < width - width_sub; j++)
		{
			if (((float*)searchImg.data)[i * width + j] < total_diff)
			{
				total_diff = ((float*)searchImg.data)[i * width + j];
				min_width = i;
				min_height = j;
			}
		}
	}
	*x = min_height;
	*y = min_width;

}


//------------------------------------------------------
//10.
int ustc_SubImgMatch_mag(Mat grayImg, Mat subImg,int* x,int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}

	int width = grayImg.cols;
	int height = grayImg.rows;

	Mat gradImg_x(height, width, CV_32FC1);
	Mat gradImg_y(height, width, CV_32FC1);
	gradImg_x.setTo(0);
	gradImg_y.setTo(0);
	int flag = ustc_CalcGrad(grayImg, gradImg_x, gradImg_y);
	Mat magImg(height, width, CV_32SC1);
	magImg.setTo(0);

	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j += 1)
		{
			int grad_x = ((float*)gradImg_x.data)[row_i * width + col_j];
			int grad_y = ((float*)gradImg_y.data)[row_i * width + col_j];
			int mag = sqrt(grad_y * grad_y + grad_x * grad_x);
			((int*)magImg.data)[row_i * width + col_j] = mag;
		}
	}

	int width_sub = subImg.cols;
	int height_sub = subImg.rows;
	Mat subImg_x(height_sub, width_sub, CV_32FC1);
	Mat subImg_y(height_sub, width_sub, CV_32FC1);
	subImg_x.setTo(0);
	subImg_y.setTo(0);

	int flag2 = ustc_CalcGrad(subImg, subImg_x, subImg_y);
	Mat magImg_sub(height_sub, width_sub, CV_32SC1);
	magImg_sub.setTo(0);
	//计算角度图
	for (int row_i = 1; row_i < height_sub - 1; row_i++)
	{
		for (int col_j = 1; col_j < width_sub - 1; col_j += 1)
		{
			int grad_x = ((float*)subImg_x.data)[row_i * width_sub + col_j];
			int grad_y = ((float*)subImg_y.data)[row_i * width_sub + col_j];
			int mag = sqrt(grad_y * grad_y + grad_x * grad_x);
			//自己找办法优化三角函数速度，并且转化为角度制，规范化到0-360
			((int*)magImg_sub.data)[row_i * width_sub + col_j] = mag;
		}
	}

	Mat searchImg(height, width, CV_32FC1);
	//匹配误差初始化
	searchImg.setTo(0);
	//遍历大图每一个像素，注意行列的起始、终止坐标
	for (int i = 0; i < height - height_sub; i++)
	{
		for (int j = 0; j < width - width_sub; j++)
		{
		
			int total_diff = 0;
			for (int x = 0; x < height_sub; x++)
			{
				for (int y = 0; y < width_sub; y++)
				{
					int row_index = i + x;
					int col_index = j + y;
					int bigImg_pix = ((int*)magImg.data)[row_index * width + col_index];
					int template_pix = ((int*)magImg_sub.data)[x * width_sub + y];
					int answer = abs(bigImg_pix - template_pix);
					total_diff += answer;
				}
			}
			//存储当前像素位置的匹配误差

			((float*)searchImg.data)[i * width + j] = total_diff;


		}
	}
	float total_diff = ((float*)searchImg.data)[0];
	int min_width = 0;
	int min_height = 0;
	for (int i = 0; i < height - height_sub; i++)
	{
		for (int j = 0; j < width - width_sub; j++)
		{
			if (((float*)searchImg.data)[i * width + j] < total_diff)
			{
				total_diff = ((float*)searchImg.data)[i * width + j];
				min_width = i;
				min_height = j;
			}
		}
	}
	*x = min_height;
	*y = min_width;
}


//-------------------------------------------------------------
//11.
int ustc_Cal(Mat subImg, int* hist)
{
	if (NULL == subImg.data || NULL == hist)
	{
		cout << "image is NULL." << endl;
		return NULL;
	}

	int width = subImg.cols;
	int height = subImg.rows;

	//计算直方图
	for (int row_i = 0; row_i < height; row_i++)
	{
		for (int col_j = 0; col_j < width; col_j += 1)
		{
			int pixVal = subImg.data[row_i * width + col_j];
			hist[pixVal]++;
		}
	}

	waitKey();


}
int ustc_SubImgMatch_hist(Mat grayImg, Mat subImg, int* x, int* y)
{
	
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return MY_FAIL;
	}
	int hist[256];
	for (int i = 0; i < 256; ++i)
	{
		hist[i] = 0;
	}
	int flag3 = ustc_Cal(subImg, hist);

	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	Mat searchImg(height, width, CV_32FC1);
	searchImg.setTo(FLT_MAX);
	int* hist_temp = new int[256];
	memset(hist_temp, 0, sizeof(int) * 256);
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{
			//清零
			memset(hist_temp, 0, sizeof(int) * 256);

			//计算当前位置直方图
			for (int x = 0; x < sub_height; x++)
			{
				for (int y = 0; y < sub_width; y++)
				{
					//大图上的像素位置
					int row_index = i + x;
					int col_index = j + y;
					int bigImg_pix = grayImg.data[row_index * width + col_index];
					hist_temp[bigImg_pix]++;
				}
			}

			//根据直方图计算匹配误差
			int total_diff = 0;
			for (int ii = 0; ii < 256; ii++)
			{
				total_diff += abs(hist_temp[ii] - hist[ii]);
			}
			//存储当前像素位置的匹配误差
			((float*)searchImg.data)[i * width + j] = total_diff;
			//cout << ((float*)searchImg.data)[i * width + j] << endl;
		}
	}
	int min = 10000000;
	int min_height = 0;
	int min_width = 0;
	delete[] hist_temp;
	for (int i = 0; i < height - sub_height; i++)
	{
		for (int j = 0; j < width - sub_width; j++)
		{

			if (min >((float*)searchImg.data)[i * width + j])
			{
				min = ((float*)searchImg.data)[i * width + j];
				min_height = j;
				min_width = i;
			}
		}
	}

	*x = min_height;
	*y = min_width;



}


