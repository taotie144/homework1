#include "SubImageMatch.h"


int ustc_ConvertBgr2Gray(Mat bgrImg, Mat& grayImg)

{
	if (NULL == bgrImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (bgrImg.channels() != 3 || grayImg.channels() != 1)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int width = bgrImg.cols;
	int height = bgrImg.rows;
	if (width != grayImg.cols || height != grayImg.rows)
	{
		cout << "image size  wrong." << endl;
		return -1;
	}
	int all = width*height;
	int b, g, r, grayVal;
	for (int line = 0; line < all; line++)
	{
		b = bgrImg.data[3 * line + 0];
		g = bgrImg.data[3 * line + 1];
		r = bgrImg.data[3 * line + 2];
		grayVal = (b * 117 + g * 601 + r * 306) >> 10;
		grayImg.data[line] = grayVal;
	}
	return 1;

}


int ustc_Threshold(Mat grayImg, Mat& binaryImg, int th)

{
	if (NULL == grayImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (th<0 || th>255)
	{
		cout << "th fail." << endl;
		return -1;
	}
	if (binaryImg.channels() != 1 || grayImg.channels() != 1)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	if (width != binaryImg.cols || height != binaryImg.rows)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int all = width*height;
	int pixlabel[256] = { 0 };
	for (int i = th; i < 256; i++)
		pixlabel[i] = 255;

	for (int line = 0; line < all; line++)
	{
		int pixVal = grayImg.data[line];
		binaryImg.data[line] = pixlabel[pixVal];
	}
	return 1;
}


int ustc_CalcGrad(Mat grayImg, Mat& gradImg_x, Mat& gradImg_y)
{
	if (NULL == grayImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (gradImg_x.channels() != 1 || grayImg.channels() != 1 || gradImg_y.channels() != 1)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	if (width != gradImg_x.cols || height != gradImg_x.rows || width != gradImg_y.cols || height != gradImg_y.rows)
	{
		cout << "image   wrong." << endl;
		return -1;
	}
	int height_true = height - 1;
	int width_true = width - 1;
	if (height_true <= 1 || width_true <= 1)
	{
		cout << "image size fail." << endl;
		return -1;
	}
	int num_x, num_y, grad_x, grad_y;
	gradImg_x.setTo(0);
	gradImg_y.setTo(0);

	for (int col_j = 1; col_j<width_true; col_j++)
	{
		for (int row_i = 1; row_i<height_true; row_i++)
		{
			num_x = row_i * width + col_j;
			grad_x =
				grayImg.data[num_x - width + 1]
				+ (grayImg.data[num_x + 1] << 1)
				+ grayImg.data[num_x + width + 1]
				- grayImg.data[num_x - width - 1]
				- (grayImg.data[num_x - 1] << 1)
				- grayImg.data[num_x + width - 1];

			((float*)gradImg_x.data)[num_x] = grad_x;
		}
	}

	for (int row_i = 1; row_i < height_true; row_i++)
	{
		for (int col_j = 1; col_j < width_true; col_j++)
		{
			num_y = row_i * width + col_j;
			grad_y =
				-grayImg.data[num_y - width + 1]
				+ grayImg.data[num_y + width + 1]
				- (grayImg.data[num_y - width] << 1)
				+ (grayImg.data[num_y + width] << 1)
				- grayImg.data[num_y - width - 1]
				+ grayImg.data[num_y + width - 1];

			((float*)gradImg_y.data)[num_y] = grad_y;
		}
	}
	return 1;
}


int ustc_CalcAngleMag(Mat gradImg_x, Mat gradImg_y, Mat& angleImg, Mat& magImg)
{
	if (NULL == gradImg_x.data || NULL == gradImg_y.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (gradImg_x.channels() != 1 || angleImg.channels() != 1 || gradImg_y.channels() != 1 || magImg.channels() != 1)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int width = gradImg_x.cols;
	int height = gradImg_x.rows;
	if (width != gradImg_y.cols || height != gradImg_y.rows || width != angleImg.cols || height != angleImg.rows || width != magImg.cols || height != magImg.rows)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	if (height <= 2 || width <= 2)
	{
		cout << "fail." << endl;
		return -1;
	}
	angleImg.setTo(0);
	magImg.setTo(0);
	
	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j += 1)
		{
			int numangle = row_i * width + col_j;
			float grad_x = ((float*)gradImg_x.data)[numangle];
			float grad_y = ((float*)gradImg_y.data)[numangle];
			float xielv = grad_y / grad_x;
			float angle = xielv - xielv*xielv*xielv / 3;
			angle = angle / 3.14 * 180;
			if (angle < 0)
				angle += 180;
			if (grad_y<0)
				angle += 180;
			
			((float*)angleImg.data)[numangle] = angle;
		}
	}

	
	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j += 1)
		{
			int nummag = row_i * width + col_j;
			float grad_x = ((float*)gradImg_x.data)[nummag];
			float grad_y = ((float*)gradImg_y.data)[nummag];
			float magval = grad_x*grad_x + grad_y*grad_y;
			float mag = sqrt(magval);
			((float*)magImg.data)[nummag] = mag;
		}
	}
	return 1;

}


int ustc_CalcHist(Mat grayImg, int* hist, int hist_len)
{
	if (NULL == grayImg.data || NULL == hist)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (256 != hist_len)
	{
		cout << "hist_len is wrong." << endl;
		return -1;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	int all = width*height;

	for (int i = 0; i < hist_len; i++)
	{
		hist[i] = 0;
	}


	for (int i = 0; i < all; i++)
	{
		int pixVal = grayImg.data[i];
		hist[pixVal]++;

	}
	return 1;
}


int ustc_SubImgMatch_gray(Mat grayImg, Mat subImg, int* x1, int* y1)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (grayImg.channels() != 1 || subImg.channels() != 1)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	int checkhang = height - sub_height;
	int checklie = width - sub_width;
	if (checkhang<0 || checklie<0)
	{
		cout << "subimage > grayimage,fail" << endl;
		return -1;
	}
	int v_size = width*height;

	Mat searchImg(height, width, CV_32FC1);
	searchImg.setTo(FLT_MAX);

	for (int i = 0; i <= checkhang; i++)
	{
		for (int j = 0; j <= checklie; j++)
		{
			int total_diff = 0;
			
			for (int xx = 0; xx < sub_height; xx++)
			{
				for (int yy = 0; yy < sub_width; yy++)
				{
					int height_index = i + xx;
					int width_index = j + yy;
					int bigImg_pix = grayImg.data[height_index * width + width_index];
					
					int template_pix = subImg.data[xx * sub_width + yy];
					int chazhi = bigImg_pix - template_pix;
					total_diff += ((chazhi ^ (chazhi >> 31)) - (chazhi >> 31));
				}
			}
			
			((float *)searchImg.data)[i * width + j] = total_diff;
		}
	}
	int min_num = 0;
	for (int i = 1; i < v_size; i++)
	{
		if (((float *)searchImg.data)[i]<((float *)searchImg.data)[min_num])
			min_num = i;

	}
	*y1 = min_num / width;
	*x1 = min_num%width;
	return 1;
}


int ustc_SubImgMatch_bgr(Mat colorImg, Mat subImg, int* x1, int* y1)
{
	if (NULL == colorImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (colorImg.channels() != 3 || subImg.channels() != 3)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int width = colorImg.cols;
	int height = colorImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	int checkhang = height - sub_height;
	int checklie = width - sub_width;
	if (checkhang<0 || checklie<0)
	{
		cout << "subimage > colorimage,fail" << endl;
		return -1;
	}
	int v_size = width*height;
	
	Mat searchImg(height, width, CV_32FC1);
	
	searchImg.setTo(FLT_MAX);


	for (int i = 0; i <= checkhang; i++)
	{
		for (int j = 0; j <= checklie; j++)
		{
			int total_diff = 0;
		
			for (int x = 0; x < sub_height; x++)
			{
				for (int y = 0; y < sub_width; y++)
				{
				
					int row_index = i + x;
					int col_index = j + y;
					int bignum = row_index * width + col_index;
					int tempnum = x * sub_width + y;
					int big_b = colorImg.data[3 * bignum + 0];
					int big_g = colorImg.data[3 * bignum + 1];
					int big_r = colorImg.data[3 * bignum + 2];
					int chazhi_b = subImg.data[3 * tempnum + 0] - big_b;
					int chazhi_g = subImg.data[3 * tempnum + 1] - big_g;
					int chazhi_r = subImg.data[3 * tempnum + 2] - big_r;

					total_diff += ((chazhi_b ^ (chazhi_b >> 31)) - (chazhi_b >> 31));
					total_diff += ((chazhi_g ^ (chazhi_g >> 31)) - (chazhi_g >> 31));
					total_diff += ((chazhi_r ^ (chazhi_r >> 31)) - (chazhi_r >> 31));
				}
			}
			
			((float*)searchImg.data)[i * width + j] = total_diff;
		}
	}
	int min_num = 0;
	for (int i = 1; i < v_size; i++)
	{
		if (((float *)searchImg.data)[i]<((float *)searchImg.data)[min_num])
			min_num = i;
	}
	*y1 = min_num / width;
	*x1 = min_num%width;
	return 1;
}


int ustc_SubImgMatch_corr(Mat grayImg, Mat subImg, int *x, int *y)
{
	
	if (NULL ==grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (grayImg.channels() != 1 || subImg.channels() != 1)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	
	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	int sub_hang = height - sub_height;
	int sub_lie = width - sub_width;
	if (sub_lie<0 || sub_hang<0)
	{
		cout << "size wrong" << endl;
		return -1;
	}

	int corr_i, corr_j;
	Mat searchImg(sub_hang, sub_lie, CV_32FC1);
	
	int v_size = sub_height*sub_width;
	for (corr_i = 0; corr_i <sub_hang; corr_i++)
		 {
		for (corr_j = 0; corr_j <sub_lie; corr_j++)
			 {
			 
				float xy = 0, x = 0, y = 0, xx = 0, yy = 0;
				for (int y1 = 0; y1 < sub_height; y1++)
				 {
					int row_index = corr_i + y1;
					int big_num = row_index * width;
					int sub_num = y1 * sub_width;
					for (int x1 = 0; x1 < sub_width; x1++)
					 {
					 
						int col_index = corr_j + x1;
						int bigImg_pix = grayImg.data[big_num + col_index];
						int template_pix = subImg.data[sub_num + x1];
						xy += bigImg_pix*template_pix;
						x += bigImg_pix;
						y += template_pix;
						xx += bigImg_pix*bigImg_pix;
						yy += template_pix*template_pix;
					 
					}
				}
				int fenzi = v_size*xy - x*y;
				int fenmu1 = sqrt(v_size*xx - x*x);
				int fenmu2 = sqrt(v_size*yy - y*y);
				int fenmu = fenmu1*fenmu2;
				int total_corr = fenzi / fenmu;
			
				((float*)searchImg.data)[corr_i *sub_lie + corr_j] = total_corr;
			 
			 }
		}
	int close1_x = 0, close1_y = 0;
	float maxcorr = ((float*)searchImg.data)[0];
	for (int i = 0; i < sub_hang; i++)
		 {
		for (int j = 0; j < sub_lie; j++)
		 {
			float a = ((float*)searchImg.data)[i * sub_lie + j];
			if (abs(a - 1) <abs(maxcorr - 1)) {
				maxcorr = a;
				close1_x = j;
				close1_y = i;
				
			}
			}
		}
	*y = close1_y;
	*x =close1_x;
	return 1;
}

int ustc_SubImgMatch_angle(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	if (grayImg.channels() != 1 || subImg.channels() != 1)
	{
		cout << "image channel  wrong." << endl;
		return -1;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	int subwidth = subImg.cols;
	int subheight = subImg.rows;

	Mat gradImg_x(height, width, CV_32FC1);
	Mat gradImg_y(height, width, CV_32FC1);
	Mat angleImg(height, width, CV_32FC1);

	Mat subgradImg_x(subheight, subwidth, CV_32FC1);
	Mat subgradImg_y(subheight, subwidth, CV_32FC1);
	Mat subangleImg(subheight, subwidth, CV_32FC1);

	int height_true = height - 1;
	int width_true = width - 1;
	if (height_true <= 1 || width_true <= 1)
	{
		cout << "image size fail." << endl;
		return -1;
	}
	int num_x, num_y, grad_x, grad_y;
	gradImg_x.setTo(0);
	gradImg_y.setTo(0);
	for (int col_j = 1; col_j<width_true; col_j++)
	{
		for (int row_i = 1; row_i<height_true; row_i++)
		{
			num_x = row_i * width + col_j;
			grad_x =
				grayImg.data[num_x - width + 1]
				+ (grayImg.data[num_x + 1] << 1)
				+ grayImg.data[num_x + width + 1]
				- grayImg.data[num_x - width - 1]
				- (grayImg.data[num_x - 1] << 1)
				- grayImg.data[num_x + width - 1];

			((float*)gradImg_x.data)[num_x] = grad_x;
		}
	}
	for (int row_i = 1; row_i < height_true; row_i++)
	{
		for (int col_j = 1; col_j < width_true; col_j++)
		{
			num_y = row_i * width + col_j;
			grad_y =
				-grayImg.data[num_y - width + 1]
				+ grayImg.data[num_y + width + 1]
				- (grayImg.data[num_y - width] << 1)
				+ (grayImg.data[num_y + width] << 1)
				- grayImg.data[num_y - width - 1]
				+ grayImg.data[num_y + width - 1];

			((float*)gradImg_y.data)[num_y] = grad_y;
		}
	}

	angleImg.setTo(0);
	for (int row_i = 1; row_i < height - 1; row_i++)
	{
		for (int col_j = 1; col_j < width - 1; col_j += 1)
		{
			int numangle = row_i * width + col_j;
			float grad_x_duru = ((float*)gradImg_x.data)[numangle];
			float grad_y_duru = ((float*)gradImg_y.data)[numangle];
			float xielv = grad_y_duru / grad_x_duru;
			float angle_duru = xielv - xielv*xielv*xielv / 3;
			angle_duru = angle_duru / 3.14 * 180;
			if (angle_duru < 0)
				angle_duru += 180;
			if (grad_y_duru<0)
				angle_duru += 180;
			((float*)angleImg.data)[numangle] = angle_duru;
		}
	}


	int subheight_true = subheight - 1;
	int subwidth_true = subwidth - 1;
	if (subheight_true <= 1 || subwidth_true <= 1)
	{
		cout << "image size fail." << endl;
		return -1;
	}
	int subnum_x, subnum_y, subgrad_x, subgrad_y;
	subgradImg_x.setTo(0);
	subgradImg_y.setTo(0);
	for (int col_j = 1; col_j<subwidth_true; col_j++)
	{
		for (int row_i = 1; row_i<subheight_true; row_i++)
		{
			subnum_x = row_i * subwidth + col_j;
			subgrad_x =
				subImg.data[subnum_x - subwidth + 1]
				+ (subImg.data[subnum_x + 1] << 1)
				+ subImg.data[subnum_x + subwidth + 1]
				- subImg.data[subnum_x - subwidth - 1]
				- (subImg.data[subnum_x - 1] << 1)
				- subImg.data[subnum_x + subwidth - 1];

			((float*)subgradImg_x.data)[subnum_x] = subgrad_x;
		}
	}
	for (int row_i = 1; row_i < subheight_true; row_i++)
	{
		for (int col_j = 1; col_j < subwidth_true; col_j++)
		{
			subnum_y = row_i * subwidth + col_j;
			subgrad_y =
				-subImg.data[subnum_y - subwidth + 1]
				+ subImg.data[subnum_y + subwidth + 1]
				- (subImg.data[subnum_y - subwidth] << 1)
				+ (subImg.data[subnum_y + subwidth] << 1)
				- subImg.data[subnum_y - subwidth - 1]
				+ subImg.data[subnum_y + subwidth - 1];

			((float*)subgradImg_y.data)[subnum_y] = subgrad_y;
		}
	}
	subangleImg.setTo(0);
	for (int row_i = 1; row_i < subheight - 1; row_i++)
	{
		for (int col_j = 1; col_j < subwidth - 1; col_j += 1)
		{
			int subnumangle = row_i * subwidth + col_j;
			float grad_x_sub = ((float*)subgradImg_x.data)[subnumangle];
			float grad_y_sub = ((float*)subgradImg_y.data)[subnumangle];
			float xielvsub = grad_y_sub / grad_x_sub;
			float angle_sub = xielvsub - xielvsub*xielvsub*xielvsub / 3;
			angle_sub = angle_sub / 3.14 * 180;
			if (angle_sub < 0)
				angle_sub += 180;
			if (grad_y_sub<0)
				angle_sub += 180;
			((float*)subangleImg.data)[subnumangle] = angle_sub;
		}
	}

	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	int checkhang = height - sub_height;
	int checklie = width - sub_width;
	if (checkhang<0 || checklie<0)
	{
		cout << "subimage > grayimage,fail" << endl;
		return -1;
	}
	int v_size = width*height;

	Mat searchImg(height, width, CV_32FC1);

	searchImg.setTo(FLT_MAX);

	for (int i = 0; i <= checkhang; i++)
	{
		for (int j = 0; j <= checklie; j++)
		{
			int total_diff = 0;

			for (int x = 1; x < sub_height - 1; x++)
			{
				for (int y = 1; y < sub_width - 1; y++)
				{
					int row_index = i + x;
					int col_index = j + y;
					float bigImg_pix = ((float*)angleImg.data)[row_index * width + col_index];
					float template_pix = ((float*)subangleImg.data)[x * sub_width + y];
					int chazhi = bigImg_pix - template_pix;
					total_diff += ((chazhi ^ (chazhi >> 31)) - (chazhi >> 31));
				}
			}

			((float *)searchImg.data)[i * width + j] = total_diff;
		}
	}
	int min_num = 0;
	for (int i = 1; i < v_size; i++)
	{
		if (((float *)searchImg.data)[i]<((float *)searchImg.data)[min_num])
			min_num = i;

	}
	*y= min_num / width;
	*x = min_num%width;
	return 1;
}

int ustc_SubImgMatch_mag(Mat grayImg, Mat subImg, int* x, int* y)
{

	if (NULL == grayImg.data || NULL == subImg.data)
	{
	cout << "image is NULL." << endl;
	return -1;
	}
	if (grayImg.channels() != 1 || subImg.channels() != 1)
	{
	cout << "image channel  wrong." << endl;
	return -1;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	int subwidth = subImg.cols;
	int subheight = subImg.rows;

	Mat gradImg_x(height, width, CV_32FC1);
	Mat gradImg_y(height, width, CV_32FC1);
	Mat magImg(height, width, CV_32FC1);

	Mat subgradImg_x(subheight, subwidth, CV_32FC1);
	Mat subgradImg_y(subheight, subwidth, CV_32FC1);
	Mat submagImg(subheight, subwidth, CV_32FC1);

	int height_true = height - 1;
	int width_true = width - 1;
	if (height_true <= 1 || width_true <= 1)
	{
	cout << "image size fail." << endl;
	return -1;
	}
	int num_x, num_y, grad_x, grad_y;
	gradImg_x.setTo(0);
	gradImg_y.setTo(0);
	for (int col_j = 1; col_j<width_true; col_j++)
	{
	for (int row_i = 1; row_i<height_true; row_i++)
	{
	num_x = row_i * width + col_j;
	grad_x =
	grayImg.data[num_x - width + 1]
	+ (grayImg.data[num_x + 1] << 1)
	+ grayImg.data[num_x + width + 1]
	- grayImg.data[num_x - width - 1]
	- (grayImg.data[num_x - 1] << 1)
	- grayImg.data[num_x + width - 1];

	((float*)gradImg_x.data)[num_x] = grad_x;
	}
	}
	for (int row_i = 1; row_i < height_true; row_i++)
	{
	for (int col_j = 1; col_j < width_true; col_j++)
	{
	num_y = row_i * width + col_j;
	grad_y =
	-grayImg.data[num_y - width + 1]
	+ grayImg.data[num_y + width + 1]
	- (grayImg.data[num_y - width] << 1)
	+ (grayImg.data[num_y + width] << 1)
	- grayImg.data[num_y - width - 1]
	+ grayImg.data[num_y + width - 1];

	((float*)gradImg_y.data)[num_y] = grad_y;
	}
	}

	magImg.setTo(0);
	for (int row_i = 1; row_i < height - 1; row_i++)
	{
	for (int col_j = 1; col_j < width - 1; col_j += 1)
	{
	int nummag = row_i * width + col_j;
	float grad_x_duru = ((float*)gradImg_x.data)[nummag];
	float grad_y_duru = ((float*)gradImg_y.data)[nummag];
	float magval = grad_x_duru*grad_x_duru + grad_y_duru*grad_y_duru;
	float mag = sqrt(magval);
	((float*)magImg.data)[nummag] = mag;
	}
	}

	int subheight_true = subheight - 1;
	int subwidth_true = subwidth - 1;
	if (subheight_true <= 1 || subwidth_true <= 1)
	{
	cout << "image size fail." << endl;
	return -1;
	}
	int subnum_x, subnum_y, subgrad_x, subgrad_y;
	subgradImg_x.setTo(0);
	subgradImg_y.setTo(0);
	for (int col_j = 1; col_j<subwidth_true; col_j++)
	{
	for (int row_i = 1; row_i<subheight_true; row_i++)
	{
	subnum_x = row_i * subwidth + col_j;
	subgrad_x =
	subImg.data[subnum_x - subwidth + 1]
	+ (subImg.data[subnum_x + 1] << 1)
	+ subImg.data[subnum_x + subwidth + 1]
	- subImg.data[subnum_x - subwidth - 1]
	- (subImg.data[subnum_x - 1] << 1)
	- subImg.data[subnum_x + subwidth - 1];

	((float*)subgradImg_x.data)[subnum_x] = subgrad_x;
	}
	}
	for (int row_i = 1; row_i < subheight_true; row_i++)
	{
	for (int col_j = 1; col_j < subwidth_true; col_j++)
	{
	subnum_y = row_i * subwidth + col_j;
	subgrad_y =
	-subImg.data[subnum_y - subwidth + 1]
	+ subImg.data[subnum_y + subwidth + 1]
	- (subImg.data[subnum_y - subwidth] << 1)
	+ (subImg.data[subnum_y + subwidth] << 1)
	- subImg.data[subnum_y - subwidth - 1]
	+ subImg.data[subnum_y + subwidth - 1];

	((float*)subgradImg_y.data)[subnum_y] = subgrad_y;
	}
	}
	submagImg.setTo(0);
	for (int row_i = 1; row_i < subheight - 1; row_i++)
	{
	for (int col_j = 1; col_j < subwidth - 1; col_j += 1)
	{
	int subnummag = row_i * subwidth + col_j;
	float grad_x_sub= ((float*)subgradImg_x.data)[subnummag];
	float grad_y_sub = ((float*)subgradImg_y.data)[subnummag];
	float submagval = grad_x_sub*grad_x_sub + grad_y_sub*grad_y_sub;
	float submag = sqrt(submagval);
	((float*)submagImg.data)[subnummag] = submag;
	}
	}

	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	int checkhang = height - sub_height;
	int checklie = width - sub_width;
	if (checkhang<0 || checklie<0)
	{
	cout << "subimage > grayimage,fail" << endl;
	return -1;
	}
	int v_size = width*height;

	Mat searchImg(height, width, CV_32FC1);

	searchImg.setTo(FLT_MAX);

	for (int i = 0; i <= checkhang; i++)
	{
	for (int j = 0; j <= checklie; j++)
	{
	int total_diff = 0;

	for (int x = 1; x < sub_height - 1; x++)
	{
	for (int y = 1; y < sub_width - 1; y++)
	{
	int row_index = i + x;
	int col_index = j + y;
	float bigImg_pix = ((float*)magImg.data)[row_index * width + col_index];
	float template_pix = ((float*)submagImg.data)[x * sub_width + y];
	int chazhi = bigImg_pix - template_pix;
	total_diff += ((chazhi ^ (chazhi >> 31)) - (chazhi >> 31));
	}
	}

	((float *)searchImg.data)[i * width + j] = total_diff;
	}
	}
	int min_num = 0;
	for (int i = 1; i < v_size; i++)
	{
	if (((float *)searchImg.data)[i]<((float *)searchImg.data)[min_num])
	min_num = i;

	}
	*y = min_num / width;
	*x = min_num%width;
	return 1;
}


int ustc_SubImgMatch_hist(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		cout << "image is NULL." << endl;
		return -1;
	}
	int width = grayImg.cols;
	int height = grayImg.rows;
	int sub_width = subImg.cols;
	int sub_height = subImg.rows;
	int checkhang = height - sub_height;
	int checklie = width - sub_width;
	if (checkhang<0 || checklie<0)
	{
		cout << "subimage > colorimage,fail" << endl;
		return -1;
	}
	int v_size = width*height;
	int hist_len = 256;
	int* hist_sub = new int[hist_len];
	memset(hist_sub, 0, sizeof(int) * hist_len);
	int subflag = ustc_CalcHist(subImg, hist_sub, hist_len);
	if (subflag<0)
	{
		cout << "sub hist fail." << endl;
		return -1;
	}
	
	Mat searchImg(height, width, CV_32FC1);
	
	searchImg.setTo(FLT_MAX);
	
	int* hist_temp = new int[hist_len];
	memset(hist_temp, 0, sizeof(int) * hist_len);

	for (int i = 0; i < checkhang; i++)
	{
		for (int j = 0; j <checklie; j++)
		{
			
			memset(hist_temp, 0, sizeof(int) * hist_len);

			
			for (int x = 0; x < sub_height; x++)
			{
				for (int y = 0; y < sub_width; y++)
				{
					
					int row_index = i + x;
					int col_index = j + y;
					int bigImg_pix = grayImg.data[row_index * width + col_index];
					hist_temp[bigImg_pix]++;
				}
			}
			
			int total_diff = 0;
			for (int ii = 0; ii < hist_len; ii++)
			{
				int chazhi = hist_temp[ii] - hist_sub[ii];
				total_diff += ((chazhi ^ (chazhi >> 31)) - (chazhi >> 31));
			}
			
			((float*)searchImg.data)[i * width + j] = total_diff;
		}
	}
	delete[] hist_temp;
	delete[] hist_sub;
	int min_num = 0;
	for (int i = 1; i < v_size; i++)
	{
		if (((float *)searchImg.data)[i]<((float *)searchImg.data)[min_num])
			min_num = i;
	}
	*y = min_num / width;
	*x = min_num%width;
	return 1;

}
