
#include"SubImageMatch.h"


int ustc_ConvertBgr2Gray(Mat bgrImg, Mat& grayImg) {
	int i=0, psize;
	uchar *Oimagep, *temp;

	if (NULL == grayImg.data || NULL == bgrImg.data )
	{
		printf("image is NULL.");
		return SUB_IMAGE_MATCH_FAIL;
	}

	temp = bgrImg.data;
	Oimagep = grayImg.data;
	 psize = (bgrImg.rows) * (bgrImg.cols);
	 
	for (; i < psize-4; i+=4) {
		*(Oimagep + i)= (*temp * 15 + *(temp + 1) * 75 + *(temp + 2) * 38) >> 7;

		*(Oimagep + i + 1) = (*(temp+3) * 15 + *(temp + 4) * 75 + *(temp + 5) * 38) >> 7;

		*(Oimagep + i + 2) = (*(temp+6) * 15 + *(temp + 7) * 75 + *(temp + 8) * 38) >> 7;

		*(Oimagep + i + 3) = (*(temp+9) * 15 + *(temp + 10) * 75 + *(temp + 11) * 38) >> 7;

		temp = temp + 12;

	}
	
	for (; i < psize; i++) {
		*(Oimagep + i) = (*temp * 15 + *(temp + 1) * 75 + *(temp + 2) * 38) >> 7;
		temp = temp + 3;
	}
	return SUB_IMAGE_MATCH_OK;
}

int ustc_CalcGrad(Mat grayImg, Mat& gradImg_x, Mat& gradImg_y) {
	int i, j, row, col;
	uchar *grayImgp;
	float *gradImg_xp = (float*)gradImg_x.data;
	float *gradImg_yp = (float*)gradImg_y.data;

	
	if (NULL == grayImg.data || NULL == gradImg_x.data || NULL == gradImg_y.data)
	{
		printf("image is NULL.");
		return SUB_IMAGE_MATCH_FAIL;
	}

	grayImgp = grayImg.data;
	row = grayImg.rows;
	col = grayImg.cols;
	for (i = 1; i < row - 1; i++) {
		for (j = 1; j < col - 1; j++) {
			gradImg_xp[i*row + j] = grayImgp[(i - 1)*row + j + 1] + grayImgp[i*row + j + 1] * 2 + grayImgp[(i + 1)*row + j + 1] - grayImgp[(i - 1)*row + j - 1] - grayImgp[i*row + j - 1] * 2 - grayImgp[(i + 1)*row + j - 1];
			gradImg_yp[i*row + j] = grayImgp[(i + 1)*row + j - 1] + grayImgp[(i + 1)*row + j] * 2 + grayImgp[(i + 1)*row + j + 1] - grayImgp[(i - 1)*row + j - 1] - grayImgp[(i - 1)*row + j] * 2 - grayImgp[(i - 1)*row + j + 1];
		}
	}
	return SUB_IMAGE_MATCH_OK;
}

int ustc_CalcAngleMag(Mat gradImg_x, Mat gradImg_y, Mat& angleImg, Mat& magImg) {
	if (NULL == angleImg.data || NULL == gradImg_x.data || NULL == gradImg_y.data||NULL==magImg.data)
	{
		printf("image is NULL.");
		return SUB_IMAGE_MATCH_FAIL;
	}

	static float atan_table[231] = {
		-89.50178839, -89.49741833,	-89.49297094,	-89.48844413,	-89.48383577,	-89.47914363,	-89.47436539,	-89.46949868,	-89.46454101,	-89.45948981,	-89.45434241,	-89.44909602,	-89.44374777,	-89.43829467,	-89.43273359,	-89.4270613,	-89.42127443,	-89.41536948,	-89.40934279,	-89.40319055,	-89.39690881,	-89.39049342,	-89.38394009,	-89.37724431,	-89.37040139,	-89.36340642,	-89.35625429,	-89.34893962,	-89.34145682,	-89.33380003,	-89.3259631, -89.31793961,	-89.3097228,	-89.30130562,	-89.29268063,	-89.28384005,	-89.2747757,	-89.26547897,	-89.2559408,	-89.24615167,	-89.23610154,	-89.22577984,	-89.2151754,	-89.20427645,	-89.19307054,	-89.18154454,	-89.16968451,	-89.15747574,	-89.1449026,	-89.13194855,	-89.118596,	-89.10482629,	-89.09061955,
		-89.07595465,	-89.06080905,	-89.04515875,	-89.02897807,	-89.0122396,	-88.99491399,	-88.97696981,	-88.95837332,	-88.93908831,	-88.91907581,	-88.89829388,	-88.87669729,	-88.85423716,	-88.83086067,	-88.80651058,	-88.78112476,	-88.75463573,	-88.72696998,	-88.69804733,	-88.66778015,	-88.63607247,	-88.60281897,	-88.56790382,	-88.53119929,	-88.49256424,	-88.4518423,	-88.40885973,	-88.36342296,	-88.31531568,	-88.26429541,	-88.21008939,	-88.15238973,	-88.09084757,	-88.02506599,	-87.95459151,	-87.8789036,	-87.79740184,	-87.70938996,	-87.61405597,	-87.51044708,	-87.3974378,	-87.27368901,	-87.13759477,	-86.9872125,	-86.82016988,	-86.63353934,	-86.42366563,	-86.18592517,	-85.91438322,	-85.60129465,	-85.23635831,	-84.80557109,	-84.28940686,
		-83.65980825,	-82.87498365,	-81.86989765,	-80.53767779,	-78.69006753,	-75.96375653,	-71.56505118,	-63.43494882,	-45,	  0,	45, 63.43494882, 71.56505118,	     75.96375653,	78.69006753,	80.53767779,	81.86989765,	82.87498365,	83.65980825,	84.28940686,	84.80557109,	85.23635831,	85.60129465,	85.91438322,	86.18592517,	86.42366563,	86.63353934,	86.82016988,	86.9872125,	87.13759477,	87.27368901,	87.3974378,	87.51044708,	87.61405597,	87.70938996,	87.79740184,	87.8789036, 87.95459151,	88.02506599,	88.09084757,	88.15238973,	88.21008939,	88.26429541,	88.31531568,	88.36342296,	88.40885973,	88.4518423,	88.49256424,	88.53119929,	88.56790382,	88.60281897,	88.63607247,	88.66778015,	88.69804733,	88.72696998,	88.75463573,
		88.78112476,	88.80651058,	88.83086067,	88.85423716,	88.87669729,	88.89829388,	88.91907581,	88.93908831,	88.95837332,	88.97696981,	88.99491399,	89.0122396,	89.02897807,	89.04515875,	89.06080905,	89.07595465,	89.09061955,	89.10482629,	89.118596,	89.13194855,	89.1449026,	89.15747574,	89.16968451,	89.18154454,	89.19307054,	89.2042764,	89.2151754, 89.22577984,	89.23610154,	89.24615167,	89.2559408,	89.26547897,	89.2747757,	89.28384005,	89.29268063,	89.30130562,	89.3097228,	89.31793961,	89.3259631,	89.33380003,	89.34145682,	89.34893962,	89.35625429,	89.36340642,	89.37040139,	89.37724431,	89.38394009,	89.39049342,	89.39690881,	89.40319055,	89.40934279,	89.41536948,	89.42127443,	89.4270613,	89.43273359,	89.43829467,
		89.44374777,	89.44909602,	89.45434241,	89.45948981,	89.46454101,	89.46949868,	89.47436539,	89.47914363,	89.48383577,	89.48844413,	89.49297094,	89.49741833,	89.50178839
	};



	int i, psize;
	float x, y, tg, Theta, temp;
	float *gradImg_xp = (float*)gradImg_x.data;
	float *gradImg_yp = (float*)gradImg_y.data;
	float *angleImgp = (float*)angleImg.data;
	float *magImgp = (float*)magImg.data;
	psize = (gradImg_x.rows) * (gradImg_x.cols);
	
	for (i = 0; i < psize; i++) {
		y = gradImg_yp[i];
		x = gradImg_xp[i];
		
		if (x && y) {
			tg = y / x;
		}
		else {
			tg = 0;
		}
		float t = x*x + y*y;
		
		float xhalf = 0.5f*t;

		int k = *(int*)&t;

		k = 0x5f3759df - (k >> 1); 

		t = *(float*)&k;

		t = t*(1.5f - xhalf*t*t); 

		magImgp[i] = 1 / t;


		Theta = (atan_table[(int)(tg)+1+115] - atan_table[(int)(tg)+115])*(tg - (int)tg) + atan_table[(int)(tg)+115];
		
		if (x && y) {
			angleImgp[i] = Theta;
		}
		else if (x && ! y) {
			angleImgp[i] = Theta + 360;
		}
		else {
			angleImgp[i] = Theta + 180;
		}
}

	return SUB_IMAGE_MATCH_OK;
}

int ustc_Threshold(Mat grayImg, Mat& binaryImg, int th) {
	int i, psize;
	uchar *grayImgp = grayImg.data;
	uchar *binaryImgp = binaryImg.data;

	
	if (NULL == grayImg.data || NULL == binaryImg.data){
		printf("image is NULL.");
		return SUB_IMAGE_MATCH_FAIL;
	}

	psize = grayImg.rows * grayImg.cols;

	for (i = 0; i < psize ; i ++) {
		binaryImgp[i] = (grayImgp[i] >= th) * 255;
	}
	return 1;
}

int ustc_CalcHist(Mat grayImg, int* hist, int hist_len) {
	int i, psize;
	uchar* grayImgp = grayImg.data;

	
	if (NULL == grayImg.data || NULL == hist)
	{
		printf("image is NULL.");
		return -1;
	}

	
	for (int i = 0; i < hist_len; i++)
	{
		hist[i] = 0;
	}
	
	psize = grayImg.rows * grayImg.cols;
	for (i = 0; i < psize; i++) {
		if (grayImgp[i] <= hist_len) {
			hist[grayImgp[i]]++;
		}
	}
	return SUB_IMAGE_MATCH_OK;
}
int ustc_SubImgMatch_gray(Mat grayImg, Mat subImg, int* x, int* y) {
	if (NULL == grayImg.data || NULL == subImg.data) {
		printf("img is NULL.");
		return SUB_IMAGE_MATCH_FAIL;
	}
	
	int row = grayImg.rows;
	int col = grayImg.cols;
	int sub_row = subImg.rows;
	int sub_col = subImg.cols;

	if (sub_row > row || sub_col > col) {
		printf("subimage fail");
		return SUB_IMAGE_MATCH_FAIL;
	}

	int min_diff = 0x4fffffff;
	int bigImg_pix;
	int template_pix;
	int total_diff = 0;
	int yy;
	int big_temp;
	int sub_temp;
	uchar *big_p=grayImg.data, *sub_p=subImg.data;
	
	for (int i = row - sub_row; i>=0 ; --i) {
		for (int j = col - sub_col; j>=0 ; --j) {
			total_diff = 0;
			


			for (int xx = 0; xx<sub_row; ++xx) {

				yy = 0;
				for (; yy < sub_col - 4; yy += 4) {
					
					big_temp = (i+xx)*col + j + yy;
					sub_temp = xx*sub_col + yy;
					bigImg_pix = big_p[big_temp];
					
					template_pix = sub_p[ sub_temp];

					total_diff += ((bigImg_pix - template_pix) ^ ((bigImg_pix - template_pix) >> 31)) - ((bigImg_pix - template_pix) >> 31);


					bigImg_pix = big_p[big_temp+ 1];
					
					template_pix = sub_p[sub_temp + 1];

					total_diff += ((bigImg_pix - template_pix) ^ ((bigImg_pix - template_pix) >> 31)) - ((bigImg_pix - template_pix) >> 31);


					bigImg_pix = big_p[big_temp+ 2];
					
					template_pix = sub_p[sub_temp + 2];

					total_diff += ((bigImg_pix - template_pix) ^ ((bigImg_pix - template_pix) >> 31)) - ((bigImg_pix - template_pix) >> 31);


					bigImg_pix = big_p[big_temp + 3];
					
					template_pix = sub_p[sub_temp + 3];

					total_diff += ((bigImg_pix - template_pix) ^ ((bigImg_pix - template_pix) >> 31)) - ((bigImg_pix - template_pix) >> 31);

				}
				for (; yy < sub_col; yy += 1) {
					
					

					bigImg_pix = big_p[(i + xx)*col + j+yy];
					
					template_pix = sub_p[xx*sub_col + yy];

					total_diff += ((bigImg_pix - template_pix) ^ ((bigImg_pix - template_pix) >> 31)) - ((bigImg_pix - template_pix) >> 31);

				}
				
				if (total_diff < min_diff) {
					*x = j;
					*y = i;
					min_diff = total_diff;
				}
			}
		}
	}
	return SUB_IMAGE_MATCH_OK;
}
int ustc_SubImgMatch_bgr(Mat colorImg, Mat subImg, int* x, int* y) {
	if (NULL == colorImg.data || NULL == subImg.data) {
		printf("img is NULL.");
		return SUB_IMAGE_MATCH_FAIL;
	}
	int xx, yy, i, j;
	int row = colorImg.rows;
	int col = colorImg.cols;
	int sub_row = subImg.rows;
	int sub_col = subImg.cols;

	if (sub_row > row || sub_col > col) {
		printf("subimage fail");
		return SUB_IMAGE_MATCH_FAIL;
	}

	int min_diff = 0x4fffffff;
	int big_b, big_g,big_r, big_bgr;
	int sub_b,sub_g, sub_r, sub_bgr;
	int total_diff;
	int big_temp;
	int sub_temp;
	uchar* sub_p, *big_p;
	sub_p = subImg.data;
	big_p = colorImg.data;
	
	for ( i = row - sub_row; i >= 0; --i) {

		for (j = col - sub_col; j>= 0; --j) {
			total_diff = 0;
			
			for (xx = 0; xx < sub_row; ++xx) {
				yy = 0;
				for (; yy < sub_col- 4; yy+=4) {
					
					
					
                    big_temp = 3 * ((i+xx)*col + j+yy);

                    sub_temp = 3 * (xx*sub_col + yy);

                   
					sub_b = sub_p[sub_temp + 0];
					sub_g = sub_p[sub_temp + 1];
					sub_r = sub_p[sub_temp + 2];

					big_b = big_p[big_temp + 0];
					big_g = big_p[big_temp + 1];
					big_r = big_p[big_temp + 2];
		
					total_diff += ((((big_b - sub_b) ^ ((big_b - sub_b) >> 31)) - ((big_b - sub_b) >> 31)) + (((big_g - sub_g) ^ ((big_g - sub_g) >> 31)) - ((big_g - sub_g) >> 31)) + (((big_r - sub_r) ^ ((big_r - sub_r) >> 31)) - ((big_r - sub_r) >> 31)));
		
					sub_b = sub_p[sub_temp + 3];
					sub_g = sub_p[sub_temp + 4];
					sub_r = sub_p[sub_temp + 5];

					big_b = big_p[big_temp + 3];
					big_g = big_p[big_temp + 4];
					big_r = big_p[big_temp + 5];

					total_diff += ((((big_b - sub_b) ^ ((big_b - sub_b) >> 31)) - ((big_b - sub_b) >> 31)) + (((big_g - sub_g) ^ ((big_g - sub_g) >> 31)) - ((big_g - sub_g) >> 31)) + (((big_r - sub_r) ^ ((big_r - sub_r) >> 31)) - ((big_r - sub_r) >> 31)));

					sub_b = sub_p[sub_temp + 6];
					sub_g = sub_p[sub_temp + 7];
					sub_r = sub_p[sub_temp + 8];

					big_b = big_p[big_temp + 6];
					big_g = big_p[big_temp + 7];
					big_r = big_p[big_temp + 8];

					total_diff += ((((big_b - sub_b) ^ ((big_b - sub_b) >> 31)) - ((big_b - sub_b) >> 31)) + (((big_g - sub_g) ^ ((big_g - sub_g) >> 31)) - ((big_g - sub_g) >> 31)) + (((big_r - sub_r) ^ ((big_r - sub_r) >> 31)) - ((big_r - sub_r) >> 31)));

					sub_b = sub_p[sub_temp + 9];
					sub_g = sub_p[sub_temp + 10];
					sub_r = sub_p[sub_temp + 11];

					big_b = big_p[big_temp + 9];
					big_g = big_p[big_temp + 10];
					big_r = big_p[big_temp + 11];

					total_diff += ((((big_b - sub_b) ^ ((big_b - sub_b) >> 31)) - ((big_b - sub_b) >> 31)) + (((big_g - sub_g) ^ ((big_g - sub_g) >> 31)) - ((big_g - sub_g) >> 31)) + (((big_r - sub_r) ^ ((big_r - sub_r) >> 31)) - ((big_r - sub_r) >> 31)));

					for (; yy < sub_col ; yy += 1) {
						



						big_temp = 3 * ((i + xx)*col + j + yy);

						sub_temp = 3 * (xx*sub_col + yy);


						sub_b = sub_p[sub_temp + 0];
						sub_g = sub_p[sub_temp + 1];
						sub_r = sub_p[sub_temp + 2];

						big_b = big_p[big_temp + 0];
						big_g = big_p[big_temp + 1];
						big_r = big_p[big_temp + 2];

						total_diff += ((((big_b - sub_b) ^ ((big_b - sub_b) >> 31)) - ((big_b - sub_b) >> 31)) + (((big_g - sub_g) ^ ((big_g - sub_g) >> 31)) - ((big_g - sub_g) >> 31)) + (((big_r - sub_r) ^ ((big_r - sub_r) >> 31)) - ((big_r - sub_r) >> 31)));
					}
				}
			}
			
			if (total_diff < min_diff) {
				*x = j;
				*y = i;
				min_diff = total_diff;
			}
		}
	}
	return SUB_IMAGE_MATCH_OK;
}


int ustc_SubImgMatch_corr(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		printf( "image is NULL." );
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (grayImg.rows < subImg.rows || grayImg.cols < subImg.cols)
	{
		printf("image is false." );
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (x == NULL || NULL == y)
	{
		x = new int(0);
		y = new int(0);
	}
	int row_f = grayImg.rows;
	int col_f = grayImg.cols;
	int row_s = subImg.rows;
	int col_s = subImg.cols;
	int fx = 0, fy = 0;
	long int gray_sqsum = 0;
	long int sub_sqsum = 0;
	long int gray_sub_sum = 0;
	float corr = 0;
	int diff[2] = { 0 };
	int rowlim = row_f - row_s + 1;
	int collim = col_f - col_s + 1;
	int i_f, j_f, i_sub, j_sub;
	uchar *pgray = grayImg.data;
	uchar *psub = subImg.data;
	int res = col_s % 4, mod4j = col_s - res;
	for (i_sub = 0; i_sub < row_s; i_sub++)
	{
		int xis = i_sub*col_s;
		for (j_sub = 0; j_sub < col_s; j_sub++)
		{
			sub_sqsum += psub[xis + j_sub] * psub[xis + j_sub];
		}
	}
	for (i_f = 0; i_f < rowlim; i_f++)
	{

		for (j_f = 0; j_f < collim; j_f++)
		{
			corr = 0;
			gray_sqsum = 0;
			gray_sub_sum = 0;
			for (i_sub = 0; i_sub < row_s; i_sub++)
			{
				int xif = (i_f + i_sub)*col_f + j_f;
				int xis1 = i_sub*col_s;
				for (j_sub = 0; j_sub < mod4j; j_sub += 4)
				{
					int numf = xif + j_sub;
					int nums = xis1 + j_sub;
					int pgray_ij = pgray[numf + 0];
					gray_sqsum += pgray_ij * pgray_ij;
					gray_sub_sum += pgray_ij * psub[nums + 0];

					pgray_ij = pgray[numf + 1];
					gray_sqsum += pgray_ij * pgray_ij;
					gray_sub_sum += pgray_ij * psub[nums + 1];

					pgray_ij = pgray[numf + 2];
					gray_sqsum += pgray_ij * pgray_ij;
					gray_sub_sum += pgray_ij * psub[nums + 2];

					pgray_ij = pgray[numf + 3];
					gray_sqsum += pgray_ij * pgray_ij;
					gray_sub_sum += pgray_ij * psub[nums + 3];

					
				}
				for (; j_sub < res; j_sub++)
				{
					gray_sqsum += pgray[xif + j_sub] * pgray[xif + j_sub];
					gray_sub_sum += pgray[xif + j_sub] * psub[xis1 + j_sub];
				}
			}
			corr =(gray_sub_sum) / (sqrt(gray_sqsum)*sqrt(sub_sqsum));
			diff[0] = corr * 1000;
			if (diff[0] > diff[1])
			{
				diff[1] = diff[0];
				fx = j_f;
				fy = i_f;
			}
		}
	}
	if (diff[1] == 0)
		return SUB_IMAGE_MATCH_FAIL;

	*x = fx;
	*y = fy;
	if (i_f == rowlim && j_f == collim)
		return SUB_IMAGE_MATCH_OK;
	else
		return SUB_IMAGE_MATCH_FAIL;
}



int ustc_SubImgMatch_angle(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		printf("image is NULL."); 
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (grayImg.rows < subImg.rows || grayImg.cols < subImg.cols)
	{
		printf("image is false.");
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (x == NULL || NULL == y)
	{
		x = new int(0);
		y = new int(0);
	}
	int i_f, j_f, i_sub, j_sub;

	
	Mat gray_x = Mat::zeros(grayImg.size(), CV_32FC1);
	Mat gray_y = Mat::zeros(grayImg.size(), CV_32FC1);
	Mat gray_angle = Mat::zeros(grayImg.size(), CV_32FC1);
	Mat gray_mag = Mat::zeros(grayImg.size(), CV_32FC1);
	ustc_CalcGrad(grayImg, gray_x, gray_y);
	ustc_CalcAngleMag(gray_x, gray_y, gray_angle, gray_mag);

	Mat sub_x = Mat::zeros(subImg.size(), CV_32FC1);
	Mat sub_y = Mat::zeros(subImg.size(), CV_32FC1);
	Mat sub_angle = Mat::zeros(subImg.size(), CV_32FC1);
	Mat sub_mag = Mat::zeros(subImg.size(), CV_32FC1);
	ustc_CalcGrad(subImg, sub_x, sub_y);
	ustc_CalcAngleMag(sub_x, sub_y, sub_angle, sub_mag);

	int row_f = grayImg.rows;
	int col_f = grayImg.cols;
	int row_s = subImg.rows;
	int col_s = subImg.cols;
	int diff[2] = { 0,2147483647 };
	int sum = 0;
	int fx = 0, fy = 0;
	float *pgray = (float*)gray_angle.data;
	float *psub = (float*)sub_angle.data;
	int rowlim = row_f - row_s + 1;
	int collim = col_f - col_s + 1;
	int res = col_s % 8, mod5j = col_s - res;
	//match
	for (i_f = 0; i_f < rowlim; i_f++)
	{

		for (j_f = 0; j_f < collim; j_f++)
		{
			diff[0] = 0;
			for (i_sub = 1; i_sub < row_s - 1; i_sub++)
			{
				int xif = (i_f + i_sub)*col_f + j_f;
				int xis = i_sub*col_s;
				for (j_sub = 1; j_sub < mod5j; j_sub = 8 + j_sub)
				{
					int numf = xif + j_sub;
					int nums = xis + j_sub;
					sum = pgray[numf + 0] - psub[nums + 0];
					sum = (sum >> 31) *(sum << 1) + sum;
					int sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;

					sum = pgray[numf + 1] - psub[nums + 1];
					sum = (sum >> 31) *(sum << 1) + sum;
					sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;

					sum = pgray[numf + 2] - psub[nums + 2];
					sum = (sum >> 31) *(sum << 1) + sum;
					sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;

					sum = pgray[numf + 3] - psub[nums + 3];
					sum = (sum >> 31) *(sum << 1) + sum;
					sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;

					sum = pgray[numf + 4] - psub[nums + 4];
					sum = (sum >> 31) *(sum << 1) + sum;
					sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;

					sum = (int)pgray[numf + 5] - (int)psub[nums + 5];
					sum = (sum >> 31) *(sum << 1) + sum;
					sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;

					sum = pgray[numf + 6] - psub[nums + 6];
					sum = (sum >> 31) *(sum << 1) + sum;
					sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;

					sum = pgray[numf + 7] - psub[nums + 7];
					sum = (sum >> 31) *(sum << 1) + sum;
					sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;
				}
				for (; j_sub < res; j_sub++)
				{
					sum = (int)pgray[xif + j_sub] - (int)psub[xis + j_sub];
					sum = (sum >> 31) *(sum << 1) + sum;
					int sum0 = 180 - sum;
					sum = (sum0 >> 31) *((sum << 1) - 360) + sum;
					diff[0] += sum;
				}
			}
			if (diff[0] < diff[1])
			{
				diff[1] = diff[0];
				fx = j_f;
				fy = i_f;
			}
		}
	}
	if (diff[1] == 2147483647)
		return SUB_IMAGE_MATCH_FAIL;

	*x = fx;
	*y = fy;
	if (i_f == rowlim && j_f == collim)
		return SUB_IMAGE_MATCH_OK;
	else
		return SUB_IMAGE_MATCH_FAIL;
}

int ustc_SubImgMatch_mag(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		printf( "image is NULL." );
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (grayImg.rows < subImg.rows || grayImg.cols < subImg.cols)
	{
		printf("image is false." );
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (x == NULL || NULL == y)
	{
		x = new int(0);
		y = new int(0);
	}
	int i_f, j_f, i_sub, j_sub;

	
	Mat gray_x = Mat::zeros(grayImg.size(), CV_32FC1);
	Mat gray_y = Mat::zeros(grayImg.size(), CV_32FC1);
	Mat gray_angle = Mat::zeros(grayImg.size(), CV_32FC1);
	Mat gray_mag = Mat::zeros(grayImg.size(), CV_32FC1);
	ustc_CalcGrad(grayImg, gray_x, gray_y);
	ustc_CalcAngleMag(gray_x, gray_y, gray_angle, gray_mag);

	Mat sub_x = Mat::zeros(subImg.size(), CV_32FC1);
	Mat sub_y = Mat::zeros(subImg.size(), CV_32FC1);
	Mat sub_angle = Mat::zeros(subImg.size(), CV_32FC1);
	Mat sub_mag = Mat::zeros(subImg.size(), CV_32FC1);
	ustc_CalcGrad(subImg, sub_x, sub_y);
	ustc_CalcAngleMag(sub_x, sub_y, sub_angle, sub_mag);

	int row_f = grayImg.rows;
	int col_f = grayImg.cols;
	int row_s = subImg.rows;
	int col_s = subImg.cols;
	int diff[2] = { 0,2147483647 };
	int sum = 0;
	int fx = 0, fy = 0;
	float *pgray = (float*)gray_mag.data;
	float *psub = (float*)sub_mag.data;
	int rowlim = row_f - row_s + 1;
	int collim = col_f - col_s + 1;
	int res = col_s % 5, mod5j = col_s - res;
	//match
	for (i_f = 0; i_f < rowlim; i_f++)
	{

		for (j_f = 0; j_f < collim; j_f++)
		{
			diff[0] = 0;
			for (i_sub = 0; i_sub < row_s - 1; i_sub++)
			{
				int xif = (i_f + i_sub)*col_f + j_f;
				int xis = i_sub*col_s;
				for (j_sub = 0; j_sub < mod5j; j_sub += 5)
				{
					int numf = xif + j_sub;
					int nums = xis + j_sub;
					sum = (int)pgray[numf + 0] - psub[nums + 0];
					sum = (sum >> 31)*(sum << 1) + sum;
					diff[0] += sum;

					sum = (int)pgray[numf + 1] - psub[nums + 1];
					sum = (sum >> 31)*(sum << 1) + sum;
					diff[0] += sum;

					sum = (int)pgray[numf + 2] - psub[nums + 2];
					sum = (sum >> 31)*(sum << 1) + sum;
					diff[0] += sum;

					sum = (int)pgray[numf + 3] - psub[nums + 3];
					sum = (sum >> 31)*(sum << 1) + sum;
					diff[0] += sum;

					sum = (int)pgray[numf + 4] - psub[nums + 4];
					sum = (sum >> 31)*(sum << 1) + sum;
					diff[0] += sum;
				}
				for (; j_sub < res; j_sub++)
				{
					sum = (int)pgray[xif + j_sub] - psub[xis + j_sub];
					sum = (sum >> 31)*(sum << 1) + sum;
					diff[0] += sum;
				}
			}
			if (diff[0] < diff[1])
			{
				diff[1] = diff[0];
				fx = j_f;
				fy = i_f;
			}
		}
	}
	if (diff[1] == 2147483647)
		return SUB_IMAGE_MATCH_FAIL;

	*x = fx;
	*y = fy;
	if (i_f == rowlim && j_f == collim)
		return SUB_IMAGE_MATCH_OK;
	else
		return SUB_IMAGE_MATCH_FAIL;
}



int ustc_SubImgMatch_hist(Mat grayImg, Mat subImg, int* x, int* y)
{
	if (NULL == grayImg.data || NULL == subImg.data)
	{
		printf("image is NULL.") ;
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (grayImg.rows < subImg.rows || grayImg.cols < subImg.cols)
	{
		printf("image is false." );
		return SUB_IMAGE_MATCH_FAIL;
	}
	if (x == NULL || NULL == y)
	{
		x = new int(0);
		y = new int(0);
	}
	uchar *pgray = grayImg.data;
	int row_f = grayImg.rows;
	int col_f = grayImg.cols;
	int row_s = subImg.rows;
	int col_s = subImg.cols;
	int size_sub = row_s*col_s;
	int rowlim = row_f - row_s + 1;
	int collim = col_f - col_s + 1;
	int subHist[256] = { 0 };
	int diff[2] = { 0,2147483647 };
	int sum, fx = 0, fy = 0;
	ustc_CalcHist(subImg, subHist, 256);
	int i_f, j_f, i_sub, j_sub;
	for (i_f = 0; i_f < rowlim; i_f++)
	{
		for (j_f = 0; j_f < collim; j_f++)
		{
			diff[0] = 0;
			int grayHist[256] = { 0 };
			for (i_sub = 0; i_sub < row_s; i_sub++)
			{
				int xif = (i_f + i_sub)*col_f + j_f;
				int xis = i_sub*col_s;
				for (j_sub = 0; j_sub < col_s; j_sub++)
				{
					grayHist[pgray[xif + j_sub]]++;
				}
			}
			for (i_sub = 0; i_sub < 256; i_sub++)
			{
				sum = grayHist[i_sub] - subHist[i_sub];
				sum = (sum >> 31)*(sum << 1) + sum;
				diff[0] += sum;
			}
			if (diff[0] < diff[1])
			{
				diff[1] = diff[0];
				fx = j_f;
				fy = i_f;
			}
		}
	}
	if (diff[1] == 2147483647)
		return SUB_IMAGE_MATCH_FAIL;
	*x = fx;
	*y = fy;
	if (i_f == rowlim&& j_f == collim)
		return SUB_IMAGE_MATCH_OK;
	else
		return SUB_IMAGE_MATCH_FAIL;
}




