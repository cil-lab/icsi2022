#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "unistd.h"
#include "string.h"
#include "icsi2022.h"

#define INF 1.0e99
#define EPS 1.0e-14
#define E 2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

double *O, *M, *y, *z;
int ini_flag, n_flag, func_flag, *SS;

void ceval(double *x, int nx, int mx, double *f, int func_num, char *path)
{
	// api for python
	// parameter:
	// nx: dimension of x
	// mx: row size
	//

	int f_num = 10, i, j;


	// check for function num and dimension of x
	if (func_num < 1 || func_num > 10)
	{
		return;
	}
	if (!(nx == 10 || nx == 20 || nx == 50))
	{
		return;
	}

	FILE *fpt;
	char Path[512];
	char FileName[256];
	free(M);
	free(O);
	free(y);
	free(z);
	y = (double *)malloc(sizeof(double) * nx);
	z = (double *)malloc(sizeof(double) * nx);



	/* Load shift_data */

	strcpy(Path, path);
	sprintf(FileName, strcat(Path, "/data/S_%d_D%d.txt"), func_num, nx);

	fpt = fopen(FileName, "r");

	if (fpt == NULL)
	{
		printf("\nError: An error occurred while attempting open file %s!\n",FileName);
		return;
	}
	// for basic functions
	if (func_num < 8)
	{
		O = (double *)malloc(nx * sizeof(double));
		if (O == NULL)
		{
			return;

		}
		for (i = 0; i < nx; i++)
		{
			fscanf(fpt, "%lf", &O[i]);
		}
	}
	else
	{
		// load shift data for composition functions
		O = (double *)malloc(nx * f_num * sizeof(double));
		if (O == NULL)
		{
			return;
		}
		for (i = 0; i < f_num - 1; i++)
		{
			for (j = 0; j < nx; j++)
			{
				fscanf(fpt, "%lf", &O[i * nx + j]);
			}
			fscanf(fpt, "%*[^\n]%*c");
		}
		for (j = 0; j < nx; j++)
		{
			fscanf(fpt, "%lf", &O[nx * (f_num - 1) + j]);
		}
	}
	fclose(fpt);

	/* Load Rotate Matrix M*/
	strcpy(Path, path);
	sprintf(FileName, strcat(Path, "/data/M_%d_D%d.txt"), func_num, nx);
	fpt = fopen(FileName, "r");
	if (fpt == NULL)
	{
		printf("\nError: An error occurred while attempting open file %s!\n",FileName);
		return;
	}
	if (func_num < 8)
	{
		M = (double *)malloc(nx * nx * sizeof(double));
		if (M == NULL)
		{
			return;
		}
		for (i = 0; i < nx * nx; i++)
		{
			fscanf(fpt, "%lf", &M[i]);
		}
	}
	else
	{
		M = (double *)malloc(f_num * nx * nx * sizeof(double));
		if (M == NULL)
		{
			return;
		}
		for (i = 0; i < f_num * nx * nx; i++)
		{
			fscanf(fpt, "%lf", &M[i]);
		}
	}
	fclose(fpt);


	// related to the init flag
	// n_flag=nx;
	// func_flag=func_num;
	// ini_flag=1;

	// }

	// switch for functions
	for (i = 0; i < mx; i++)
	{
		switch (func_num)
		{
		case 1:
			ellips_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=1000.0;
			break;
		case 2:
			bent_cigar_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=1000.0;
			break;
		case 3:
			rosenbrock_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=200.0;
			break;
		case 4:
			rastrigin_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=200.0;
			break;

		case 5:
			schwefel_func(&x[i * nx], &f[i], nx, 1,O, 1,M);
			f[i]+=300;
			break;
		case 6:
			alple_func(&x[i * nx], &f[i], nx, 1, O, 1,M);
			f[i]+=300;
			break;
		case 7:
			happycat_func(&x[i * nx], &f[i], nx, 1,O, 0,M);
			f[i]+=500;
			break;
		case 8:
			cf01(&x[i * nx], &f[i], nx, 1,O, 0,M);
			break;
		case 9:
			cf02(&x[i * nx], &f[i], nx, 1,O, 0,M);
			break;
		case 10:
			cf03(&x[i * nx], &f[i], nx, 1,O, 0,M);
			break;

		default:
			f[i] = 0.0;
			break;
		}
	}
}

void bent_cigar_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Bent_Cigar */
{
	double index;
	int i;

	if(nx==50)
		index=4.0;
	else
		index=6.0;
	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag);
	f[0] = z[0] * z[0];
	for (i = 1; i < nx; i++)
	{
		f[0] += pow(10.0, index) * z[i] * z[i];
	}
}

void weierstrass_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Weierstrass's  */
{
	int i, j, k_max;
	double sum, sum2, a, b;
	a = 0.5;
	b = 3.0;
	k_max = 20;
	f[0] = 0.0;

	shift_and_rotate(x, z, nx, Os, Mr, 0.5 / 100.0, s_flag, r_flag); /* shift and rotate */

	for (i = 0; i < nx; i++)
	{
		sum = 0.0;
		sum2 = 0.0;
		for (j = 0; j <= k_max; j++)
		{
			sum += pow(a, j) * cos(2.0 * PI * pow(b, j) * (z[i] + 0.5));
			sum2 += pow(a, j) * cos(2.0 * PI * pow(b, j) * 0.5);
		}
		f[0] += sum;
	}
	f[0] -= nx * sum2;
}

void ellips_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Ellipsoidal */
{
	int i;
	f[0] = 0.0;
	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
	for (i = 0; i < nx; i++)
	{
		f[0] += pow(10.0, 6.0 * i / (nx - 1)) * z[i] * z[i];
	}
}


void rosenbrock_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Rosenbrock's */
{
	int i;
	double tmp1, tmp2;
	f[0] = 0.0;
	shift_and_rotate(x, z, nx, Os, Mr, 2.048 / 100.0, s_flag, r_flag); /* shift and rotate */
	z[0] += 1.0;											  // shift to orgin
	for (i = 0; i < nx - 1; i++)
	{
		z[i + 1] += 1.0; // shift to orgin
		tmp1 = z[i] * z[i] - z[i + 1];
		tmp2 = z[i] - 1.0;
		f[0] += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	}
}

void schwefel_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Schwefel's  */
{
	int i;
	double tmp;
	f[0] = 0.0;

	shift_and_rotate(x, z, nx, Os, Mr, 1000.0 / 100.0, s_flag, r_flag); /* shift and rotate */

	for (i = 0; i < nx; i++)
	{
		z[i] += 4.209687462275036e+002;
		if (z[i] > 500)
		{
			f[0] -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
			tmp = (z[i] - 500.0) / 100;
			f[0] += tmp * tmp / nx;
		}
		else if (z[i] < -500)
		{
			f[0] -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
			tmp = (z[i] + 500.0) / 100;
			f[0] += tmp * tmp / nx;
		}
		else
			f[0] -= z[i] * sin(pow(fabs(z[i]), 0.5));
	}
	f[0] += 4.189828872724338e+002 * nx;
}

void rastrigin_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Rastrigin's  */
{
	int i;
	f[0] = 0.0;

	shift_and_rotate(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag, r_flag); /* shift and rotate */

	for (i = 0; i < nx; i++)
	{
		f[0] += (z[i] * z[i] - 5 * cos(2.0 * PI * z[i]) + 5.0);
	}
}

void shift(double *x, double *xshift, int nx, double *Os)
{
	int i;
	for (i = 0; i < nx; i++)
	{
		xshift[i] = x[i] - Os[i];
	}
}

void rotate(double *x, double *xrot, int nx, double *Mr)
{
	int i, j;
	for (i = 0; i < nx; i++)
	{
		xrot[i] = 0;
		for (j = 0; j < nx; j++)
		{
			xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
		}
	}
}

void shift_and_rotate(double *x, double *sr_x, int nx, double *Os, double *Mr, double sh_rate, int s_flag, int r_flag) /* shift and rotate */
{
	int i;
	if (s_flag == 1)
	{
		if (r_flag == 1)
		{
			shift(x, y, nx, Os);
			for (i = 0; i < nx; i++) // shrink to the orginal search range
			{
				y[i] = y[i] * sh_rate;
			}
			rotate(y, sr_x, nx, Mr);
		}
		else
		{
			shift(x, sr_x, nx, Os);
			for (i = 0; i < nx; i++) // shrink to the orginal search range
			{
				sr_x[i] = sr_x[i] * sh_rate;
			}
		}
	}
	else
	{

		if (r_flag == 1)
		{
			for (i = 0; i < nx; i++) // shrink to the orginal search range
			{
				y[i] = x[i] * sh_rate;
			}
			rotate(y, sr_x, nx, Mr);
		}
		else
			for (i = 0; i < nx; i++) // shrink to the orginal search range
			{
				sr_x[i] = x[i] * sh_rate;
			}
	}
}


void grie_rosen_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Griewank-Rosenbrock  */
{
	int i;
	double temp, tmp1, tmp2;
	f[0] = 0.0;

	shift_and_rotate(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

	z[0] += 1.0; // shift to orgin
	for (i = 0; i < nx - 1; i++)
	{
		z[i + 1] += 1.0; // shift to orgin
		tmp1 = z[i] * z[i] - z[i + 1];
		tmp2 = z[i] - 1.0;
		temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
		f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	}
	tmp1 = z[nx - 1] * z[nx - 1] - z[0];
	tmp2 = z[nx - 1] - 1.0;
	temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
	;
	f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
}

void discus_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Discus */
{
	int i;
	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
	f[0] = pow(10.0, 6.0) * z[0] * z[0];
	for (i = 1; i < nx; i++)
	{
		f[0] += z[i] * z[i];
	}
}


void ackley_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Ackley's  */
{
	int i;
	double sum1, sum2;
	sum1 = 0.0;
	sum2 = 0.0;

	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i = 0; i < nx; i++)
	{
		sum1 += z[i] * z[i];
		sum2 += cos(2.0 * PI * z[i]);
	}
	sum1 = -0.2 * sqrt(sum1 / nx);
	sum2 /= nx;
	f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
}
void happycat_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
/* original global optimum: [-1,-1,...,-1] */
{
	int i;
	double alpha, r2, sum_z;
	alpha = 1.0 / 8.0;

	shift_and_rotate(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

	r2 = 0.0;
	sum_z = 0.0;
	for (i = 0; i < nx; i++)
	{
		z[i] = z[i] - 1.0; // shift to orgin
		r2 += z[i] * z[i];
		sum_z += z[i];
	}
	f[0] = pow(fabs(r2 - nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
}

void katsuura_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Katsuura  */
{
	int i, j;
	double temp, tmp1, tmp2, tmp3;
	f[0] = 1.0;
	tmp3 = pow(1.0 * nx, 1.2);

	shift_and_rotate(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

	for (i = 0; i < nx; i++)
	{
		temp = 0.0;
		for (j = 1; j <= 32; j++)
		{
			tmp1 = pow(2.0, j);
			tmp2 = tmp1 * z[i];
			temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
		}
		f[0] *= pow(1.0 + (i + 1) * temp, 10.0 / tmp3);
	}
	tmp1 = 10.0 / nx / nx;
	f[0] = f[0] * tmp1 - tmp1;
}

void griewank_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Griewank's  */
{
	int i;
	double s, p;
	s = 0.0;
	p = 1.0;

	shift_and_rotate(x, z, nx, Os, Mr, 600.0 / 100.0, s_flag, r_flag); /* shift and rotate */

	for (i = 0; i < nx; i++)
	{
		s += z[i] * z[i];
		p *= cos(z[i] / sqrt(1.0 + i));
	}
	f[0] = 1.0 + s / 4000.0 - p;
}

void escaffer6_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;
	double temp1, temp2;

	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	f[0] = 0.0;
	for (i = 0; i < nx - 1; i++)
	{
		temp1 = sin(sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
		temp1 = temp1 * temp1;
		temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);
		f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
	}
	temp1 = sin(sqrt(z[nx - 1] * z[nx - 1] + z[0] * z[0]));
	temp1 = temp1 * temp1;
	temp2 = 1.0 + 0.001 * (z[nx - 1] * z[nx - 1] + z[0] * z[0]);
	f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
}

void alple_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;

	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	f[0] = 0.0;
	for (i = 0; i < nx ; i++)
	{
		f[0] += fabs(z[i]*sin(z[i])+0.1*z[i]);
	}
}

void dixon_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;

	shift_and_rotate(x, z, nx, Os, Mr, 10.0/100.0, s_flag, r_flag); /* shift and rotate */

	f[0] = (z[0]-1)*(z[0]-1);
	for (i = 1; i < nx ; i++)
	{
		f[0] += (i+1)*(2*z[i]*z[i]-z[i-1])*(2*z[i]*z[i]-z[i-1]);
	}
}

void ex3_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;
	double tmp;

	shift_and_rotate(x, z, nx, Os, Mr, 10.0/100.0, s_flag, r_flag); /* shift and rotate */

	f[0]=0.0;
	for (i = 0; i < nx ; i++)
	{
		if(z[i]==0)
			tmp=0;
		else
			tmp=pow(z[i],6.0)*(sin(1/z[i])+2);
		
		f[0] += tmp;
	}
}

void logexp_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;
	double tmp;

	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	f[0]=0.0;
	tmp=0.0;
	for (i = 0; i < nx ; i++)
	{
		tmp+=z[i]*z[i];
		f[0] += tmp;
	}
	f[0]=log(2*exp(-1e-4*tmp));
}


void invert_cos_wave_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;
	double tmp;

	shift_and_rotate(x, z, nx, Os, Mr, 5.0/100.0 , s_flag, r_flag); /* shift and rotate */

	f[0]=0.0;
	// for (i = 0; i < nx ; i++)
	// {
	// 	z[i]=z[i]/20.0;
	// }
	for (i = 0; i < nx-1 ; i++)
	{
		tmp+=z[i]*z[i]+z[i+1]*z[i+1]+0.5*z[i]*z[i+1];
		f[0] += -(exp(-tmp/8)*cos(4*sqrt(tmp)));
	}
	f[0]+=9;
}

void patho_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;
	double tmp1,tmp2;

	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); 

	f[0]=0.0;

	for (i = 0; i < nx-1 ; i++)
	{
		tmp1=sqrt(100*z[i]*z[i]+z[i+1]*z[i+1]);
		tmp2=1+1e-3*(z[i]*z[i]-2*z[i]*z[i+1]+z[i+1]*z[i+1]);
		f[0] += 0.5+((pow(sin(tmp1),2)-0.5)/tmp2);
	}
	f[0]*=1e3;
}
void salomon_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i;
	double tmp;

	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); 

	f[0]=0.0;
	tmp=0.0;
	for (i = 0; i < nx ; i++)
	{
		tmp=z[i]*z[i];
	}
	tmp=sqrt(tmp);
	f[0] = 1-cos(2*PI*tmp)+0.1*tmp;
	// f[0]*=1e3;
}
void sargan_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i,j;
	double *tmp;
	tmp = (double *)malloc(nx * sizeof(double));
	double sum=0.0;

	shift_and_rotate(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); 

	for(i=0;i<nx;i++)
	{
		sum+=z[i];
	}
	
	for(i=0;i<nx;i++)
	{
		tmp[i]=sum-z[i];
	}


	f[0]=0.0;
	for (i = 0; i < nx ; i++)
	{
		f[0]+=nx*(z[i]*z[i]+0.4*z[i]*tmp[i]);
	}
	// f[0]*=1e3;
}
void wavy_func(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr)
{
	int i,j;
	double *tmp;
	double k=10.0;
	tmp = (double *)malloc(nx * sizeof(double));
	double sum=0.0;

	shift_and_rotate(x, z, nx, Os, Mr, PI/100.0, s_flag, r_flag); 
	f[0]=0.0;
	// for(i=0;i<nx;i++)
	// {
	// 	z[i]=z[i]*PI/100;
	// }
	for(i=0;i<nx;i++)
	{
		f[0]+=1-cos(k*z[i])*exp(-0.5*z[i]*z[i]);
	}
	f[0]*=100.0;
	
}


// for caculate the composition function

void cf_cal2(double *x, double *f, double *bias, double *fit, int f_num)
{
	int i, j;
	double *w;
	double w_max = 0, w_sum = 0;
	double sigma = 2;
	w = (double *)malloc(f_num * sizeof(double));
	f[0] = 1;
	for (i = 0; i < f_num; i++)
	{
		fit[i] += bias[i];
		f[0] = f[0] * fit[i];
	}

	f[0] = 1e-3 * (exp(sigma * log(f[0] + 1)) - 1);
	free(w);
}

void cf01(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Composition Function 1 */
{
	int i, f_num = 4;
	double fit[4];
	double bias[4] = {0, 1e3, 1e2, 1e2};

	i = 0;
	ackley_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e4;
	i = 1;
	ellips_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e+10;
	i = 2;
	griewank_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e4;
	i = 3;
	rastrigin_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 3e1;
	cf_cal2(x, f, bias, fit, f_num);
}

void cf02(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Composition Function 2 */
{
	int i, f_num = 4;
	double fit[4];
	double bias[4] = {0, 1e1, 1e2, 1e3};
	i = 0;
	alple_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e3;
	i = 1;
	katsuura_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = 1*fit[i] / 5e1;
	i = 2;
	rosenbrock_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e1;
	i = 3;
	rastrigin_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e3;
	cf_cal2(x, f, bias, fit, f_num);
}

void cf03(double *x, double *f, int nx, int s_flag, double *Os, int r_flag, double *Mr) /* Composition Function 3 */
{
	int i, f_num = 5;
	double fit[5];
	double bias[5] = {0, 1e2, 1e2, 1e2, 1e2};
	i = 0;
	happycat_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e+5;
	i = 1;
	grie_rosen_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e+5;
	i = 2;
	schwefel_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e+5;
	i = 3;
	weierstrass_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e+8;
	i = 4;
	ellips_func(x, &fit[i], nx, s_flag ,&Os[i * nx],r_flag ,&Mr[i * nx * nx]);
	fit[i] = fit[i] / 1e+10;
	cf_cal2(x, f, bias, fit, f_num);
}

