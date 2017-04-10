
#include "grishagin_function.h"

unsigned char GrishaginFunction::matcon[10][45] = {
	{ 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0 },
	{ 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1 },
	{ 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0 },
	{ 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 },
	{ 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1 },
	{ 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0 },
	{ 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
	{ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1 },
	{ 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1 },
	{ 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0 },
};

void GrishaginFunction::mapScalar_To_2DSpace(double x, double point[DIMENSION])
{
	mapd(x, PRECISION, point, DIMENSION, KEY);
	linear_transform(point);
}

void GrishaginFunction::linear_transform(double point[DIMENSION]) {
	for (int i = 0; i < DIMENSION; i++)
	{
		point[i] = point[i] + 0.5;
	}
}

GrishaginFunction::GrishaginFunction(int _sequential_number) : MathFunction()
{
	this->sequential_number = _sequential_number;
	set_random(sequential_number);
}

GrishaginFunction::~GrishaginFunction() {}

void GrishaginFunction::gen(unsigned char k[], unsigned char k1[], int kap1, int kap2)
{
	int jct, i, j;

	jct = 0;
	for (i = kap2; i >= kap1; i--){
		j = (k[i] + k1[i] + jct) / 2;
		k[i] = k[i] + k1[i] + jct - j * 2;
		jct = j;
	}
	if (jct != 0)
		for (i = kap2; i >= kap1; i--){
			j = (k[i] + jct) / 2;
			k[i] = k[i] + jct - j * 2;
			jct = j;
		}
}

double GrishaginFunction::rndm20(unsigned char k[])
{
	int i;
	unsigned char k1[45];
	double de2, rndm;

	for (i = 0; i<38; i++)
		k1[i] = k[i + 7];
	for (i = 38; i<45; i++)
		k1[i] = 0;
	for (i = 0; i<45; i++)
		k[i] = abs(k[i] - k1[i]);
	for (i = 27; i<45; i++)
		k1[i] = k[i - 27];
	for (i = 0; i<27; i++)
		k1[i] = 0;

	gen(k, k1, 9, 44);
	gen(k, k1, 0, 8);

	rndm = 0.;
	de2 = 1.;
	for (i = 0; i<36; i++){
		de2 = de2 / 2;
		rndm = rndm + k[i + 9] * de2;
	}
	return (rndm);
}

void GrishaginFunction::set_random(int nf)
{
	int lst, i, j, i1, i2, i3;

	if (nf<1 || nf>100)
		nf = 1;
	lst = 10;
	i1 = (nf - 1) / lst;
	i2 = i1*lst;
	for (j = 0; j<45; j++)
		icnf[j] = matcon[i1][j];
	if (i2 != (nf - 1)){
		i3 = nf - 1 - i2;
		for (j = 1; j <= i3; j++)
			for (i = 0; i<196; i++)
				rndm20(icnf);
	}
	for (j = 0; j<7; j++)
		for (i = 0; i<7; i++){
			af[i][j] = 2.*rndm20(icnf) - 1.;
			cf[i][j] = 2.*rndm20(icnf) - 1.;
		}
	for (j = 0; j<7; j++)
		for (i = 0; i<7; i++){
			bf[i][j] = 2.*rndm20(icnf) - 1.;
			df[i][j] = 2.*rndm20(icnf) - 1.;
		}
}

double GrishaginFunction::random_func(double* y, int n)
{
	int i, j;
	double d1, d2, sx1, cx1, sy1, cy1;

	d1 = M_PI*y[0];
	d2 = M_PI*y[1];
	sx1 = sin(d1);
	cx1 = cos(d1);
	sy1 = sin(d2);
	cy1 = cos(d2);
	snx[0] = sx1;
	snx[0] = sx1;
	csx[0] = cx1;
	sny[0] = sy1;
	csy[0] = cy1;
	for (i = 0; i<6; i++){
		snx[i + 1] = snx[i] * cx1 + csx[i] * sx1;
		csx[i + 1] = csx[i] * cx1 - snx[i] * sx1;
		sny[i + 1] = sny[i] * cy1 + csy[i] * sy1;
		csy[i + 1] = csy[i] * cy1 - sny[i] * sy1;
	}
	d1 = 0;
	d2 = 0;
	for (i = 0; i<7; i++)
		for (j = 0; j<7; j++){
			d1 = d1 + af[i][j] * snx[i] * sny[j] + bf[i][j] * csx[i] * csy[j];
			d2 = d2 + cf[i][j] * snx[i] * sny[j] - df[i][j] * csx[i] * csy[j];
		}
	return -(sqrt(d1*d1 + d2*d2));
}

double GrishaginFunction::getValue(double arg)
{
	double point[2];
	mapd(arg, PRECISION, point, DIMENSION, KEY);
	linear_transform(point);
	return random_func(point, DIMENSION);
}

int GrishaginFunction::getFunctionSequentialNumber()
{
	return this->sequential_number;
}
