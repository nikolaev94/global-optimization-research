
#ifndef GRISHAGIN_FUNCTION_H
#define GRISHAGIN_FUNCTION_H

#define _USE_MATH_DEFINES

#include <cmath>

#include <memory>

#include "map/map.h"
#include "math_function.h"

class GrishaginFunction : public MathFunction
{
private:
	static const int DIMENSION = 2, KEY = 1, PRECISION = 10;

	static unsigned char matcon[10][45];
	static double rand_minimums[];

	static void linear_transform(double point[DIMENSION]);

	unsigned char icnf[45];
	double af[7][7], bf[7][7], cf[7][7], df[7][7];
	double snx[7], csx[7], sny[7], csy[7];

	int sequential_number = 0;

	void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);
	double rndm20(unsigned char k[]);

	void set_random(int nf);
	double random_func(double* y, int n);

public:

	typedef std::shared_ptr<GrishaginFunction> GrishaginFunctionPtr;

	GrishaginFunction(int _sequential_number = 0);
	~GrishaginFunction();

	static void mapScalar_To_2DSpace(double x, double point[DIMENSION]);

	int getFunctionSequentialNumber();
	double getValue(double arg) override;
};

#endif
