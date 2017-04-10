
/*
#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "all.h"
#include "function.h"

class Integral
	: public Function {
private:
	const int N = 100;
	const double LEFT = -M_PI / 2, RIGTH = M_PI / 2;
	double** a, **b, *alpha;

	double getIntegrandFunctionValue(double x, double t, int j) const {
		double f = 0;
		for (int i = 0; i < 14; i++) {
			f += a[j][i] * sin(2 * M_PI * (i + 1) * (alpha[j] - x) * t)
				+ b[j][i] * cos(2 * M_PI * (i + 1) * x * (alpha[j] - x) * t);
		}
		return f;
	}

	double integrate(double x, int j) const {
		double h = (RIGTH - LEFT) / N;
		double sum1 = 0.0, sum2 = 0.0;
		for (int i = 0; i < N; i++) {
			sum1 += getIntegrandFunctionValue(x, LEFT + h * i + h / 2, j);
		}
		for (int i = 1; i < N; i++) {
			sum2 += getIntegrandFunctionValue(x, LEFT + h * i, j);
		}
		return h / 6.0 * (getIntegrandFunctionValue(x, LEFT, j)
			+ getIntegrandFunctionValue(x, RIGTH, j) + 4.0 * sum1 + 2.0 * sum2);
	}

public:

	Integral() : a(nullptr), b(nullptr), alpha(nullptr) {}

	Integral(int funcsNum) : Function(funcsNum) {
		a = new double*[funcsNum];
		b = new double*[funcsNum];
		alpha = new double[funcsNum];
		for (int i = 0; i < funcsNum; i++) {
			a[i] = new double[14];
			b[i] = new double[14];
		}
	}

	~Integral() {
		for (int i = 0; i < funcsNum; i++) {
			delete[] a[i];
			delete[] b[i]; 
		}
		delete[] a;
		delete[] b;
		delete[] alpha;
		a = nullptr;
		b = nullptr;
		alpha = nullptr;
	}

	double calculateObjectiveValue(double x) const override {
		return integrate(x, 0);
	}

	double calculateConstraintValue(double x, int j) const override {
		return integrate(x, j) - minimax - delta;
	}

	void setData(double* buf, double delta) override {
		if (!a || !b || !alpha)
			return;
		for (int i = 0; i < funcsNum; i++) {
			std::copy(&buf[i * 29], &buf[i * 29 + 14],
				stdext::checked_array_iterator<double*>(a[i], 14));
			std::copy(&buf[i * 29 + 14], &buf[i * 29 + 28],
				stdext::checked_array_iterator<double*>(b[i], 14));
			std::copy(&buf[i * 29 + 28], &buf[(i + 1) * 29], 
				stdext::checked_array_iterator<double*>(&alpha[i], 1));
		}
		if (funcsNum == 1) {
			return;
		}
		int iend = static_cast<int>(1.0 / STEP) + 1;
		double min = DBL_MAX;
		for (int i = 0; i < iend; i++) {
			double x = i * STEP, max = -DBL_MAX;
			for (int j = 1; j < funcsNum; j++) {
				double tmp = calculateConstraintValueFast(x, j);
				if (tmp > max) {
					max = tmp;
				}
			}
			if (max < min) {
				min = max;
			}
		}
		this->minimax = min;
		this->delta = delta;
	}

	double calculateObjectiveValueFast(double x) const override {
		double r = 0.0, l = 0.0;
		if (fabs(alpha[0] - x) < DBL_EPSILON) {
			for (int i = 0; i < 14; i++) {
				r += b[0][i] * RIGTH;
				l += b[0][i] * LEFT;
			}
		} else if (fabs(x) < DBL_EPSILON) {
			for (int i = 0; i < 14; i++) {
				r += -a[0][i] * cos(2 * M_PI * (i + 1) * (alpha[0] - x) * RIGTH)
					/ (2 * M_PI * (i + 1)  * (alpha[0] - x)) + b[0][i] * RIGTH;
				l += -a[0][i] * cos(2 * M_PI * (i + 1) * (alpha[0] - x) * LEFT)
					/ (2 * M_PI * (i + 1)  * (alpha[0] - x)) + b[0][i] * LEFT;
			}
		} else {
			for (int i = 0; i < 14; i++) {
				r += (-a[0][i] * cos(2 * M_PI * (i + 1) * (alpha[0] - x) * RIGTH) * x
					+ b[0][i] * sin(2 * M_PI * (i + 1) * x * (alpha[0] - x) * RIGTH))
					/ (2 * M_PI * (i + 1) * x * (alpha[0] - x));
				l += (-a[0][i] * cos(2 * M_PI * (i + 1) * (alpha[0] - x) * LEFT) * x
					+ b[0][i] * sin(2 * M_PI * (i + 1) * x * (alpha[0] - x) * LEFT))
					/ (2 * M_PI * (i + 1) * x * (alpha[0] - x));
			}
		}
		return r - l;
	}

	double calculateConstraintValueFast(double x, int j) const override {
		double r = 0.0, l = 0.0;
		if (fabs(alpha[0] - x) < DBL_EPSILON) {
			for (int i = 0; i < 14; i++) {
				r += b[j][i] * RIGTH;
				l += b[j][i] * LEFT;
			}
		} else if (fabs(x) < DBL_EPSILON) {
			for (int i = 0; i < 14; i++) {
				r += -a[j][i] * cos(2 * M_PI * (i + 1) * (alpha[j] - x) * RIGTH)
					/ (2 * M_PI * (i + 1)  * (alpha[j] - x)) + b[j][i] * RIGTH;
				l += -a[j][i] * cos(2 * M_PI * (i + 1) * (alpha[j] - x) * LEFT)
					/ (2 * M_PI * (i + 1)  * (alpha[j] - x)) + b[j][i] * LEFT;
			}
		} else {
			for (int i = 0; i < 14; i++) {
				r += (-a[j][i] * cos(2 * M_PI * (i + 1) * (alpha[j] - x) * RIGTH) * x
					+ b[j][i] * sin(2 * M_PI * (i + 1) * x * (alpha[j] - x) * RIGTH))
					/ (2 * M_PI * (i + 1) * x * (alpha[j] - x));
				l += (-a[j][i] * cos(2 * M_PI * (i + 1) * (alpha[j] - x) * LEFT) * x
					+ b[j][i] * sin(2 * M_PI * (i + 1) * x * (alpha[j] - x) * LEFT))
					/ (2 * M_PI * (i + 1) * x * (alpha[j] - x));
			}
		}
		return r - l - minimax - delta;
	}
};

#endif
*/
