
/*
#ifndef HILL_H
#define HILL_H

#include "all.h"
#include "function.h"

class Hill :
	public Function {
private:
	double** a, **b;

public:
	Hill() :a(nullptr), b(nullptr) {}
	Hill(int funcsNum) :Function(funcsNum) {
		a = new double*[funcsNum];
		b = new double*[funcsNum];
		for (int i = 0; i < funcsNum; i++) {
			a[i] = new double[14];
			b[i] = new double[14];
		}
	}

	void setData(double* buf, double delta) override {
		if (!a || !b)
			return;
		for (int i = 0; i < funcsNum; i++) {
			std::copy(&buf[i * 28],
				&buf[i * 28 + 14], 
				stdext::checked_array_iterator<double*>(a[i], 14));
			std::copy(&buf[i * 28 + 14],
				&buf[(i + 1) * 28],
				stdext::checked_array_iterator<double*>(b[i], 14));
		}
		int iend = static_cast<int>(1.0 / STEP) + 1;
		double min = DBL_MAX;
		for (int i = 0; i < iend; i++) {
			double x = i * STEP, max = -DBL_MAX;
			for (int j = 1; j < funcsNum; j++) {
				double tmp = calculateConstraintValue(x, j);
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

	~Hill() {
		for (int i = 0; i < funcsNum; i++) {
			delete[] a[i];
			delete[] b[i];
		}
		delete[] a;
		delete[] b;
		a = nullptr;
		b = nullptr;
	}

	double calculateObjectiveValue(double x) const override {
		double f = 0;
		for (int i = 0; i < 14; i++) {
			f += a[0][i] * sin(2 * M_PI * (i + 1) * x)
				+ b[0][i] * cos(2 * M_PI * (i + 1) * x);
		}
		return f;
	}

	double calculateConstraintValue(double x, int j) const override {
		double f = 0;
		for (int i = 0; i < 14; i++) {
			f += a[j][i] * sin(2 * M_PI * (i + 1) * x)
				+ b[j][i] * cos(2 * M_PI *(i + 1) * x);
		}
		return f - minimax - delta;
	}
};

#endif
*/
