/*
#ifndef FUNCTION_H
#define FUNCTION_H

class Function {
protected:
	const double STEP = 0.00001;
	int funcsNum;
	double delta, minimax;

public:
	Function() : funcsNum(0), minimax(0.0), delta(0.0) {}
	Function(int funcsNum) : funcsNum(funcsNum), minimax(0.0), delta(0.0) {}

	virtual ~Function() {}

	virtual double calculateObjectiveValue(double x) const = 0;
	virtual double calculateConstraintValue(double x, int constraintIdx) const = 0;

	virtual void setData(double* buf, double delta) {}
	
	virtual double calculateObjectiveValueFast(double x) const {
		return calculateObjectiveValue(x);
	}

	virtual double calculateConstraintValueFast(double x, int constraintIdx) const {
		return calculateConstraintValue(x, constraintIdx);
	}
	
	int getFuncsNum() const {
		return funcsNum;
	}
};

#endif
*/
