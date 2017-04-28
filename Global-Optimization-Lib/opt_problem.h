
#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H

#include <memory>

#include "math_function.h"

class OptProblem
{

public:
	typedef std::shared_ptr<OptProblem> OptProblemPtr;

	virtual double getObjectiveValue(double arg) const = 0;
	virtual double getContraintValue(std::size_t number, double arg) const = 0;
	virtual std::size_t getConstraintsNumber() const = 0;
	virtual double getReferenceMinError(double scalar) = 0;
	virtual double getReferenceMinimum() = 0;
	virtual unsigned getDimention() = 0;
};

#endif
