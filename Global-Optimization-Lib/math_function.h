
#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include <memory>


class MathFunction
{
public:
	typedef std::shared_ptr<MathFunction> MathFunctionPtr;

	virtual double getValue(double arg) = 0;
};

#endif
