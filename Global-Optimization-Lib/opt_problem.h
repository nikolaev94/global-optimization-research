
#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H

#include <memory>

#include "all.h"
#include "math_function.h"

class OptProblem
{
private:
	// MathFunction* objective = nullptr;
	// std::vector<MathFunction*> constrains;

	/*double solution;
	void estimateSolution();*/

protected:
	//MathFunction::MathFunctionPtr objective;
	//std::vector<MathFunction::MathFunctionPtr> constrains;
public:
	typedef std::shared_ptr<OptProblem> OptProblemPtr;

	// OptProblem(MathFunction::MathFunctionPtr);
	// virtual ~OptProblem();
	// OptProblem(const OptProblem& src);

	virtual double getObjectiveValue(double arg) const = 0;
	virtual double getContraintValue(std::size_t number, double arg) const = 0;
	virtual std::size_t getConstraintsNumber() const = 0;
	virtual double getReferenceMinError(double scalar) = 0;
	virtual unsigned getDimention() = 0;

	// virtual void addConstraintFunc(MathFunction::MathFunctionPtr constraint) = 0;
};

#endif
