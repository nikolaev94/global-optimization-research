
#ifndef GKLS_PROBLEM_H
#define GKLS_PROBLEM_H

#include "opt_problem.h"
#include "gkls_function.h"


class GklsProblem : public OptProblem
{
private:
	GklsFunction::GklsFunctionPtr objective;

public:

	GklsProblem(GklsFunction::GklsFunctionPtr);
	~GklsProblem();
	GklsProblem(const GklsProblem&);

	virtual double getObjectiveValue(double arg) const override;
	virtual double getContraintValue(std::size_t number, double arg) const override;
	virtual std::size_t getConstraintsNumber() const override;

	virtual double getReferenceMinError(double scalar) override;
	virtual double getReferenceMinimum() override;
	virtual unsigned getDimention() override;
};

#endif
