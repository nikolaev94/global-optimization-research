
#ifndef GRISHAGIN_PROBLEM_H
#define GRISHAGIN_PROBLEM_H

#include "grishagin_function.h"
#include "opt_problem.h"

class GrishaginProblem : public OptProblem
{
private:
	static double rand_minimums[];
	static double getDistance(double[2], double[2]);

	GrishaginFunction::GrishaginFunctionPtr objective;

public:

	GrishaginProblem(GrishaginFunction::GrishaginFunctionPtr);
	~GrishaginProblem();
	double getReferenceMinError(double scalar) override;
	unsigned getDimention() override;

	double getObjectiveValue(double arg) const override;
	double getContraintValue(std::size_t number, double arg) const override;
	std::size_t getConstraintsNumber() const override;
};

#endif
