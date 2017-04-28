#include "gkls_problem.h"

GklsProblem::GklsProblem(GklsFunction::GklsFunctionPtr obj) : objective(obj) {}

GklsProblem::GklsProblem(const GklsProblem& src)
{
	objective = src.objective;
}

GklsProblem::~GklsProblem() {}

double GklsProblem::getObjectiveValue(double arg) const
{
	return objective->getValue(arg);
}

double GklsProblem::getReferenceMinError(double scalar)
{
	double* trial_point = new double [GklsFunction::GKLS_dim];
	GklsFunction::mapScalar_To_nDSpace(scalar, trial_point);

	double* min_point = new double[GklsFunction::GKLS_dim];

	objective->getFunctionMinimum(min_point);

	double distance = GklsFunction::GKLS_norm(trial_point, min_point);

	delete[] min_point;
	delete[] trial_point;

	return distance;
}

double GklsProblem::getReferenceMinimum()
{
	double* min_point = new double[GklsFunction::GKLS_dim];

	objective->getFunctionMinimum(min_point);

	delete[] min_point;

	return 0.0;
}



double GklsProblem::getContraintValue(std::size_t number, double arg) const
{
	return 0.0;
}

std::size_t GklsProblem::getConstraintsNumber() const
{
	return 0;
}

unsigned GklsProblem::getDimention()
{
	return GklsFunction::GKLS_dim;
}
