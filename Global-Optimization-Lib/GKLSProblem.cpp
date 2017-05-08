#include "GKLSProblem.h"

GKLSProblem::~GKLSProblem()
{
}


GKLSProblem::GKLSProblem(GKLSFunctionPtr in_objective) : objective(in_objective)
{
	generator.SetObjective(in_objective.get());

	problem = generator.GenerateProblem(0);
}


GKLSProblem::GKLSProblem(GKLSFunctionPtr in_objective,
	const std::vector<GKLSFunctionPtr>& in_constrains) : objective(in_objective)
{
	generator.SetObjective(in_objective.get());

	for (const auto& function : in_constrains)
	{
		constrains.push_back(function);

		generator.AddConstraint(function.get(), CONSTRAINT_SHIFT_PARAMETER);
	}

	if (in_constrains.empty())
	{
		problem = generator.GenerateProblem(DEFAULT_GENERATOR_FLAGS);
	}
	else
	{
		problem = generator.GenerateProblem(CONSTRAINTED_PROBLEM_GENERATOR_FLAGS);
	}
}


void GKLSProblem::linearTransform(double point[])
{
	double* left_bounds = new double[problem.GetDimension()];

	double* right_bounds = new double[problem.GetDimension()];

	problem.GetBounds(left_bounds, right_bounds);

	for (unsigned int i = 0; i < problem.GetDimension() ; i++)
	{
		point[i] = point[i] * (right_bounds[i] - left_bounds[i])
			+ (right_bounds[i] + left_bounds[i]) / 2;
	}

	delete[] left_bounds;
	delete[] right_bounds;
}


void GKLSProblem::mapScalarToNDimSpace(double x, double point[])
{
	mapd(x, PRECISION, point, problem.GetDimension(), KEY);

	linearTransform(point);
}


double GKLSProblem::getObjectiveValue(double scalar)
{
	double* point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, point);

	double value = problem.CalculateFunction(point, 0);

	delete[] point;

	return value;
}


double GKLSProblem::getContraintValue(std::size_t number, double scalar)
{
	double* point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, point);

	double value = problem.CalculateFunction(point, number);

	delete[] point;

	return value;
}


std::size_t GKLSProblem::getConstraintsNumber() const
{
	if (problem.GetConstraintsNumber() < 0)
	{
		return 0;
	}
	return problem.GetConstraintsNumber();
}


double GKLSProblem::getEuclideanDistance(double *lhs, double *rhs)
{
	unsigned int i;
	double norm = 0.0;
	for (i = 0; i < problem.GetDimension(); i++)
		norm += (lhs[i] - rhs[i])*(lhs[i] - rhs[i]);
	return sqrt(norm);
}


double GKLSProblem::getReferenceMinError(double scalar)
{
	double* trial_point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, trial_point);

	double* min_point = new double[problem.GetDimension()];

	problem.GetOptimumPoint(min_point);

	double distance = getEuclideanDistance(trial_point, min_point);

	delete[] min_point;

	delete[] trial_point;

	return distance;
}


double GKLSProblem::getReferenceMinimum()
{
	double* min_point = new double[problem.GetDimension()];

	problem.GetOptimumPoint(min_point);

	double min_scalar = 0.0;

	xyd(&min_scalar, PRECISION, min_point, problem.GetDimension());

	return min_scalar;
}


unsigned int GKLSProblem::getDimention()
{
	return problem.GetDimension();
}

