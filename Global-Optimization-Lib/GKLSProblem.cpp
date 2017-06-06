
#include "GKLSProblem.h"


GKLSProblem::GKLSProblem(GKLSFunctionP in_objective) : objective(in_objective)
{
	generator.SetObjective(in_objective.get());

	problem = generator.GenerateProblem(DEFAULT_GENERATOR_FLAGS);
}


GKLSProblem::GKLSProblem(GKLSFunctionP in_objective,
	const std::vector<GKLSFunctionP>& in_constrains) : objective(in_objective)
{
	generator.SetObjective(in_objective.get());

	for (const auto& function : in_constrains)
	{
		constrains.push_back(function);

		generator.AddConstraint(function.get(), CONSTRAINT_SHIFT_PARAMETER, 1, 5.0);
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
	/*double* left_bounds = new double[problem.GetDimension()];

	double* right_bounds = new double[problem.GetDimension()];*/


	for (int i = 0; i < problem.GetDimension(); i++)
	{
		point[i] = point[i] * 2.0;
	}


	/*problem.GetBounds(left_bounds, right_bounds);

	for (unsigned int i = 0; i < problem.GetDimension(); i++)
	{
		point[i] = point[i] * (right_bounds[i] - left_bounds[i])
			+ (right_bounds[i] + left_bounds[i]) / 2;
	}*/

	/*delete[] left_bounds;

	delete[] right_bounds;*/
}


void GKLSProblem::reverseLinearTransform(double point[])
{

	/*double* left_bounds = new double[problem.GetDimension()];

	double* right_bounds = new double[problem.GetDimension()];*/


	for (int i = 0; i < problem.GetDimension(); i++)
	{
		point[i] = point[i] * 0.5;
	}

	/*delete[] left_bounds;
	delete[] right_bounds;*/

	/*double* left_bounds = new double[problem.GetDimension()];

	double* right_bounds = new double[problem.GetDimension()];

	problem.GetBounds(left_bounds, right_bounds);

	for (unsigned int i = 0; i < problem.GetDimension(); i++)
	{
		point[i] = (point[i] - (left_bounds[i] + right_bounds[i]) / 2.0)
			/ (right_bounds[i] - left_bounds[i]);
	}

	delete[] left_bounds;
	delete[] right_bounds;*/
}


void GKLSProblem::mapScalarToNDimSpace(double x, double point[])
{
	mapd(x, PRECISION, point, problem.GetDimension(), KEY);

	linearTransform(point);
}


void GKLSProblem::mapNDimVectorToScalar(double point[], double *x)
{
	reverseLinearTransform(point);

	xyd(x, PRECISION, point, problem.GetDimension());
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


double GKLSProblem::getDistance(double *lhs, double *rhs)
{
	double sum = 0.0;

	for (int i = 0; i < problem.GetDimension(); i++)
	{
		sum += std::pow(lhs[i] - rhs[i], 2.0);
	}

	return sqrt(sum);
}


double GKLSProblem::getReferenceMinError(double scalar)
{
	double* min_point = new double[problem.GetDimension()];

	problem.GetOptimumPoint(min_point);

	double* trial_point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, trial_point);

	double distance = getDistance(trial_point, min_point);

	delete[] trial_point;

	delete[] min_point;

	return distance;
}


double GKLSProblem::getReferenceMinimum()
{
	double* min_point = new double[problem.GetDimension()];

	problem.GetOptimumPoint(min_point);

	double min_scalar = 0.0;

	mapNDimVectorToScalar(min_point, &min_scalar);

	delete[] min_point;

	return min_scalar;
}


unsigned int GKLSProblem::getDimention()
{
	return problem.GetDimension();
}

void GKLSProblem::mapScalarToVector(double scalar, std::vector<double>& out_point)
{
	double* point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, point);

	for (int i = 0; i < problem.GetDimension(); i++)
	{
		out_point.push_back(point[i]);
	}

	delete[] point;
}


void GKLSProblem::initContourData(ContourData& contour_data)
{
	double* leftBound = new double[problem.GetDimension()];
	double* rightBound = new double[problem.GetDimension()];

	problem.GetBounds(leftBound, rightBound);

	double step_x = abs(rightBound[0] - leftBound[0]) / contour_data.grid_size;

	double step_y = abs(rightBound[1] - leftBound[1]) / contour_data.grid_size;

	for (unsigned int i = 0; i < contour_data.grid_size; i++)
	{
		double yi = leftBound[1] + i * step_y;

		contour_data.put_y(yi);
	}

	for (unsigned int j = 0; j < contour_data.grid_size; j++)
	{
		double xj = leftBound[0] + j * step_y;

		contour_data.put_x(xj);
	}


	for (const auto x_comp : contour_data.x_values)
	{
		std::size_t row = 0;
		for (const auto y_comp : contour_data.y_values)
		{
			double node_point[] = { x_comp, y_comp };

			contour_data.put_z(row++, problem.CalculateFunction(node_point, 0));
		}
	}

	delete[] leftBound;
	delete[] rightBound;
}
