
#include "VagrisProblem.h"


VagrisProblem::VagrisProblem(VagrisFunctionP in_objective) : objective(in_objective)
{
	generator.SetObjective(in_objective.get());

	problem = generator.GenerateProblem(CONSTRAINTED_PROBLEM_GENERATOR_FLAGS);
}


VagrisProblem::VagrisProblem(VagrisFunctionP in_objective,
	const std::vector<VagrisFunctionP>& in_constrains) : objective(in_objective)
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


void VagrisProblem::linearTransform(double point[])
{
	for (int i = 0; i < problem.GetDimension(); i++)
	{
		point[i] = point[i] + 0.5;
	}
}


void VagrisProblem::reverseLinearTransform(double point[])
{
	for (int i = 0; i < problem.GetDimension(); i++)
	{
		point[i] = point[i] - 0.5;
	}
}


void VagrisProblem::mapScalarToNDimSpace(double scalar, double point[])
{
	mapd(scalar, PRECISION, point, problem.GetDimension(), KEY);

	linearTransform(point);
}


void VagrisProblem::mapNDimVectorToScalar(double point[], double* scalar)
{
	reverseLinearTransform(point);

	xyd(scalar, PRECISION, point, problem.GetDimension());
}


double VagrisProblem::getObjectiveValue(double scalar)
{
	double* point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, point);

	double value = problem.CalculateFunction(point, 0);

	delete[] point;

	return value;
}


double VagrisProblem::getContraintValue(std::size_t number, double scalar)
{
	double* point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, point);

	double value = problem.CalculateFunction(point, number);

	delete[] point;

	return value;
}


std::size_t VagrisProblem::getConstraintsNumber() const
{
	return problem.GetConstraintsNumber();
}


double VagrisProblem::getEuclideanDistance(double lhs[], double rhs[])
{
	double sum = 0.0;

	for (int i = 0; i < problem.GetDimension(); i++)
	{
		sum += std::pow(lhs[i] - rhs[i], 2.0);
	}

	return sqrt(sum);
}


double VagrisProblem::getReferenceMinError(double scalar)
{
	double* trial_point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, trial_point);

	double* min_point = new double[problem.GetDimension()];

	problem.GetOptimumPoint(min_point);

	double min_point_scalar = 0.0;

	double distance = getEuclideanDistance(trial_point, min_point);

	delete[] min_point;

	delete[] trial_point;

	return distance;
}


double VagrisProblem::getReferenceMinimum()
{
	double* min_point = new double[problem.GetDimension()];

	problem.GetOptimumPoint(min_point);

	double min_scalar = 0.0;

	mapNDimVectorToScalar(min_point, &min_scalar);

	delete[] min_point;

	return min_scalar;
}


unsigned int VagrisProblem::getDimention()
{
	return problem.GetDimension();
}


void VagrisProblem::mapScalarToVector(double scalar, std::vector<double>& out_point)
{
	double* point = new double[problem.GetDimension()];

	mapScalarToNDimSpace(scalar, point);

	for (int i = 0; i < problem.GetDimension(); i++)
	{
		out_point.push_back(point[i]);
	}

	delete[] point;
}


void VagrisProblem::initContourData(ContourData& contour_data)
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


	std::size_t row = 0;
	for (const auto x_comp : contour_data.x_values)
	{
		for (const auto y_comp : contour_data.y_values)
		{
			double node_point[] = { x_comp, y_comp };

			contour_data.put_z(row, problem.CalculateFunction(node_point, 0));
		}
		++row;
	}

	delete[] leftBound;
	delete[] rightBound;
}

