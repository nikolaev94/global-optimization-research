
#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H

#include <memory>
#include <vector>

#include "math_function.h"

class OptProblem
{
public:
	typedef std::shared_ptr<OptProblem> OptProblemP;

	struct ContourData
	{
		unsigned int grid_size;
		std::vector<double> x_values;
		std::vector<double> y_values;
		std::vector<std::vector<double>> z_values;

		ContourData(unsigned int in_grid_size) : grid_size(in_grid_size)
		{
			x_values.reserve(grid_size);
			y_values.reserve(grid_size);

			z_values.assign(grid_size, std::vector<double>());
		}

		void put_x(double in_x)
		{
			x_values.push_back(in_x);
		}

		void put_y(double in_y)
		{
			y_values.push_back(in_y);
		}


		void put_z(std::size_t row, double in_z)
		{
			z_values[row].push_back(in_z);
		}
	};

	virtual double getObjectiveValue(double arg) = 0;
	virtual double getContraintValue(std::size_t number, double arg) = 0;
	virtual std::size_t getConstraintsNumber() const = 0;
	virtual double getReferenceMinError(double scalar) = 0;
	virtual double getReferenceMinimum() = 0;
	virtual unsigned int getDimention() = 0;

	virtual void mapScalarToVector(double scalar, std::vector<double>& out_point) = 0;
	virtual void initContourData(ContourData&) = 0;
};

#endif
