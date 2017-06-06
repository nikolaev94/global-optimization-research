
#ifndef VAGRIS_PROBLEM_H
#define VAGRIS_PROBLEM_H

#include <memory>
#include <vector>

#include "opt_problem.h"

#include "examin\ConstrainedProblem.hpp"
#include "examin\ConstrainedProblemGenerator.hpp"
#include "examin\problemPar.hpp"

#include "Grishagin\grishagin_function.hpp"

#include "map\map.h"

class VagrisProblem : public OptProblem
{
public:
	typedef std::shared_ptr<vagrish::GrishaginFunction> VagrisFunctionP;
private:
	static const int KEY = 1;
	static const int PRECISION = 14;

	VagrisFunctionP objective;

	std::vector<VagrisFunctionP> constrains;

	TConstrainedProblemGenerator<vagrish::GrishaginFunction> generator;

	TConstrainedProblem<vagrish::GrishaginFunction> problem;

	const int CONSTRAINTED_PROBLEM_GENERATOR_FLAGS =
		IMPROVE_OBJECTIVE | TOTAL_DELTA;

	const int DEFAULT_GENERATOR_FLAGS = 0;

	const double CONSTRAINT_SHIFT_PARAMETER = 0.1f;

	void mapScalarToNDimSpace(double scalar, double point[]);

	void mapNDimVectorToScalar(double point[], double* scalar);

	void linearTransform(double point[]);

	void reverseLinearTransform(double point[]);

	double getEuclideanDistance(double lhs[], double rhs[]);

public:

	VagrisProblem(VagrisFunctionP);

	VagrisProblem(VagrisFunctionP, const std::vector<VagrisFunctionP>&);

	virtual double getObjectiveValue(double scalar) override;

	virtual double getContraintValue(std::size_t number, double scalar) override;

	virtual std::size_t getConstraintsNumber() const override;

	virtual double getReferenceMinError(double scalar) override;

	virtual double getReferenceMinimum() override;

	virtual unsigned int getDimention() override;

	virtual void mapScalarToVector(double scalar, std::vector<double>& out_point) override;

	virtual void initContourData(ContourData&) override;
};

#endif
