
#ifndef GKLS_PROBLEM_H
#define GKLS_PROBLEM_H

#include <vector>

#include "opt_problem.h"

#include "map\map.h"

#include "examin\ConstrainedProblem.hpp"
#include "examin\ConstrainedProblemGenerator.hpp"
#include "examin\problemPar.hpp"

#include "gkls\gkls_function.hpp"

class GKLSProblem : public OptProblem
{
public:
	typedef std::shared_ptr<gkls::GKLSFunction> GKLSFunctionP;
private:
	static const int KEY = 1;
	static const int PRECISION = 10;

	const double CONSTRAINT_SHIFT_PARAMETER = 0.01f;

	const int CONSTRAINTED_PROBLEM_GENERATOR_FLAGS =
		IMPROVE_OBJECTIVE | TOTAL_DELTA; //| ZOOM | SHIFT;

	const int DEFAULT_GENERATOR_FLAGS = 0;

	GKLSFunctionP objective;

	std::vector<GKLSFunctionP> constrains;

	TConstrainedProblemGenerator<gkls::GKLSFunction> generator;

	TConstrainedProblem<gkls::GKLSFunction> problem;

	void linearTransform(double point[]);

	void mapScalarToNDimSpace(double scalar, double point[]);

	double getEuclideanDistance(double* lhs, double* rhs);

public:
	static const gkls::GKLSFuncionType DEFAULT_GKLS_FUNCTION_TYPE = gkls::GKLSFuncionType::TD2;

	GKLSProblem(GKLSFunctionP);

	GKLSProblem(GKLSFunctionP, const std::vector<GKLSFunctionP>&);

	virtual double getObjectiveValue(double scalar) override;

	virtual double getContraintValue(std::size_t number, double scalar) override;

	virtual std::size_t getConstraintsNumber() const override;

	virtual double getReferenceMinError(double scalar) override;

	virtual double getReferenceMinimum() override;

	virtual unsigned getDimention() override;
};

#endif
