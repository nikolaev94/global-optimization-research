#ifndef USER_INPUT_H
#define USER_INPUT_H

#include <iostream>
#include <string>

#include <solver.h>

enum FunctionClass
{
	GKLS_SIMPLE, GKLS_HARD, VAGRIS
};

struct UserParameters
{
private:
	unsigned int dimensions;
	unsigned int num_workers;
	double method_parameter;
	double precision;
	FunctionClass function_class;
	Solver::SolvingMethod solving_method;

	void show_usage();

public:
	UserParameters() : dimensions(2), num_workers(2), method_parameter(5.0),
		precision(0.1), function_class(FunctionClass::GKLS_SIMPLE), solving_method(Solver::DYNAMIC) {}

	void parse_arguments_from_command_line(int argc, char* argv[]);
	unsigned int get_dimensions() const;
	FunctionClass get_function_class() const;
	unsigned int get_workers_number() const;
	double get_method_parameter() const;
	Solver::SolvingMethod get_solving_method() const;
	double get_precision() const;
};

#endif
